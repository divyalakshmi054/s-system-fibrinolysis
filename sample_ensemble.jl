# load the include -
include("Include.jl")
using Plots

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Fibrinolysis.bst")

# build the default model structure -
model = build(path_to_model_file)

# load the training data -
_PATH_TO_DATA = joinpath(pwd(),"data")
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Composition-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# which visit?
visit = 4

#let's filter visit 4s since we look to train using that visit
visit_df = filter(:Visit => x->(x==visit), training_df) 

# size of training set -
(R,C) = size(visit_df)

# main simulation -
SF = 1e9
for i ∈ 1:R

    # build new model -
    dd = deepcopy(model)

    # setup static -
    sfa = dd.static_factors_array
    sfa[1] = 0.0                    # 1 tPA     SET tPA conc here!
    sfa[2] = 0.5                    # 2 PAI1; calculated from literature
    sfa[3] = visit_df[i,:TAFI]   # 3 TAFI
    sfa[4] = visit_df[i,:AT]     # 4 AT  
    tpa_int = Int64(sfa[1])
    
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = dd.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[i, :II]      # 1 FII
    xₒ[2] = visit_df[i, :Fbgn]    # 2 FI / Fbgn
    xₒ[3] = (1e-14)*SF               # 3 FIIa
    xₒ[6] = visit_df[i, :Plgn]    # 4 Plgn
    dd.initial_condition_array = xₒ

    #update α -
    α = dd.α
    α[1] = 4.0
    α[2] = 0.5
    α[3] = 10.0
    α[4] = 0.03
    α[5] = 0.02

    #update G -
    G = dd.G
    
    FII_idx = findfirst(x->x=="FII",dd.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",dd.total_species_list)
    AT_idx = findfirst(x->x=="AT",dd.total_species_list)
    FI_idx = findfirst(x->x=="FI",dd.total_species_list)
    tPA_idx = findfirst(x->x=="tPA",dd.total_species_list)
    Plgn_idx = findfirst(x->x=="Plgn",dd.total_species_list)
    PAI1_idx = findfirst(x->x=="PAI1",dd.total_species_list)
    Plasmin_idx = findfirst(x->x=="Plasmin",dd.total_species_list)
    TAFI_idx = findfirst(x->x=="TAFI",dd.total_species_list)
    FIa_idx = findfirst(x->x=="FIa",dd.total_species_list)

    # adjusting parameters for r1
    G[FII_idx, 1] = 0.9
    G[FIIa_idx, 1] = 0.9

    # adjusting parameters for r2
    G[FIIa_idx, 2] = 1.1
    G[AT_idx, 2] = 0.05

    # adjusting parameters for r3
    G[FIIa_idx, 3] = 0.5
    G[FI_idx, 3] = 2.0

    # adjusting parameters for r4
    G[tPA_idx,4] = 0.9    
    G[Plgn_idx,4] = 0.8
    G[PAI1_idx,4] = 0.9

    # adjusting parameters for r5
    G[Plasmin_idx,5] = 0.8    
    G[TAFI_idx,5] = 0.1    
    G[FIa_idx,5] = 0.45

    # run the model -
    global (T,U) = evaluate_w_delay(dd,tspan=(0.0,180.0))
    data = [T U]
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,4],sfa[1],U[:,3],xₒ[1])
    
    # dump -
    _PATH_TO_TMP = joinpath(pwd(),"tmp")
    path_to_sim_data = joinpath(_PATH_TO_TMP, "SIM-visit-$(visit)-Fib-$(tpa_int)-nM-tPA-run-$(i).csv")
    CSV.write(path_to_sim_data, Tables.table(hcat(data,CF),header=vcat("Time",dd.list_of_dynamic_species,"CF")))

    # figures -
    _PATH_TO_FIGS = joinpath(pwd(),"figs")
    path_to_CFfigs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_CF_run$(i).png")
    path_to_thrombin_figs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_thrombin_run$(i).png")
    path_to_fibrin_figs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_fibrin_run$(i).png")
    path_to_CF_ensemble_figs = joinpath(_PATH_TO_FIGS, "tPA_$(tpa_int)nM_visit_$(visit)_CF_runs.png")

    Plots.savefig(Plots.plot(T, CF, xticks=0.0:10:180, xlabel="Time (min)", ylabel="CF (mm)", title="Clot firmness vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_CFfigs)
    Plots.savefig(Plots.plot(T, U[:,3], xticks=0.0:10:180,xlabel="Time (min)", ylabel="FIIa (nM)", title="[Thrombin] vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_thrombin_figs)
    #Plots.savefig(Plots.plot(T, U[:,4], xticks=0.0:10:180, xlabel="Time (min)", ylabel="FIa (nM)", title="[Fibrin] vs. time, visit $(visit), [tPA] = $(tpa_int)nM"), path_to_fibrin_figs)
    
end