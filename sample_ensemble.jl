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

# size of training set -
(R,C) = size(training_df)

# main simulation -
SF = 1e9
for i ∈ 1:10

    # build new model -
    dd = deepcopy(model)

    # setup static -
    sfa = dd.static_factors_array
    sfa[1] = 0.0                    # 1 tPA
    sfa[2] = 0.5                    # 2 PAI1; calculated from literature
    sfa[3] = training_df[i,:TAFI]   # 3 TAFI
    sfa[4] = training_df[i,:AT]     # 4 AT   
    
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = dd.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = training_df[i, :II]      # 1 FII
    xₒ[2] = training_df[i, :Fbgn]    # 2 FI / Fbgn
    xₒ[3] = (1e-14)*SF               # 3 FIIa
    xₒ[6] = training_df[i, :Plgn]    # 4 Plgn
    dd.initial_condition_array = xₒ

    #update α -
    α = dd.α
    α[2] = 0.7
    #α[4] = 0.05
    #α[5] = 0.025

    #update G -
    G = dd.G
    
    # what is the index of AT?
    idx = findfirst(x->x=="AT",dd.total_species_list)
    G[idx, 2] = 0.15

    # what is the index of FIIa?
    #idx = findfirst(x->x=="FIIa",dd.total_species_list)
    #G[idx, 3] = 0.025

    # what is the index of FI?
    idx = findfirst(x->x=="FI",dd.total_species_list)
    G[idx, 3] = 0.5

    # what is the index of TAFI?
   # idx = findfirst(x->x=="TAFI",dd.total_species_list)
   # G[idx,5] = 0.01

    # what is the index of fibrin?
   # idx = findfirst(x->x=="FIa",dd.total_species_list)
   # G[idx,5] = 0.05

    # run the model -
    global (T,U) = evaluate_w_delay(dd,tspan=(0.0,180.0))
    data = [T U]
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,4],sfa[1],U[:,3],xₒ[1])
    
    # dump -
    _PATH_TO_TMP = joinpath(pwd(),"tmp")
    path_to_sim_data = joinpath(_PATH_TO_TMP, "SIM-TF-NO-TM-SYN1K-$(i).csv")
    CSV.write(path_to_sim_data, Tables.table(hcat(data,CF),header=vcat("Time",dd.list_of_dynamic_species,"CF")))

    # figures -
    _PATH_TO_FIGS = joinpath(pwd(),"figs")
    path_to_CFfigs = joinpath(_PATH_TO_FIGS, "CFplot$(i).png")
    Plots.savefig(Plots.plot(T,CF), path_to_CFfigs)
    path_to_thrombin_figs = joinpath(_PATH_TO_FIGS, "thrombinplot$(i).png")
    Plots.savefig(Plots.plot(T,U[:,3]), path_to_thrombin_figs)
    path_to_fibrin_figs = joinpath(_PATH_TO_FIGS, "fibrinplot$(i).png")
    Plots.savefig(Plots.plot(T,U[:,4]), path_to_fibrin_figs)
    
end