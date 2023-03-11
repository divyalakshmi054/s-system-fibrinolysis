# load the include -
include("Include.jl")
using Plots

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Fibrinolysis.bst")

# build the default model structure -
model = build(path_to_model_file)

# load the training data -
_PATH_TO_DATA = joinpath(pwd(),"data")
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Synthetic-Composition.csv")
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
    sfa[1] = 8.0                    # 1 tPA
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
    xₒ[9] = 0.125                    # 9 CF
    dd.initial_condition_array = xₒ

    #update α -
    α = dd.α
    α[1] = 0.7
    α[2] = 0.7
    α[3] = 0.06
    α[4] = 0.06
    α[5] = 0.02

    #update G -
    G = dd.G
    
    # what is the index of AT?
    idx = findfirst(x->x=="AT",dd.total_species_list)
    G[idx, 2] = 0.15

    # what is the index of FIIa?
    idx = findfirst(x->x=="FIIa",dd.total_species_list)
    G[idx, 3] = 0.025

    # what is the index of TAFI?
    idx = findfirst(x->x=="TAFI",dd.total_species_list)
    G[idx,5] = -0.1    
    
    #PAI1 -
    idx = findfirst(x->x=="PAI1",dd.total_species_list)
    G[idx,4] = -0.1

    #tPA -
    idx = findfirst(x->x=="tPA",dd.total_species_list)
    G[idx,4] = 0.045

    #Plasmin -
    idx = findfirst(x->x=="Plasmin",dd.total_species_list)
    G[idx,5] = 0.02

    # run the model -
    global (T,U) = evaluate_w_delay(dd,tspan=(0.0,120.0))
    data = [T U]
    
    
    # dump -
    _PATH_TO_TMP = joinpath(pwd(),"tmp")
    path_to_sim_data = joinpath(_PATH_TO_TMP, "SIM-TF-NO-TM-SYN1K-$(i).csv")
    CSV.write(path_to_sim_data, Tables.table(data,header=vcat("Time",dd.list_of_dynamic_species)))
    #_PATH_TO_FIGS = joinpath(pwd(),"figs")
    #path_to_figs = joinpath(_PATH_TO_FIGS, "plot$(i).pdf")
    #Plots.savefig(Plots.plot(T,U[:,9]), path_to_figs)
end