# load the include -
include("Include.jl")

# build performance function -
function performance(κ, model::BSTModel, visit_df::DataFrame,i::Int64)

    # main simulation -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = 0.0                    # 1 tPA     SET tPA conc here!
    sfa[2] = 0.5                    # 2 PAI1    calculated from literature
    sfa[3] = visit_df[i,:TAFI]      # 3 TAFI
    sfa[4] = visit_df[i,:AT]        # 4 AT  
     
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[i, :II]      # 1 FII
    xₒ[2] = visit_df[i, :Fbgn]    # 2 FI / Fbgn
    xₒ[3] = (1e-14)*SF            # 3 FIIa
    xₒ[6] = visit_df[i, :Plgn]    # 4 Plgn
    model.initial_condition_array = xₒ
    
    #get the parameters -
    tmp_alpha = κ[1:end-12]
    g = κ[end-11:end]

    # set new parameters -
    model.α = tmp_alpha;

    # set G values -
    G = model.G;

    FII_idx = findfirst(x->x=="FII",model.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",model.total_species_list)
    AT_idx = findfirst(x->x=="AT",model.total_species_list)
    FI_idx = findfirst(x->x=="FI",model.total_species_list)
    tPA_idx = findfirst(x->x=="tPA",model.total_species_list)
    Plgn_idx = findfirst(x->x=="Plgn",model.total_species_list)
    PAI1_idx = findfirst(x->x=="PAI1",model.total_species_list)
    Plasmin_idx = findfirst(x->x=="Plasmin",model.total_species_list)
    TAFI_idx = findfirst(x->x=="TAFI",model.total_species_list)
    FIa_idx = findfirst(x->x=="FIa",model.total_species_list)

    # adjusting parameters for r1
    G[FII_idx, 1] = g[1]
    G[FIIa_idx, 1] = g[2]
 
    # adjusting parameters for r2
    G[FIIa_idx, 2] = g[3]
    G[AT_idx, 2] = g[4]
 
    # adjusting parameters for r3
    G[FIIa_idx, 3] = g[5]
    G[FI_idx, 3] = g[6]
 
    # adjusting parameters for r4
    G[tPA_idx,4] = g[7]    
    G[Plgn_idx,4] = g[8]
    G[PAI1_idx,4] = g[9]
 
    # adjusting parameters for r5
    G[Plasmin_idx,5] = g[10]   
    G[TAFI_idx,5] = g[11]  
    G[FIa_idx,5] = g[12]

    model.G = G;

    # solve -
    (T,U) = evaluate_w_delay(model, tspan = (0.0, 180.0))
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,4],sfa[1],U[:,3],xₒ[1])
    clot_firmness = transpose(CF)
    idx = findfirst(x->x==90,T)

    # test -
    return integrate(T[1:idx],clot_firmness[1:idx])    # AUC
end

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

α = model.α
α[1] = 4.0
α[2] = 0.5
α[3] = 10.0
α[4] = 0.03
α[5] = 0.02

#update G -
g = [0.9, 0.9, 1.1, 0.05, 0.5, 2.0, 0.9, 0.8, 0.9, 0.8, 0.1, 0.45] # look at sample_ensemble.jl for specific G values

# fusion -
κ = vcat(α,g)

NP = length(κ)

L = zeros(NP)
U = zeros(NP)
for pᵢ ∈ 1:(NP)
    L[pᵢ] = 0.1*κ[pᵢ]
    U[pᵢ] = 10.0*κ[pᵢ]
end
#L[end] = -3.0;
#U[end] = 0.0;

# setup call to Morris method -
F(κ) =  performance(κ, model, visit_df, 1)
m = morris(F, L, U);

# m = gsa(F, Morris(num_trajectory=10000), [[L[i],U[i]] for i in 1:(NP)]);

# dump sensitivity data to disk -
# jldsave("Sensitivity-Morris-P1-V4.jld2"; mean=m.means, variance=m.variances)