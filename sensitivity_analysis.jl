# load the include -
include("Include.jl")

# build performance function -
function performance(κ, model::BSTModel, visit_df::DataFrame, i::Int64)

    # main simulation -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = 4.0                    # 1 tPA     SET tPA conc here!
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
    tmp_alpha = κ[1:5]
    g = κ[6:end]

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
    G[TAFI_idx,5] = -1*g[11]  
    G[FIa_idx,5] = g[12]

    # put it back -
    model.G = G;

    # solve -
    (T,U) = evaluate_w_delay(model, tspan = (0.0, 180.0))
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,4],sfa[1],U[:,3],xₒ[1])
    idx = findfirst(x->x==90,T)

    # test -
    return integrate(T[1:idx], CF[1:idx])    # AUC
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

a = [0.5, 0.5, 0.5, 0.02, 0.03]

#update G -
# G = model.G
g = [0.9, 0.9, 1.1, 0.05, 0.5, 2.0, 0.9, 0.8, 0.9, 0.75, 0.75, 0.9] # look at sample_ensemble.jl for specific G values

# fusion -
parameters = vcat(a,g)

np = length(parameters)

L = zeros(np)
U = zeros(np)
for pᵢ ∈ 1:(np)
    L[pᵢ] = 0.1*parameters[pᵢ]
    U[pᵢ] = 1.0*parameters[pᵢ]
end
#L[end] = -3.0;
#U[end] = 0.0;

patient_index = 1;
samples = 10000;

# setup call to Morris method -
F(parameters) =  performance(parameters, model, visit_df, patient_index)
# m = morris(F, L, U, number_of_samples=10000);
m = gsa(F, Morris(num_trajectory=samples), [[L[i],U[i]] for i in 1:np], relative_scale = true);
means = transpose(m.means)
means_star =  transpose(m.means_star)
variances = transpose(m.variances)
results_array = hcat(means, means_star, variances)

# dump sensitivity data to disk -
 CSV.write(joinpath(pwd(),"data","Sensitivity-Morris-test-$(samples).csv"), Tables.table(results_array), header = vcat("mean", "mean_star","variance"))