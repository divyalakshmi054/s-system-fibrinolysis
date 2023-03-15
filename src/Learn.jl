# it's training time -

function learn_optim(index::Int, model::BSTModel, training_df::DataFrame; 
    pₒ::Union{Nothing,Array{Float64,1}} = nothing)

    # main simulation loop -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = 0.0                        # 1 tPA
    sfa[2] = 0.5                        # 2 PAI1; calculated from literature
    sfa[3] = training_df[index,:TAFI]   # 3 TAFI
    sfa[4] = training_df[index,:AT]     # 4 AT   
     
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = training_df[index, :II]      # 1 FII
    xₒ[2] = training_df[index, :Fbgn]    # 2 FI / Fbgn
    xₒ[3] = (1e-14)*SF                   # 3 FIIa
    xₒ[6] = training_df[index, :Plgn]    # 6 Plgn
    model.initial_condition_array = xₒ
 
   # load the training data -
    path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Clot-Parameters.csv")
    full_df = CSV.read(path_to_training_data, DataFrame)

    # let's filter visit 4 & no TPA -
    visit_df = filter(:Visit => x->(x=="V4"), full_df)
    tPA_df = filter(:TPA=>x->(x=="N"), visit_df)
    learn_df = filter(:CT => x->(x!=1000), tPA_df)


    # dimensions -
    (R,C) = size(training_df)


    # what is the output array?
    Y = Array{Float64,1}(undef,6)
    Y[1] = learn_df[index, :CT]
    Y[2] = learn_df[index, :CFT]
    Y[3] = learn_df[index, :MCF]
    Y[4] = learn_df[index, :MCFt]
    Y[5] = learn_df[index,:alpha]
    Y[6] = learn_df[index, :A30]

    # setup initial parameter values and bounds array -
    κ = [
            
            # default: hand fit set -
            1.0         0.01 10.0       ; # 1
            0.2         0.01 10.0       ; # 2
            1.0         0.01 10.0       ; # 3
            0.025       0.01 10.0       ; # 4
            0.01        0.01 10.0       ; # 5
            0.9         0.01 10.0       ; # 6
            0.7         0.01 10.0       ; # 7
            0.15        0.01 10.0       ; # 8
            0.75        0.01 10.0       ; # 9
            0.75        0.01 10.0       ; # 10
            0.1         0.01 10.0       ; # 11
            0.2         0.01 10.0       ; # 12
        ];

    # set default set as the start -
    if (isnothing(pₒ) == true)
        P = length(κ[:,1])
        σ = 0.1 # move up to 10%
        pₒ = κ[:,1].*(1 .+ σ*rand(-1:1,P))
    end

    # setup the objective function -
    inner_optimizer = NelderMead()
    obj_function(p) =  loss_scalar(p, Y, model)
    results = optimize(obj_function, κ[:,2], κ[:,3], pₒ, Fminbox(inner_optimizer),
        Optim.Options(time_limit = 600, show_trace = true, show_every = 10, iterations=100))
    
    # grab the best parameters -
    p_best = Optim.minimizer(results)
    
    # run the sim w/the best parameters -
    # 1 - 9 : α vector
    model.α = p_best[1:5]
    G = model.G
    idx = findfirst(x->x=="FII",dd.total_species_list)
    G[idx, 1] = p_best[6]
    idx = findfirst(x->x=="FIIa",dd.total_species_list)
    G[idx, 1] = p_best[7]
    idx = findfirst(x->x=="AT", model.total_species_list)   # 6
    G[idx, 2] = p_best[8]
    idx = findfirst(x->x=="FIIa", model.total_species_list) # 7
    G[idx, 3] = p_best[9]
    idx = findfirst(x->x=="FI", model.total_species_list)   # 8
    G[idx, 3] = p_best[10]
    idx = findfirst(x->x=="TAFI",model.total_species_list)  # 9
    G[idx,5] = p_best[11]
    idx = findfirst(x->x=="FIa",model.total_species_list)   # 10
    G[idx,5] = p_best[12]

    # run the model -
    (T,U) = evaluate(model)
    CF = Array{Float64,1}
    CF = amplitude(T,U[:,4],sfa[1],U[:,3],xₒ[1])
    Xₘ = hcat(CF)
    Yₘ = model_output_vector(T, CF) # properties of the CF curve 
    
    return (p_best, T, Xₘ, Yₘ, Y)
end