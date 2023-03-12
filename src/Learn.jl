# it's training time -

function learn_optim(index::Int, model::BSTModel, training_df::DataFrame; 
    pₒ::Union{Nothing,Array{Float64,1}} = nothing)

    # main simulation loop -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = 8.0                        # 1 tPA
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
    xₒ[9] = 0.125                        # 9 CF
    model.initial_condition_array = xₒ
 
    # let's load up our clot parameters sheet -
    _PATH_TO_DATA = joinpath(pwd(),"data")
    parameters_df = CSV.read(joinpath(_PATH_TO_DATA,"Training-Clot-Parameters.csv"),DataFrame)

    # dimensions -
    (R,C) = size(training_df)


    # what is the output array?
    Y = Array{Float64,1}(undef,5)
    Y[1] = parameters_df[index, :CT]
    Y[2] = parameters_df[index, :CFT]
    Y[3] = parameters_df[index, :MCF]
    Y[4] = parameters_df[index, :alpha]
    Y[5] = parameters_df[index, :A30]

# ==== CHANGES ARE ABOVE, STILL WORKING ON THIS FILE ====
# ==== BELOW FROM JV'S COAGULATION CODE, HAVE NOT CHANGED FOR FIBRINOLYTIC MODEL =====

    # setup initial parameter values and bounds array -
    κ = [
            
            # default: hand fit set -
            0.061   0.01 10.0   ; # 1
            1.0     0.01 10.0   ; # 2
            1.0     0.01 10.0   ; # 3
            1.0     0.01 10.0   ; # 4
            1.0     0.01 10.0   ; # 5
            1.0     0.01 10.0   ; # 6
            1.0     0.01 10.0   ; # 7
            1.0     0.01 10.0   ; # 8
            0.70    0.01 10.0   ; # 9
            0.11    0.01 10.0   ; # 10
            0.045   0.01 10.0   ; # 11
            0.065   0.01 10.0   ; # 12
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
    model["α"] = p_best[1:9]
    G = model["G"]
    idx = indexin(model, "FVIIa")   # 10
    G[idx, 4] = p_best[10]
    idx = indexin(model, "AT")      # 11
    G[idx, 9] = p_best[11]
    idx = indexin(model, "TFPI")    # 12
    G[idx, 1] = -1*p_best[12]

    # run the model -
    (T,U) = evaluate(model)
    Xₘ = hcat(U...)
    Yₘ = model_output_vector(T, Xₘ[9,:]) # properties of the Thrombin curve 
    
    return (p_best, T, Xₘ, Yₘ, Y)
end