# let's get some properties of the CF v t curve -

# load the include -
include("Include.jl")

function clot_properties()::Array{Float64,1}
    # load the training data -
    _PATH_TO_DATA = joinpath(pwd(),"data")
    path_to_training_data = joinpath(_PATH_TO_DATA, "Training-TEG.csv")
    training_df = CSV.read(path_to_training_data, DataFrame)

    # dimensions -
    (R,C) = size(training_df)

    T = training_df[:,1]                            # define time
    data_vector = Array{Float64}(undef,(5,C-1))     # output

    for i ∈ 2:C     #omitting column 1: time
        
        X = training_df[1:end,i]
        output_vector = Array{Float64}(undef,(5,1))
        
        #properties -
        CT = clot_time(x->x>=2.0, T, X)
        CFT = clot_formation_time(x->x>=20.0, T, X) - CT
        MCF = maximum(X)
        alpha = alpha_angle(CFT,CT,T,X)
        area_under_curve = auc(T, X)

        push!(output_vector,CFT)
        push!(output_vector,CT)
        push!(output_vector,MCF)
        push!(output_vector,alpha)
        push!(output_vector,area_under_curve)

        # storage -
        data_vector[:,i] = output_vector
    end

    # let's try & dump -
    #data_output_header = ["CFT","CT","MCF","alpha","AUC"]
   #CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "Clot-Properties.csv"), 
    #Tables.table(data_output); header = data_output_header)
end

function clot_time(rule::Function, T::Array{Int64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(rule, X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end

function clot_formation_time(rule::Function, T::Array{Int64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(rule, X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end

function alpha_angle(CFT::Float64,CT::Float64, T::Array{Int64,1}, X::Array{Float64,1})::Float64
    idx_CFT = findfirst(x->x==CFT)
    idx_CT = findfirst(x->x==CT)
    slope = (X[idx_CFT] - X[idx_CT])/(CFT - CT)
    alpha = atan(slope)
end

function auc(T::Array{Int64,1}, X::Array{Int64,1})::Float64    
    return integrate(T, X)
end