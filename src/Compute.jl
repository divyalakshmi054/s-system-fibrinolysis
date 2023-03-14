# let's get some properties of the CF v t curve to help with our training -

# load the include -
include("Include.jl")

function clot_properties()
    # load the training data -
    _PATH_TO_DATA = joinpath(pwd(),"data")
    path_to_training_data = joinpath(_PATH_TO_DATA, "Training-TEG.csv")
    unclean_training_df = CSV.read(path_to_training_data, DataFrame, copycols = true)

    # dimensions -
    (R,C) = size(unclean_training_df)
    training_df = Array{Float64}(undef,(R,C))

    # some data is 'missing' - let's replace that with 0 for now
    for i ∈ 1:R
        for j ∈ 1:C
            if(ismissing(unclean_training_df[i,j]))
                training_df[i,j] = 0
            else
                training_df[i,j] = unclean_training_df[i,j]
            end
        end
    end

    T = training_df[:,1]                            # define time
    data_vector = Array{Float64}(undef,(6,(C-1)))   # 5 properties, C-1 patients

    for i ∈ 1:(C-1)                                 #omitting column 1: time
        X = training_df[1:end,i+1]
        data_vector[:,i] = model_output_vector(T,X)
        print(data_vector[:,i],"\n")
    end
    print(size(data_vector))
    #time to dump, finally! -
    data_output_header = ["CT", "CFT", "MCF","MCFt", "alpha","A30"]
    CSV.write(joinpath(_PATH_TO_DATA,"Training-Clot-Parameters.csv"),Tables.table(transpose(data_vector));header = data_output_header)
end

function model_output_vector(T::Array{Float64,1},X::Array{Float64,1})::Array{Float64,1}
    output_vector = Array{Float64,1}()

    #properties -
    CT = clot_time(x->x>=2.0, T, X)
    CFT = clot_formation_time(x->x>=20.0, T, X) - CT
    MCF = maximum(X)
    MCF_t = MCF_time(MCF,T,X)
    alpha_slope = compute_alpha_slope(CFT, CT, T, X)*180/pi
    a_30 = compute_amp_30(CT,T,X)

    #it's go time -
    push!(output_vector,CT)
    push!(output_vector,CFT)
    push!(output_vector,MCF)
    push!(output_vector,MCF_t)
    push!(output_vector,alpha_slope)
    push!(output_vector,a_30)
    return output_vector
end

function clot_time(rule::Function, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(rule, X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end

function clot_formation_time(rule::Function, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(rule, X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end

function compute_alpha_slope(CFT::Float64, CT::Float64, T::Array{Float64,1}, X::Array{Float64,1})::Float64
    if(CT==1000.0)
        return 1000.0
    end
    idx_CFT = findfirst(x->x==CFT+CT, T)
    idx_CT = findfirst(x->x==CT, T)
    slope = (X[idx_CFT] - X[idx_CT])/(CFT)
    alpha_angle = atan(slope)
    return alpha_angle
end

function compute_amp_30(CT::Float64, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(x->x>=(CT+30),T)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return X[idx]
end

function MCF_time(MCF::Float64, T::Array{Float64,1}, X::Array{Float64,1})::Float64

    # filter -
    idx = findfirst(x->x==MCF,X)

    if (isnothing(idx) == true)
        return 1000.0
    end

    # return -
    return T[idx]
end