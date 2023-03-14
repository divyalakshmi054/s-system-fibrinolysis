function amplitude(time::Array{Float64,1},fibrin::Array{Float64,1},tPA::Float64,thrombin::Array{Float64,1},prothrombin::Float64)

    Rt = length(time)
    a0 = 0.125                              # base amplitude
    kx = 5100 - 615*tPA                     # tPA function
    td = 5                                  # delay parameter

    CF = Array{Float64}(undef,(Rt,))        # clotting firmness
    

    for i ∈ 1:Rt
        fx = fibrin[i]
        FIIa = thrombin[i]

        if(time[i]>td)
            tau = 0.0035*(1-(FIIa/prothrombin))
            S = 76.2+0.4*tPA
            a1 = S*(1 - exp(-tau*(time[i]-td)))
        else
            a1 = 0 
        end
    
        x = a0 + a1*(fx^2/(fx^2+kx^2))
        CF[i] = x

    end

    return CF
end