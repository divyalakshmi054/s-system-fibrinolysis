function amplitude(time::Float64,fibrin::Float64,tPA::Float64,FIIa::Float64)::Float64
    
    a0 = 1.5                    # base amplitude
    kx = 5100 - 615*tPA         # tPA function
    fx = fibrin                 # fibrin species
    td = 10                     # delay parameter

    if(time>td)
        tau = 0.0035*(1-(FIIa/1400))
        S = 21.6+1.08*tPA
        a1 = S*(1 - exp(-tau*(time-td)))
    else
       a1 = 0 
    end
    
    CF = a0 + a1*(fx^2/(fx^2+kx^2))
    return CF
end