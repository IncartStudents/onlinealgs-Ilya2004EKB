module savgol

using Random
using SpecialMatrices
using LinearAlgebra
using Pkg
using JSON3
using DataFrames
using Plots
using Statistics

mutable struct SavGolFilter{T}
    buf::Vector{T}
    k::Int
    need_restart::Bool
    coefs::Vector{Float32}
    function SavGolFilter{T}(window::Int) where T
        coefs = vec(calc_coef(window)[1,:])
        new(fill(T(0), window), 1, true, coefs)
    end
end


function calc_coef(order::Int)
    z = Int.(collect((1 - order)/2:1:(order - 1)/2))
    Jp= Vandermonde(z)
    J = Jp[:, 1:order-1]
    C = inv((J'*J))*J'
    return(C)
end
        
            
function exe(obj::SavGolFilter{T}, x::T) where T
    buf, k, coefs= obj.buf, obj.k, obj.coefs
    window = length(buf) 
    if obj.need_restart 
        fill!(buf, x)
        obj.need_restart = false
    end 
    if k < window
        buf[k] = x
    else 
        buf[window] = x 
    end
    fltrd = sum(coefs .* buf[1:window])
    for i in 2:(window)
        buf[i-1]=buf[i]
    end
    obj.k +=1 
    return fltrd
end

(obj::SavGolFilter)(x) = exe(obj, x)

end