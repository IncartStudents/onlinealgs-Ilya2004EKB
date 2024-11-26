using Random
using Plots
using SpecialMatrices
using LinearAlgebra

mutable struct SavGolFilter{T}
    buf::Vector{T}
    k::Int
    need_restart::Bool
    coefs::Vector{Float64}
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
        buf[1] = x 
        obj.need_restart = false
    end 
    if k < window
        buf[obj.k] = x
        fltrd = 0
    else 
        buf[window] = x 
        fltrd = sum(coefs .* buf[1:window])
        for i in 2:(window)
            buf[i-1]=buf[i]
        end
    end
    obj.k +=1 
    return fltrd

end

function signal(len::Float64, f_min::Int, f_max::Int, SNR::Float64)
    fs = 1000               
    T = 1/fs               
    t = collect(0:T:len)
    f_num = 10
    freqs = [rand(f_min:f_max) for _ in 1:f_num]  
    amps = [rand(0.1:0.1:1.0) for _ in 1:f_num]  
    signal = sum(amps[i] .* sin.(2π * freqs[i] .* t) for i in 1:f_num)
    # signal = [(i-50)^3 + (i-50)^2 + 4*(i-50) + 4 for i in 1:length(t)] # полином второй степени
    signal = [(i-50)^2 + 4*(i-50) + 4 for i in 1:length(t)] # полином третей степени
    noise = (1/SNR) .* randn(length(t)) 
    sig_noise = signal + noise
    return t, sig_noise
end

(obj::SavGolFilter)(x) = exe(obj, x)

f_min = 1 
f_max = 100
SNR = 0.05
s_length = 0.1
fs=1000
T = 1 / fs  
window = 11
latency = window ÷ 2

flt = SavGolFilter{Float64}(window)

t, sig = signal(s_length,f_min,f_max,SNR)

out = fill(0.0, size(sig))
for i in 1:length(sig)
    i=Int64(i)
    x = sig[i]
    y = flt(x)
    if i ≥ window
        out[i-latency] = y
    end
end
plot(sig, label="Signal", color=:blue, lw=2)
print(calc_coef(5)[1,:])
plot!(out, label="Fltrd", color=:red, lw=2)