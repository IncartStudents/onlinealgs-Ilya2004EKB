using Random
using Plots

mutable struct SavGolFilter{T}
    buf::Vector{T}
    k::Int
    need_restart::Bool
    coefs::Vector{Float64}
    function SavGolFilter{T}(window::Int) where T
        coefs = [-2, 3, 6, 7, 6, 3, -2] ./ 21
        new(fill(T(0), window), 1, true, coefs)
    end
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
    # signal = sum(amps[i] .* sin.(2π * freqs[i] .* t) for i in 1:f_num)
    signal = [(i-50)^3 + (i-50)^2 + 4*(i-50) + 4 for i in 1:length(t)] # полином второй степени
    signal = [(i-50)^2 + 4*(i-50) + 4 for i in 1:length(t)] # полином третей степени
    noise = (1/SNR) .* randn(length(t)) 
    sig_noise = signal + noise
    return t, sig_noise
end

(obj::SavGolFilter)(x) = exe(obj, x)

f_min = 1 
f_max = 100
SNR = 1.0
s_length = 0.1
fs=1000
T = 1 / fs  
window = 7
latency = window ÷ 2

flt = SavGolFilter{Float64}(window)

t, sig = signal(s_length,f_min,f_max,SNR)

# sig = [64.0, 49.0, 36.0, 25.0, 16.0, 9.0, 4.0, 1.0, 0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0] # y = xˆ2 + 4x + 4
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
plot!(out, label="Fltrd", color=:red, lw=1)