using Random
using Plots

mutable struct SavGolFilter{T}
    buf::Vector{T}
    k::Int
    need_restart::Bool
    function SavGolFilter{T}(window::Int) where T
        new(fill(T(0), window-1), 1, true)
    end
end



function signal(len::Float64, f_min::Int, f_max::Int, SNR::Float64)
    fs = 1000               
    T = 1 / fs               
    t = collect(0:T:len)
    f_num = 10
    freqs = [rand(f_min:f_max) for _ in 1:f_num]  
    amps = [rand(0.1:0.1:1.0) for _ in 1:f_num]  
    signal = sum(amps[i] .* sin.(2π * freqs[i] .* t) for i in 1:f_num)
    noise = (1/SNR) .* randn(length(t))
    sig_noise = signal + noise
    return t, sig_noise
end


f_min = 1 
f_max = 100
SNR = 1.0
s_length = 0.1
t, sig = signal(s_length,f_min,f_max,SNR)
plot(t, sig, label="Signal", color=:blue, lw=2)