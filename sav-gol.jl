using Random
using Plots

mutable struct SavGolFilter{T}
    buf::Vector{T}
    k::Int
    need_restart::Bool
    coefs::Vector{Int64}
    function SavGolFilter{T}(window::Int) where T
        coefs = [-3, 12, 17, 12, -3] ./ 35
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
        obj.k +=1 
        return 0
    else 
        popfirst!(buf)
        for i in 2:window-1
            buf[i-1]=buf[i]
        buf[window] = x 
        fltrd = [sum(coefs .* buf[1:window])]
        obj.k +=1 
    return fltrd

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

(obj::SavGolFilter)(x) = exe(obj, x)

f_min = 1 
f_max = 100
SNR = 1.0
s_length = 0.1
fs=1000
T = 1 / fs  

flt = SavGolFilter{Float64}(window)

t, sig = signal(s_length,f_min,f_max,SNR)
out = fill(0.0, size(sig))

for i in 1:ceil(s_length/T)
    x = sig[i]
    y = flt(x)
    out[i] = y
plot(t, sig, label="Signal", color=:blue, lw=2)
plot!(t, out, label="Fltrd", color=:red, lw=1)