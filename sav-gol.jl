using Random
using Plots
using SpecialMatrices
using LinearAlgebra
using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
# include("/Users/mac/Downloads/Chrono.jl-master/src/Chrono.jl")
# using Duration
using CSV
using DataFrames

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

function readfile(file_path::AbstractString, ch::Int)
    data = CSV.read(file_path, DataFrame; header=1)
    signal = data[!, ch+1]
    t = data[!, 1]
    return signal, t
end


function calc_coef(order::Int)
    z = Int.(collect((1 - order)/2:1:(order - 1)/2))
    Jp= Vandermonde(z)
    J = Jp[:, 1:order-1]
    C = inv((J'*J))*J'
    return(C)
end
        
function naive(obj::SavGolFilter{T}, signal, window) where T
    coefs = obj.coefs
    # coefs = [-3, 12, 17, 12, -3] ./ 35
    x=signal
    y=zeros(Float64, length(x))
    interv=Int((window-1)/2)
    for k in (interv+1:length(x)-interv)
        y[k] = sum(coefs .* x[k-interv:k+interv])
    end
    return y
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
    end
    fltrd = sum(coefs .* buf[1:window])
    for i in 2:(window)
        buf[i-1]=buf[i]
    end
    obj.k +=1 
    return fltrd
end

function pseudo_online(sig,window)
    len = length(sig)
    out = zeros(Float64, len)
    latency = window ÷ 2
    for i in 1:len
        x = sig[i]
        y = flt(x)
        if i ≥ window
            out[i-latency] = y
        end
    end
    return out
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
    signal = [(i-50)^5 + (i-50)^4 - (i-50)^3 + (i-50)^2 + 4*(i-50) + 4 for i in 1:length(t)] # полином пятой степени
    noise = (1/SNR) .* randn(length(t)) 
    sig_noise = signal + noise
    return t, sig_noise
end

(obj::SavGolFilter)(x) = exe(obj, x)


# https://physionet.org/content/auditory-eeg/1.0.0/Raw_Data/#files-panel ссылка на датасет
# Recording Tools
# OpenBCI Ganglion Board, 200 Hz sampling rate, four channels: T7, F8, Cz, and P4.
file_path = "/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/s01_ex01_s01.txt"


# https://physionet.org/content/bidmc/1.0.0/bidmc_csv/#files-panel
# Physiological signals, such as the PPG, impedance respiratory signal, and electrocardiogram (ECG). These are sampled at 125 Hz.
# Physiological parameters, such as the heart rate (HR), respiratory rate (RR), and blood oxygen saturation level (SpO2). These are sampled at 1 Hz.
# Fixed parameters, such as age and gender
# Manual annotations of breaths.

file_path = "/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/bidmc_01_Signals.csv"
channel = 5

f_min = 1 
f_max = 100
SNR = 30
s_length = 10000
sig_eeg,t = readfile(file_path, channel)

noise = (1/SNR) .* randn(length(t)) 
sig_noise = sig_eeg + noise

fs=125
T = 1 / fs  
# t = collect(0:T:length(sig_eeg))

window = 17
latency = window ÷ 2

flt = SavGolFilter{Float64}(window)

# t, sig = signal(s_length,f_min,f_max,SNR)
sig = sig_noise[1:s_length]
t = t[1:s_length]

out = fill(0.0, size(sig))
for i in 1:length(s_length)
    x = sig[i]
    y = flt(x)
    if i ≥ window
        out[i-latency] = y
    end
end
out2=naive(flt, sig, window)
plot(t,sig, label="Signal", color=:blue, lw=2)
plot!(t,out, label="Fltrd", color=:red, lw=2)
plot!(t, out2, label="Fltrd_naive", color=:orange, lw=1)


time_naive = @elapsed naive(flt, sig, window)
time_online = @elapsed pseudo_online(sig_noise, window)
println("Naive : $time_naive seconds")
println("online : $time_online seconds")
