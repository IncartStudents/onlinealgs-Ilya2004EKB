module Test

include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/sav-gol.jl")
include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/reader.jl")
include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/online.jl")

using .savgol
using .Reader
using .Online
using Random
using SpecialMatrices
using LinearAlgebra
using Pkg
using JSON3
using DataFrames
using Plots
using Statistics

function signal_generator(len::Int, f_min::Int, f_max::Int, SNR::Float64)
    fs = 1000               
    T = 1/fs               
    tmax=len * T - T
    t = collect(0:T:tmax)
    f_num = 10
    freqs = [rand(f_min:f_max) for _ in 1:f_num]  
    amps = [rand(0.1:0.1:1.0) for _ in 1:f_num]  
    signal1 = sum(amps[i] .* sin.(2π * freqs[i] .* t) for i in 1:f_num)
    # signal = [(i-50)^3 + (i-50)^2 + 4*(i-50) + 4 for i in 1:length(t)] # полином второй степени
    signal2 = [(i-50)^5 + (i-50)^4 - (i-50)^3 + (i-50)^2 + 4*(i-50) + 4 for i in 1:length(t)] # полином пятой степени
    noise = (1/SNR) .* randn(length(t)) 
    sig_noise1 = signal1 + noise
    sig_noise2 = signal2 + noise
    return t, sig_noise1, sig_noise2
end

function pseudo_online(window, len, sig) 
    flt = savgol.SavGolFilter{Float64}(window)
    latency = window ÷ 2
    out = zeros(Float64, len)
    for i in 1:latency
        push!(out, 0.0)
        flt(sig[i])
    end
    for i in latency+1:len
        out[i-latency] = flt(sig[i])
    end
    return out
end

function naive(obj::savgol.SavGolFilter{T}, signal, window) where T
    coefs = obj.coefs
    x=signal
    y=zeros(Float64, length(x))
    interv=Int((window-1)/2)
    for k in (interv+1:length(x)-interv)
        y[k] = sum(coefs .* x[k-interv:k+interv])
    end
    return y
end

println("TEST MODE")
function test_stage1(window, len, SNR)
    flt = savgol.SavGolFilter{Float64}(window)
    println("Stage 1. Тестовые сигналы ")
    t, sig1, sig2 = signal_generator(len, 1, 100, SNR)
    println("-----------------")
    println()
    println("сигнал 1. Случайный сигнал")
    error = 0
    out_naive = naive(flt, sig1, window)

    latency = window ÷ 2
    out_online = fill(0.0, len)
    error = 0 
    for i in 1:len
        x = sig1[i]
        y = flt(x)
        if i ≥ window
            out_online[i-latency] = y
            accuracy = (out_online[i-latency]/out_naive[i-latency])
            if abs(accuracy) > 1.05 || abs(accuracy) < 0.95 
                error += 1
            end 
        end
    end
    accuracy_rate1 = 1 - (error/len)
    println("Сигнал 1. Точность:  ", accuracy_rate1)
    display(plot(t, sig1, label="Signal", color=:blue, lw=2, show = true))
    display(plot!(t, out_naive, label="naive", color=:red, lw=2))
    display(plot!(t, out_online, label="online", color=:orange, lw=1))


    println("Сигнал 2. Полином")
    flt = savgol.SavGolFilter{Float64}(window)
    error = 0
    out_naive = naive(flt, sig2, window)
    out_online = fill(0.0, len)
    error = 0 
    for i in 1:len
        x = sig2[i]
        y = flt(x)
        if i ≥ window
            out_online[i-latency] = y
            accuracy = (out_online[i-latency]/out_naive[i-latency])
            if abs(accuracy) > 1.05 || abs(accuracy) < 0.95 
                error += 1
            end 
        end
    end
    accuracy_rate2 = 1 - (error/len)
    println("Сигнал 2. Точность:  ", accuracy_rate2)
    return accuracy_rate1, accuracy_rate2
end

function test_stage2(window, start, len, channel, filepath)
    println("Stage 2. Быстродействие")
    println("-----------------")
    println()
    samples = start:start + len -1 
    name, signal,t =Reader.readbin(filepath, channel, samples) 
    sig = signal .- mean(signal)
    flt = savgol.SavGolFilter{Float64}(window)
    time_naive = @elapsed naive(flt, sig, window)
    # time_online = @elapsed pseudo_online(window, len, sig)
    time_online = @elapsed Online.online(len, window, filepath, channel, "speedtest.bin", 20, start)
    println("Naive : $time_naive seconds")
    println("online : $time_online seconds")
    return time_naive, time_online
end

function test_stage3(window, len, interrupt_point, channel, filepath, filepath_result, batch_len, start_point)
    println("Stage 3. Корректность, длинный сигнал")
    println("-----------------")
    println()
    Online.online(len, window, filepath, channel, filepath_result, batch_len, start_point)

    filtered1 = Reader.openbin_results(1:1199, "data1.bin")
    filtered2 = Reader.openbin_results(1:1199, "data2.bin")
    rawdata = Reader.openbin_results(1:1199, "RAW.bin")

    error = 0

    for k in 1:length(filtered1)-1
        delta = abs(filtered1[k]/filtered2[k])
        if delta != 1 
            error += 1 
        end
    end
    accuracy_rate3 = 1 - (error/length(filtered1))

    println("Корректность, длинный сигнал. Точность: ", accuracy_rate3)
    display(plot(filtered1, label="Fltrd1", color=:blue, lw=2, show = true))
    display(plot!(filtered2, label="Fltrd2", color=:red, lw=2, show = true))
    display(plot!(rawdata, label="Signal", color=:red, lw=1, show = true))
    return nothing
end

end