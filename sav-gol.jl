module savgol

using Random
using SpecialMatrices
using LinearAlgebra
using Pkg
using JSON3
using CSV
using DataFrames
using Plots
using Statistics

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

mutable struct batch{T}
    data::Vector{T}   
    k::Int              
    marker::Int
    batch_size::Int     
    function batch{T}(batch_size::Int) where T
        new(Vector{T}(), 0, 1, batch_size)
    end
end

function readhdr(filepath::AbstractString)
    str = open(filepath, "r") do io 
        len = stat(filepath).size
        bytes = Vector{UInt8}(undef, len)
        readbytes!(io, bytes, len)
        if bytes[1:3] == [0xEF, 0xBB, 0xBF]
            bytes = bytes[4:end]
        end
        str = Array{Char}(bytes) |> x-> String(x)
    end
    io = IOBuffer(str)
    lines = readlines(io)
    lines = rstrip.(lines)
    delim = (' ', '\t')
    ln = split(lines[1], delim)
    Nch = parse(Int, ln[1])
    fs = parse(Float64, ln[2])
    lsbs =  parse.(Float64, split(lines[4], delim))
    names = String.(split(lines[3], delim)) 
    return Nch, fs,lsbs, names
end

function savestate(batch)
    open("state.json", "w") do file
        JSON3.write(file, batch)
    end
    if (batch.marker / batch.batch_size) % 20 == 0
        println("Batch N: ", (batch.marker ÷ batch.batch_size), " – SUCCESS!")
    end
end

function openstate(len, reset_batch_flag = 0)
    if reset_batch_flag != 0
        batch = JSON3.read("state.json", Batch{Float64})
    else
        batch = batch{Float64}(len)
    end
    return batch
end

function batch_append(batch::batch{T}, value) where T
        push!(batch.data, value)
        batch.k += 1
        if batch.k == batch.batch_size
            batch.marker = batch.marker + batch.batch_size
            savestate(batch)
            batch.data = Vector{T}()
            batch.k = 0
        end
end


function readbin(filepath, ch, range::Union{Nothing, UnitRange{Int}} = nothing)
    Nch, fs, lsbs, names = readhdr(filepath * ".hdr")
    offset = (range!== nothing) ? range.start - 1 : 0
    elsize = sizeof(Int32) * Nch
    byteoffset = offset * elsize  
    maxlen = (filesize(filepath .* ".bin")-byteoffset) ÷ elsize
    len = (range !== nothing) ? min(maxlen, length(range)) : maxlen
    if len <=0 
        data = Matrix{Int32}(undef, Nch, 0)
    else 
        data = Matrix{Int32}(undef, Nch, len)
        open(filepath .* ".bin", "r") do io 
            seek(io, byteoffset)
            read!(io, data)
        end
    end
    channels = [(data[i, :] .* lsbs[i]) for i in 1:Nch] |> Tuple 
    sig = channels[ch]
    T =1/fs
    t = range .* T 
    return names[ch], sig, t
end


# function readfile(file_path::AbstractString, ch::Int)
#     data = CSV.read(file_path, DataFrame; header = 1)
#     signal = data[!, ch+1]
#     return signal
# end


function calc_coef(order::Int)
    z = Int.(collect((1 - order)/2:1:(order - 1)/2))
    Jp= Vandermonde(z)
    J = Jp[:, 1:order-1]
    C = inv((J'*J))*J'
    return(C)
end
        
function naive(obj::SavGolFilter{T}, signal, window) where T
    coefs = obj.coefs
    x=signal
    y=zeros(Float64, length(x))
    interv=Int((window-1)/2)
    for k in (interv+1:length(x)-interv)
        y[k] = sum(coefs .* x[k-interv:k+interv])
    end
    return y
end
            
function exe(obj::SavGolFilter{T}, x::T, batch) where T
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
    batch_append(batch, fltrd)
    return fltrd
end

(obj::SavGolFilter)(x) = exe(obj, x, batch)

function online(window, batch_len)
    openstate(batch_len, 0)
    start = batch.marker
    batch_len = batch.batch_size
    samples = start:start + batch_len -1 
    batch = batch{Float64}(batch_len)
    name, signal,t = readbin(filepath, channel, samples) 
    signal = signal .- mean(signal)
    out = zeros(Float64, batch_len)

    flt = SavGolFilter{Float64}(window)
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

function pseudo_online(window, len, sig) 
    flt = SavGolFilter{Float64}(window)
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


function signal_generator(len::Int, f_min::Int, f_max::Int, SNR::Float64)
    fs = 1000               
    T = 1/fs               
    t = collect(0:T:len)
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

function main(test = 0, batch_len = 200, naive_incl = 0, window = 9)
    flt = SavGolFilter{Float64}(window)

    if test != 0
        println("TEST MODE")
        println("Stage 1. Тестовые сигналы ")
        len = 300
        t, sig1, sig2 = signal_generator(len,1,100,10.0)
        println("-----------------")
        println()
        println("сигнал 1. Полином")
        error = 0
        out_naive = naive(flt, sig2, window)

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
        accuracy_rate = 1 - (error/len)
        println("Сигнал 1. Точность:  ", accuracy_rate)
        println()
        
        plot(t, sig1, label="Signal", color=:blue, lw=2)
        plot(t, out_naive, label="naive", color=:red, lw=2)
        plot!(t, out_online, label="online", color=:orange, lw=1)

        println()
        println("Сигнал 2. Случайный сигнал")

        flt = SavGolFilter{Float64}(window)
        error = 0
        out_naive = naive(flt, sig1, window)
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
        accuracy_rate = 1 - (error/len)
        println("Сигнал 2. Точность:  ", accuracy_rate)

        println("Stage 2. Быстродействие")
        println("-----------------")
        println()
        filepath = "/Users/mac/Downloads/data/All/all_MX120161018125923"
        start = 1
        len = 500000
        channel = 20
        samples = start:start + len -1 
        name, signal,t = readbin(filepath, channel, samples) 
        sig = signal .- mean(signal)

        time_naive = @elapsed naive(flt, sig, window)
        time_online = @elapsed pseudo_online(window, len, sig)
        println("Naive : $time_naive seconds")
        println("online : $time_online seconds")
    end
end

end