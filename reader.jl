module Reader

include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/sav-gol.jl")

using .savgol
using Random
using SpecialMatrices
using LinearAlgebra
using Pkg
using JSON3
using DataFrames
using Plots
using Statistics

mutable struct Batch{T}
    k::Int              
    marker::Int
    batch_size::Int
    function Batch{T}(batch_size::Int) where T
        new(0, 1, batch_size)
    end
end


function readbin(filepath, ch, range::Union{Nothing, UnitRange{Int}} = nothing)
    Nch, fs, lsbs, names= readhdr(filepath * ".hdr")[1:4]
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
    ln = split(lines[2], delim)
    limit = parse(Int, ln[2])

    lsbs =  parse.(Float64, split(lines[4], delim))
    names = String.(split(lines[3], delim)) 
    return Nch, fs,lsbs, names, limit
end

function openbin_result(filepath_result)
    filtered = []
    open(filepath_result, "r") do io
        read!(io, filtered)
    end
    return filtered
end

function openbin_results(range::Union{Nothing, UnitRange{Int}} = nothing, filename::AbstractString = nothing)
    offset = (range!== nothing) ? range.start - 1 : 0
    elsize = sizeof(Float64) 
    byteoffset = offset * elsize  
    maxlen = (filesize(filename)-byteoffset) ÷ elsize
    len = (range !== nothing) ? min(maxlen, length(range)) : maxlen
    if len <=0 
        data = Vector{Float64}(undef, 0)
    else 
        data = Vector{Float64}(undef, len)
        open(filename, "r") do io 
            seek(io, byteoffset)
            data = reinterpret(Float64, read(io))
        end
    end
    return data
end


function writebin_result(filepath_result, values::Vector{Float64})
    open(filepath_result, "a+") do io
        write(io, values)
    close(io)
    end
end


function savestate(batch::Batch{Float64})
    open("state.json", "w") do file
        JSON3.write(file, batch)
    end
    if (batch.marker / batch.batch_size) % 20 == 0
        println("Batch N: ", (batch.marker ÷ batch.batch_size), " – SUCCESS!")
    end
end

function savestate_filter(savgol)
    open("filter.json", "w") do file
        JSON3.write(file, savgol)
    end
end

function openstate_filter()
    try
        if !isfile("filter.json") 
            create_emptyjson()
            return 0 
        end
        json_string = open("filter.json", "r") do io
            read(io, String)
        end
        filt = JSON3.read(json_string, Dict{String, Any})
        if isempty(filt)  
            return 0  
        else
            coefs = filt["coefs"]
            k = filt["k"]
            buf = filt["buf"]
            savgoley = savgol.SavGolFilter{Float64}(length(buf))
            savgoley.k = k
            savgoley.buf = buf
            savgoley.coefs = coefs
            savgoley.need_restart = false
            return savgoley
        end
    catch e
        println(e)
        return 0  
    end
end

function create_emptyjson()
    filename = "state.json"
    JSON3.open(filename, "w") do file
        JSON3.close(file)
        println("SUCCES!!!!")
    end
    return nothing
    end

function openstate()
    try
        if !isfile("state.json") 
            println("Создание JSON")
            create_emptyjson()
            return 0 
        end
        json_string = open("state.json", "r") do io
            read(io, String)
        end
        postavka = JSON3.read(json_string, Dict{String, Any})
        if isempty(postavka) 
            println("файл пуст")
            return 0  
        else
            batch_size = postavka["batch_size"]
            k = postavka["k"]
            marker = postavka["marker"]
            batch = Batch{Float64}(batch_size)
            batch.k = k
            batch.marker = marker
            println("загрузка состояния - успешно")
            return batch
        end
    catch e
        println(e)
        return 0  
    end
end


end