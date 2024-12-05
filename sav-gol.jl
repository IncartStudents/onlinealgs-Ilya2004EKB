module savgol

using Random
using SpecialMatrices
using LinearAlgebra
using Pkg
using JSON3
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
    T=1/fs
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
        fill!(buf, x)
        obj.need_restart = false
    end 
    if k < window
        buf[obj.k] = x
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


end