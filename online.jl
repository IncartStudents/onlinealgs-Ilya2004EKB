module Online
include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/sav-gol.jl")
include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/reader.jl")

using .savgol
using .Reader
using Plots
using Statistics
using JSON3
function online(batch_len, window, filepath, channel, filepath_result, batches, start_point)
    batch = Reader.openstate()
    if batch == 0
        batch = Reader.Batch{Float64}(batch_len)
        batch.marker = start_point
    end
    savgoley = Reader.openstate_filter()
    if savgoley == 0 
        flt = savgol.SavGolFilter{Float64}(window)
    else
        flt = savgoley
    end
    i = 1
    limit= Reader.readhdr(filepath * ".hdr")[5]
    # while i <= limit ÷ batch_len
    while i <= batches
        count = 1
        start = batch.marker
        batch_len = batch.batch_size
        samples = start:start + batch_len - 1 
        name, signal,t= Reader.readbin(filepath, channel, samples) 
        Reader.writebin_result("RAW1.bin", signal)
        # signal = signal .- mean(signal)
        filtered = Vector{Float64}()
        i+=1
        try
            while count <= batch_len
                push!(filtered, flt(signal[count]))
                batch.marker += 1 
                count += 1 
            end
            batch.k +=1
            Reader.writebin_result(filepath_result, filtered)
            Reader.savestate(batch)
            Reader.savestate_filter(flt)
        catch e
            if isa(e, InterruptException)
                Reader.savestate(batch)
                Reader.savestate_filter(flt)
                Reader.writebin_result(filepath_result, filtered)
                println("Прерывание")
            else
                rethrow(e)
            end
            return nothing
        end
    end
end

end