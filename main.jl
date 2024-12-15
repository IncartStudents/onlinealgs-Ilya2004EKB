

include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/reader.jl")
include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/test.jl")
include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/sav-gol.jl")
include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/online.jl")

using .Reader, .savgol, .Test, .Online, Plots


function main(test, batches, naive_incl, window, start, channel, len_stage1, len_stage2, len_stage3, SNR, filepath, filepath_result, batch_len)
    if test != 0
        Test.test_stage1(window, len_stage1, SNR)
        Test.test_stage2(window, start, len_stage2, channel, filepath)
        Test.test_stage3(window, len_stage3, interrupt_point, channel, filepath, filepath_result, batches, start)
    else 
        Online.online(batch_len, window, filepath, channel, filepath_result, batches, start)
        sig_len = batch_len * batches
        latency = window ÷ 2 
        filtered = Reader.openbin_results(1:399, "data.bin")
        rawdata = Reader.openbin_results(1:399, "RAW.bin")
        display(plot(rawdata[1:sig_len], label="Signal", color=:blue, lw=2, show = true))
        display(plot!(filtered[latency:sig_len], label="Fltrd", color=:red, lw=1, show = true))
    end
    return nothing
end

test = 1
batches = 20
batch_len
include_naive = 1 
window = 9
start_point = 10000
channel = 20
test1_len = 300
test2_len = 500000
test3_len = 30
interrupt_point = batch_len * 1000 + 50 
SNR = 10.0
filepath = ("/Users/mac/Downloads/data/All/all_MX120161018125923")
filepath_result = "data.bin" 

main(test, batches, include_naive, window, start_point, channel, test1_len, test2_len, test3_len, SNR, filepath, filepath_result, batch_len)