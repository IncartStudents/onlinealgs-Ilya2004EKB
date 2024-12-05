include("/Users/mac/Documents/GitHub/onlinealgs-Ilya2004EKB/sav-gol.jl")
using .savgol
using Plots
using Statistics

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
file_name = "bidmc_01_Signals.csv"
channel = 5

f_min = 1 
f_max = 100
SNR = 30
# s_length = 200
# signal,t = readfile(file_path, channel)
# noise = (1/SNR) .* randn(length(t))
# sig_noise = sig_eeg + noise

filepath = "/Users/mac/Downloads/data/All/all_MX120161018125923"
start = 100
len = 200
samples = start:start + len -1 

name, signal,t = savgol.readbin(filepath, 20, samples) 
signal = signal .- mean(signal)
println(t)
println(length(t))
window = 9
latency = window ÷ 2

(obj::savgol.SavGolFilter)(x) = savgol.exe(obj, x)

flt = savgol.SavGolFilter{Float64}(window)

# t, sig = signal(s_length,f_min,f_max,SNR)
out2=savgol.naive(flt, signal, window)
out = fill(0.0, size(signal))
for i in 1:len
    # String = Vector{Float64}()
    # read!(file_name, String) 
    # println(String)
    x = signal[i]
    y = flt(x)
    if i ≥ window
        out[i-latency] = y
    end
end
# println(t)
# println(signal)
plot(t, signal, label="Signal", color=:blue, lw=2)
plot!(t, out, label="Fltrd", color=:red, lw=2)
plot!(t, out2, label="Fltrd_naive", color=:orange, lw=1)




# time_naive = @elapsed naive(flt, sig, window)
# time_online = @elapsed pseudo_online(sig_noise, window)
# println("Naive : $time_naive seconds")
# println("online : $time_online seconds")
