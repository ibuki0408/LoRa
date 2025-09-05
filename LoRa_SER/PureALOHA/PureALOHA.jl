using Random
using FFTW
using Statistics
using Printf
using DelimitedFiles

# ===== LoRa基本関数 =====
function lora_symbol(sf::Int, bw::Float64, m::Int)
    M  = 2^sf
    tc = 1.0 / bw
    ts = tc * M
    x = Vector{ComplexF64}(undef, M)
    @inbounds for k in 1:M
        x[k] = exp(1im * 2π * (bw/(2*ts)) * ((k-1) + m)^2 * tc^2)
    end
    return x
end

function demod_lora(sf::Int, bw::Float64, x::Vector{ComplexF64})
    M  = 2^sf
    tc = 1.0 / bw
    ts = tc * M
    d = Vector{ComplexF64}(undef, M)
    @inbounds for k in 1:M
        d[k] = x[k] * exp(-1im * 2π * (bw/(2*ts)) * (k-1)^2 * tc^2)
    end
    Y = fft(d)
    return argmax(abs.(Y)) - 1
end

@inline function add_awgn!(y::Vector{ComplexF64}, snr_dB::Float64)
    var = 10^(-snr_dB/10)
    noise = sqrt(var/2)
    @inbounds @simd for i in eachindex(y)
        y[i] += ComplexF64(randn()*noise, randn()*noise)
    end
    return y
end

# ===== Pure ALOHAシミュレーション関数 =====
function purealoha_ser(sf::Int, bw::Float64, num_devices::Int;
                       payload_range::Tuple{Int,Int}=(4,32),
                       sim_time::Float64=1.0,
                       snr_dB::Float64=0.0,
                       rng=Random.default_rng())
    M = 2^sf
    Tsym = M / bw
    fs = bw  # サンプル周波数（簡易）

    # 各端末の送信開始時間とペイロード長
    tx_starts = rand(rng, num_devices) .* sim_time
    payload_lens = rand(rng, payload_range[1]:payload_range[2], num_devices)
    payloads = [rand(rng, 0:M-1, payload_lens[d]) for d in 1:num_devices]
    durations = payload_lens .* Tsym

    # 受信信号バッファ
    total_samples = Int(ceil(sim_time*fs)) + maximum(payload_lens)*M
    rx_signal = zeros(ComplexF64, total_samples)

    # 波形重畳
    for d in 1:num_devices
        for s in 1:payload_lens[d]
            sym_wave = lora_symbol(sf, bw, payloads[d][s])
            start_idx = Int(floor((tx_starts[d] + (s-1)*Tsym)*fs)) + 1
            end_idx = min(start_idx + length(sym_wave) - 1, length(rx_signal))
            if start_idx <= length(rx_signal) && end_idx >= start_idx
                rx_signal[start_idx:end_idx] .+= sym_wave[1:end_idx-start_idx+1]
            end
        end
    end

    # AWGN付加
    add_awgn!(rx_signal, snr_dB)

    # 復調と誤り判定
    errors = zeros(Int, num_devices)
    for d in 1:num_devices
        for s in 1:payload_lens[d]
            sym_wave = lora_symbol(sf, bw, payloads[d][s])
            start_idx = Int(floor((tx_starts[d] + (s-1)*Tsym)*fs)) + 1
            end_idx = min(start_idx + length(sym_wave) - 1, length(rx_signal))
            if start_idx <= length(rx_signal) && end_idx >= start_idx
                ywin = rx_signal[start_idx:end_idx]
                m_hat = demod_lora(sf, bw, ywin)
                if m_hat != payloads[d][s]
                    errors[d] += 1
                end
            end
        end
    end

    return mean(errors ./ payload_lens)  # 端末平均SER
end

# ===== 複数端末数・SNRスイープ＋CSV出力 =====
function run_purealoha_sweep(sf::Int, bw::Float64,
    device_list::Vector{Int},
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(4,32),
    sim_time::Float64=1.0,
    iter::Int=1000,
    save_path::String="purealoha_ser.csv")
rng = MersenneTwister(1234)
snrs = collect(snr_min:snr_step:snr_max)
ser_mat = zeros(Float64, length(snrs), length(device_list)+1)
ser_mat[:,1] .= snrs

for (j, nd) in enumerate(device_list)
@info "Simulating num_devices = $nd"
for (i, snr) in enumerate(snrs)
ser_accum = 0.0
for it in 1:iter
ser_accum += purealoha_ser(sf, bw, nd;
                  payload_range=payload_range,
                  sim_time=sim_time,
                  snr_dB=snr,
                  rng=rng)
end
ser_mat[i, j+1] = ser_accum / iter
@printf("num_devices=%2d  SNR=%5.1f dB  SER=%.6f\n", nd, snr, ser_mat[i,j+1])
end
end

# CSV書き出し
header = ["SNR_dB"; ["devices=$(nd)" for nd in device_list]...]
open(save_path, "w") do io
println(io, join(header, ","))
writedlm(io, ser_mat, ',')
end
@info "結果をCSVに保存しました: $save_path"

return snrs, ser_mat
end


# ===== 実行例 =====
sf = 7
bw = 125e3
device_list = [2,8]
snr_min = -20.0
snr_max = 0.0
snr_step = 0.5
iter = 1000
snrs, ser_mat = run_purealoha_sweep(sf, bw, device_list,
                                    snr_min, snr_max, snr_step,
                                    payload_range=(4,32),
                                    sim_time=100.0,
                                    iter=iter,
                                    save_path="LoRa_SER/PureALOHA/results-PureALOHA/PureALOHA_iter$(iter)_sf$(sf)_dev$(join(device_list,'-')).csv")
