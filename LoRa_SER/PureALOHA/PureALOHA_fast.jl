using Random
using FFTW
using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra

# ===== LoRa基本関数 =====
@inline function lora_symbol(sf::Int, bw::Float64, m::Int)
    M  = 2^sf
    tc = 1.0 / bw
    ts = tc * M
    x = Vector{ComplexF64}(undef, M)
    @inbounds @simd for k in 1:M
        x[k] = exp(1im * 2π * (bw/(2*ts)) * ((k-1) + m)^2 * tc^2)
    end
    return x
end

# 事前計算されたチャープ補償係数
const CHIRP_COMPENSATION_CACHE = Dict{Tuple{Int,Float64}, Vector{ComplexF64}}()
@inline function get_chirp_compensation(sf::Int, bw::Float64)
    return get!(CHIRP_COMPENSATION_CACHE, (sf, bw)) do
        M = 2^sf
        tc = 1.0 / bw
        ts = tc * M
        exp.(-1im * 2π * (bw/(2*ts)) * ((0:M-1).^2) * tc^2)
    end
end

@inline function demod_lora(sf::Int, bw::Float64, x::AbstractVector{ComplexF64})
    comp_coeff = get_chirp_compensation(sf, bw)
    d = x .* comp_coeff
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

# ===== CSなし版シンボル送信（Pure ALOHA的） =====
function aloha_ser(sf::Int, bw::Float64, num_devices::Int;
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

    # 受信信号バッファ
    total_samples = Int(ceil(sim_time*fs)) + maximum(payload_lens)*M
    rx_signal = zeros(ComplexF64, total_samples)

    # 事前計算されたシンボル波形キャッシュ
    symbol_cache = Dict{Int, Vector{ComplexF64}}()

    # 各端末の送信ループ
    tx_start_indices = [Int[] for _ in 1:num_devices]
    for d in 1:num_devices
        t = tx_starts[d]
        for s in 1:payload_lens[d]
            sym_idx = payloads[d][s]
            if !haskey(symbol_cache, sym_idx)
                symbol_cache[sym_idx] = lora_symbol(sf, bw, sym_idx)
            end
            sym_wave = symbol_cache[sym_idx]
            sym_len = length(sym_wave)

            # 送信（CSチェックなし、重なってもそのまま送る）
            start_idx = Int(floor(t*fs)) + 1
            end_idx = start_idx + sym_len - 1
            if end_idx <= total_samples
                @inbounds @simd for i in 1:sym_len
                    rx_signal[start_idx+i-1] += sym_wave[i]
                end
                push!(tx_start_indices[d], start_idx)
            end

            t += Tsym
        end
    end

    # AWGN付加
    add_awgn!(rx_signal, snr_dB)

    # 復調と誤り判定
    errors = zeros(Int, num_devices)
    processed_counts = [length(tx_start_indices[d]) for d in 1:num_devices]
    for d in 1:num_devices
        for s in 1:processed_counts[d]
            sym_idx = payloads[d][s]
            sym_wave = symbol_cache[sym_idx]
            start_idx = tx_start_indices[d][s]
            end_idx = start_idx + length(sym_wave) - 1
            if start_idx <= total_samples && end_idx <= total_samples
                ywin = @view rx_signal[start_idx:end_idx]
                m_hat = demod_lora(sf, bw, ywin)
                if m_hat != sym_idx
                    errors[d] += 1
                end
            end
        end
    end

    denom = [max(1, c) for c in processed_counts]
    return mean(errors ./ denom)
end

# ===== 複数端末数・SNRスイープ＋CSV出力 =====
function run_aloha_sweep(sf::Int, bw::Float64,
    device_list::Vector{Int},
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(4,32),
    sim_time::Float64=1.0,
    iter::Int=1000,
    save_path::String="results/ALOHA_SER.csv")

    rng = MersenneTwister(1234)
    snrs = collect(snr_min:snr_step:snr_max)
    ser_mat = zeros(Float64, length(snrs), length(device_list)+1)
    ser_mat[:,1] .= snrs

    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end

    for (j, nd) in enumerate(device_list)
        @info "Simulating num_devices = $nd"
        for (i, snr) in enumerate(snrs)
            ser_accum = 0.0
            for it in 1:iter
                ser_accum += aloha_ser(sf, bw, nd;
                                       payload_range=payload_range,
                                       sim_time=sim_time,
                                       snr_dB=snr,
                                       rng=rng)
            end
            ser_mat[i, j+1] = ser_accum / iter
            @printf("ALOHA: num_devices=%2d  SNR=%5.1f dB  SER=%.6f\n", nd, snr, ser_mat[i,j+1])
        end
    end

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
iter = 10000
snrs, ser_mat = run_aloha_sweep(sf, bw, device_list,
                                snr_min, snr_max, snr_step,
                                payload_range=(4,32),
                                sim_time=1.0,
                                iter=iter,
                                save_path="LoRa_SER/PureALOHA/results-ALOHA/PureALOHA_SER_sf$(sf)_iter$(iter)_dev$(join(device_list,'-')).csv")
