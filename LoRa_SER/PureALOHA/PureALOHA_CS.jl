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

# ===== CSMA版シンボル送信（チャネル占有タイムライン方式） =====
function csma_ser(sf::Int, bw::Float64, num_devices::Int;
                  payload_range::Tuple{Int,Int}=(4,32),
                  sim_time::Float64=1.0,
                  snr_dB::Float64=0.0,
                  backoff_max::Float64=0.5,
                  rng=Random.default_rng())

    M = 2^sf
    Tsym = M / bw
    fs = bw  # サンプル周波数（簡易）

    # 各端末の送信開始時間とペイロード長
    tx_starts = rand(rng, num_devices) .* sim_time
    payload_lens = rand(rng, payload_range[1]:payload_range[2], num_devices)
    payloads = [rand(rng, 0:M-1, payload_lens[d]) for d in 1:num_devices]

    # 受信信号バッファとチャネル占有タイムライン
    total_samples = Int(ceil(sim_time*fs)) + maximum(payload_lens)*M
    rx_signal = zeros(ComplexF64, total_samples)
    channel_busy = falses(total_samples)

    # 事前計算されたシンボル波形キャッシュ
    symbol_cache = Dict{Int, Vector{ComplexF64}}()

    # 各端末の送信ループ
    # 実際の各シンボル開始インデックスを記録（デコード時に使用）
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

            # 送信前にチャネル占有タイムラインを確認
            start_idx = Int(floor(t*fs)) + 1
            end_idx = start_idx + sym_len - 1

            while end_idx <= total_samples && any(channel_busy[start_idx:end_idx])
                # チャネル占有 → バックオフ
                t += rand(rng) * backoff_max
                start_idx = Int(floor(t*fs)) + 1
                end_idx = start_idx + sym_len - 1
            end

            # 送信確定 → チャネルを占有状態にする
            if end_idx <= total_samples
                channel_busy[start_idx:end_idx] .= true
                @inbounds @simd for i in 1:sym_len
                    rx_signal[start_idx+i-1] += sym_wave[i]
                end
                # 実際の開始位置を記録
                push!(tx_start_indices[d], start_idx)
            end

            t += Tsym
        end
    end

    # AWGN付加
    add_awgn!(rx_signal, snr_dB)

    # 復調と誤り判定
    errors = zeros(Int, num_devices)
    # 実際に送信できたシンボル数（バッファ外に出た分は除外）
    processed_counts = [length(tx_start_indices[d]) for d in 1:num_devices]
    for d in 1:num_devices
        for s in 1:processed_counts[d]
            sym_idx = payloads[d][s]
            sym_wave = symbol_cache[sym_idx]
            # 送信時に記録した実開始位置を使用
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

    # 端末ごとの実処理シンボル数で割る（ゼロ回避）
    denom = [max(1, c) for c in processed_counts]
    return mean(errors ./ denom)
end

# ===== 複数端末数・SNRスイープ＋CSV出力 =====
function run_csma_sweep(sf::Int, bw::Float64,
    device_list::Vector{Int},
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(4,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.5,
    iter::Int=1000,
    save_path::String="results/CSMA_SER.csv")

    rng = MersenneTwister(1234)
    snrs = collect(snr_min:snr_step:snr_max)
    ser_mat = zeros(Float64, length(snrs), length(device_list)+1)
    ser_mat[:,1] .= snrs

    # ディレクトリがなければ作成
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end

    for (j, nd) in enumerate(device_list)
        @info "Simulating num_devices = $nd"
        for (i, snr) in enumerate(snrs)
            ser_accum = 0.0
            for it in 1:iter
                ser_accum += csma_ser(sf, bw, nd;
                                      payload_range=payload_range,
                                      sim_time=sim_time,
                                      snr_dB=snr,
                                      backoff_max=backoff_max,
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
device_list = [1000]
snr_min = -10.0
snr_max = 0.0
snr_step = 0.5
iter = 1000
snrs, ser_mat = run_csma_sweep(sf, bw, device_list,
                               snr_min, snr_max, snr_step,
                               payload_range=(4,32),
                               sim_time=1.0,
                               backoff_max=0.01,
                               iter=iter,
                               save_path="LoRa_SER/PureALOHA/results-CS/CSMA_SER_sf$(sf)_iter$(iter)_dev$(join(device_list,'-')).csv")
