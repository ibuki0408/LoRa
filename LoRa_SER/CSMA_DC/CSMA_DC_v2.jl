using Random
using FFTW
using Statistics
using Printf
using DelimitedFiles

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

# 事前計算されたチャープ補償係数キャッシュ
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

# ===== CSMA + PER（1端末=1パケット、誤り1つでもPER=1）+ DC制約
function csma_per(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(32,64),
    sim_time::Float64=1.0,
    snr_dB::Float64=0.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    rng=Random.default_rng())

    M = 2^sf
    Tsym = M / bw
    fs = bw  # サンプル周波数（簡易）

    # 各端末の送信開始時間とペイロード長・データ
    tx_starts = rand(rng, num_devices) .* sim_time
    payload_lens = rand(rng, payload_range[1]:payload_range[2], num_devices)
    payloads = [rand(rng, 0:M-1, payload_lens[d]) for d in 1:num_devices]

    # 受信信号バッファとチャネル占有タイムライン
    total_samples = Int(ceil(sim_time*fs)) + maximum(payload_lens)*M
    rx_signal = zeros(ComplexF64, total_samples)
    channel_busy = falses(total_samples)

    # 事前計算されたシンボル波形キャッシュ
    symbol_cache = Dict{Int, Vector{ComplexF64}}()

    # 送信予定時刻順に端末をソート
    device_order = sortperm(tx_starts)

    # DC: 各端末の次回送信可能時刻
    next_time = zeros(Float64, num_devices)

    # 各端末の送信ルーチン（CSMA + DC）
    tx_start_indices = [Int[] for _ in 1:num_devices]
    for d_idx in device_order
        t = tx_starts[d_idx]

        # DC制約: 次回送信可能時刻まで待機
        if t < next_time[d_idx]
            t = next_time[d_idx]
        end

        # この端末のペイロード長（パケット長）
        L = payload_lens[d_idx]
        first_start_idx = 0

        # パケットの各シンボルをCSMAで順に送る
        for s in 1:L
            sym_idx = payloads[d_idx][s]
            if !haskey(symbol_cache, sym_idx)
                symbol_cache[sym_idx] = lora_symbol(sf, bw, sym_idx)
            end
            sym_wave = symbol_cache[sym_idx]
            sym_len = length(sym_wave)

            start_idx = Int(floor(t*fs)) + 1
            end_idx = start_idx + sym_len - 1

            # チャネル占有中はバックオフ
            while end_idx <= total_samples && any(channel_busy[start_idx:end_idx])
                t += rand(rng) * backoff_max
                start_idx = Int(floor(t*fs)) + 1
                end_idx = start_idx + sym_len - 1
            end

            if end_idx <= total_samples
                channel_busy[start_idx:end_idx] .= true
                @inbounds @simd for i in 1:sym_len
                    rx_signal[start_idx+i-1] += sym_wave[i]
                end
                if s == 1
                    first_start_idx = start_idx
                end
                push!(tx_start_indices[d_idx], start_idx)
            end

            t += Tsym
        end

        # DC: オフ時間を設定
        T_on = L * Tsym
        t_packet_end = (first_start_idx == 0) ? t : ((first_start_idx - 1) / fs + T_on)
        if dc > 0.0
            T_off = T_on * (1 - dc) / dc
            next_time[d_idx] = t_packet_end + T_off
        end
    end

    # AWGN付加
    add_awgn!(rx_signal, snr_dB)

    # 復調とパケットエラーフラグ
    packet_error = falses(num_devices)
    processed_counts = [length(tx_start_indices[d]) for d in 1:num_devices]
    for d in 1:num_devices
        # 1シンボルでも誤りがあれば、その端末のパケットはエラー
        for s in 1:processed_counts[d]
            sym_idx = payloads[d][s]
            sym_wave = symbol_cache[sym_idx]
            start_idx = tx_start_indices[d][s]
            end_idx = start_idx + length(sym_wave) - 1
            if start_idx <= total_samples && end_idx <= total_samples
                ywin = @view rx_signal[start_idx:end_idx]
                m_hat = demod_lora(sf, bw, ywin)
                if m_hat != sym_idx
                    packet_error[d] = true
                    break
                end
            end
        end
    end

    return mean(Float64.(packet_error))
end

# ===== SNRスイープ（PER vs SNR）＋CSV出力（num_devicesごとに列を分ける）
function run_csma_per_sweep(sf::Int, bw::Float64,
    device_list::Vector{Int},
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(32,64),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    iter::Int=200,
    seed::Int=1234,
    save_path::String="LoRa_SER/CSMA_DC/CSMA_PER_sweep.csv")

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    per_mat = zeros(Float64, length(snrs), length(device_list)+1)
    per_mat[:,1] .= snrs

    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end

    for (j, nd) in enumerate(device_list)
        @info "Simulating num_devices = $nd (PER)"
        for (i, snr) in enumerate(snrs)
            acc = 0.0
            for it in 1:iter
                acc += csma_per(sf, bw, nd;
                                payload_range=payload_range,
                                sim_time=sim_time,
                                snr_dB=snr,
                                backoff_max=backoff_max,
                                dc=dc,
                                rng=rng)
            end
            per = acc / iter
            per_mat[i, j+1] = per
            @printf("num_devices=%d  SNR=%6.1f dB  PER=%.6f\n", nd, snr, per)
        end
    end

    # CSV書き出し
    header = ["SNR_dB"; ["devices=$(nd)" for nd in device_list]...]
    open(save_path, "w") do io
        println(io, join(header, ","))
        writedlm(io, per_mat, ',')
    end
    @info "結果をCSVに保存しました: $save_path"

    return snrs, per_mat
end

# ===== 実行例（このまま実行可能） =====
sf = 7
bw = 125e3
device_list = [1000]                # 必要に応じて [2,4,8,16] などに変更
snr_min = -10.0
snr_max = 0.0
snr_step = 0.5
iter = 100                          # 反復回数（増やすと精度向上・時間増）
snrs, per_mat = run_csma_per_sweep(sf, bw, device_list,
                                   snr_min, snr_max, snr_step;
                                   payload_range=(32,64),
                                   sim_time=1.0,
                                   backoff_max=0.01,
                                   dc=0.01,
                                   iter=iter,
                                   seed=1234,
                                   save_path="LoRa_SER/CSMA_DC/results/CSMA_DC_PER_sf$(sf)_iter$(iter)_dev$(join(device_list,'-')).csv")