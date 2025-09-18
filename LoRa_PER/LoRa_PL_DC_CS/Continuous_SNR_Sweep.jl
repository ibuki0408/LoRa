using Random, FFTW, Statistics, Printf, DelimitedFiles, LinearAlgebra, StatsBase

# ===== パラメータ設定 =====
# エリアサイズ(km)
const area_size = 0.5

# パスロス係数
const α = 4.0
const β = 9.5
const γ = 4.5

# 搬送波周波数(MHz)
const f_c = 923.2

# シャドウイング標準偏差[dB]
const shadowing_std = 3.48

# 送信電力(dBm)
const Tx_dB = 13.0

# 雑音指数[dB]
const noise_figure = 10

# 帯域幅(Hz)
const band_width = 125e3

# 雑音スペクトラム密度(dBm/Hz)
const noise_power_spectrum_density = -174

# SF値（固定）
const SF = 7

# SNR閾値（SF7用）
const SNR_threshold = -7.5

# ===== LoRa 基本関数 =====
@inline function lora_symbol(sf::Int, bw::Float64, m::Int)
    M  = 2^sf
    tc = 1.0 / bw
    ts = M * tc
    x = Vector{ComplexF64}(undef, M)
    @inbounds for k in 1:M
        x[k] = exp(1im * 2π * (bw/(2*ts)) * ((k-1) + m)^2 * tc^2)
    end
    return x
end

const CHIRP_COMPENSATION_CACHE = Dict{Tuple{Int,Float64}, Vector{ComplexF64}}()

function get_chirp_compensation(sf::Int, bw::Float64)
    return get!(CHIRP_COMPENSATION_CACHE, (sf,bw)) do
        M = 2^sf
        tc = 1.0 / bw
        ts = tc * M
        exp.(-1im * 2π * (0:M-1).^2 * bw / (2*ts) * tc^2)
    end
end

# 復調
function demod_lora(sf::Int, bw::Float64, x::AbstractVector{ComplexF64})
    comp_coeff = get_chirp_compensation(sf, bw)
    d = x .* comp_coeff
    Y = fft(d)
    return argmax(abs.(Y)) - 1
end

# AWGN付加
function add_awgn!(y::AbstractVector{ComplexF64}, snr_dB::Float64)
    var = 10^(-snr_dB/10)
    σ = sqrt(var/2)
    @inbounds for i in eachindex(y)
        y[i] += ComplexF64(randn()*σ, randn()*σ)
    end
    return y
end

# ===== 伝搬モデル関数 =====
function distance(x1, y1, x2, y2)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

function path_loss(distance_km::Float64)
    return 10 * α * log10(distance_km) + β + 10 * γ * log10(f_c)
end

# シャドウイング値生成関数
function shadowing_value(rng)
    return randn(rng) * shadowing_std
end

# 受信電力計算関数
function received_power(tx_power_dBm::Float64, distance_km::Float64, shadowing_dB::Float64)
    pl = path_loss(distance_km)
    return tx_power_dBm - pl - shadowing_dB
end

# SNR計算関数
function calculate_snr(distance_km::Float64, shadowing_dB::Float64)
    received_pwr = received_power(Tx_dB, distance_km, shadowing_dB)
    snr = received_pwr - (noise_power_spectrum_density + 10*log10(band_width) + noise_figure)
    return snr
end

# dBm/mW 変換
@inline function dBm_to_mW(p_dBm::Float64)
    return 10.0^((p_dBm) / 10.0)
end

@inline function mW_to_dBm(p_mW::Float64)
    return p_mW > 0 ? 10.0 * log10(p_mW) : -Inf
end

# ===== ポアソン点過程による端末配置 =====
function generate_poisson_positions(num_devices::Int, rng)
    positions = Vector{Tuple{Float64,Float64}}(undef, num_devices)
    
    for i in 1:num_devices
        theta = rand(rng) * 2π
        r = sqrt(rand(rng)) * area_size
        x = r * cos(theta)
        y = r * sin(theta)
        positions[i] = (x, y)
    end
    
    return positions
end

# ===== 継続的送信版CSMA + DC + PER =====
function continuous_csma_dc_per(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=10.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    packet_interval::Float64=2.0,
    fixed_snr::Union{Float64, Nothing}=nothing,
    cs_threshold_dBm::Union{Float64, Nothing}=nothing,
    use_shadowing::Bool=true,
    rng=Random.default_rng())

    M  = 2^sf
    Ts = M / bw
    fs = bw

    # ポアソン点過程で端末配置生成
    node_positions = generate_poisson_positions(num_devices, rng)
    
    # シャドウイング値生成
    shadowing_values = use_shadowing ? [shadowing_value(rng) for _ in 1:num_devices] : zeros(num_devices)
    
    # 各端末のSNR計算（固定SNRまたは距離ベース）
    if fixed_snr !== nothing
        # 固定SNRの場合はシャドウイングを無効化して正確なSNRを保証
        device_snrs = fill(fixed_snr, num_devices)
        shadowing_values = zeros(num_devices)
    else
        if use_shadowing
            device_snrs = [calculate_snr(distance(pos[1], pos[2], 0.0, 0.0), shadowing_values[i]) for (i, pos) in enumerate(node_positions)]
        else
            device_snrs = [calculate_snr(distance(pos[1], pos[2], 0.0, 0.0)) for pos in node_positions]
        end
    end
    
    # 継続的送信用のデータ構造
    device_next_tx_time = rand(rng, num_devices) .* packet_interval
    device_packet_count = zeros(Int, num_devices)
    device_success_count = zeros(Int, num_devices)
    
    # 送信シンボル記録用のデータ構造
    device_payloads = Vector{Vector{Int}}[]
    device_packet_starts = Vector{Int}[]
    for d in 1:num_devices
        push!(device_payloads, Vector{Int}[])
        push!(device_packet_starts, Int[])
    end
    
    # 受信バッファ
    total_samples = Int(ceil(sim_time * fs * 2))
    rx_signal = zeros(ComplexF64, total_samples)
    channel_busy = falses(total_samples)
    
    # シンボルキャッシュ
    sym_cache = Dict{Int, Vector{ComplexF64}}()
    
    # 送信スケジュール記録
    scheduled_windows = Tuple{Int,Int,Int}[]
    
    # 受信電力ベースのキャリアセンス関数
    sensed_power_at_device = function(d::Int, n0::Int, n1::Int)
        if cs_threshold_dBm === nothing
            return -Inf
        end
        total_mW = 0.0
        xd, yd = node_positions[d]
        for (s, e, tx) in scheduled_windows
            if tx != d && !(e < n0 || s > n1)
                xt, yt = node_positions[tx]
                dist_km = distance(xd, yd, xt, yt)
                if use_shadowing
                    p_dBm = received_power(Tx_dB, dist_km, shadowing_values[tx])
                else
                    p_dBm = received_power(Tx_dB, dist_km)
                end
                total_mW += dBm_to_mW(p_dBm)
            end
        end
        return mW_to_dBm(total_mW)
    end
    
    # 継続的送信シミュレーション
    current_time = 0.0
    time_step = 0.001
    
    while current_time < sim_time
        for d in 1:num_devices
            if device_next_tx_time[d] <= current_time
                # パケット送信を試行
                payload_len = rand(rng, payload_range[1]:payload_range[2])
                payload = [rand(rng, 0:M-1) for _ in 1:payload_len]
                
                tx_start_time = current_time
                tx_success = true
                first_start_sample = 0
                
                # パケット内の各シンボルを送信
                for s in 1:payload_len
                    m = payload[s]
                    if !haskey(sym_cache, m)
                        sym_cache[m] = lora_symbol(sf, bw, m)
                    end
                    wave = sym_cache[m]
                    
                    n0 = Int(floor(tx_start_time * fs)) + 1
                    n1 = n0 + M - 1
                    
                    if n1 > total_samples
                        tx_success = false
                        break
                    end
                    
                    # CSMA: busyなら待機＋バックオフ
                    backoff_count = 0
                    max_backoffs = 10
                    
                    while backoff_count < max_backoffs && (
                        (cs_threshold_dBm === nothing && any(channel_busy[n0:n1])) ||
                        (cs_threshold_dBm !== nothing && sensed_power_at_device(d, n0, n1) >= cs_threshold_dBm)
                    )
                        backoff_time = rand(rng) * backoff_max
                        tx_start_time += backoff_time
                        n0 = Int(floor(tx_start_time * fs)) + 1
                        n1 = n0 + M - 1
                        backoff_count += 1
                        
                        if n1 > total_samples
                            tx_success = false
                            break
                        end
                    end
                    
                    if !tx_success || backoff_count >= max_backoffs
                        tx_success = false
                        break
                    end
                    
                    # 送信実行
                    if n0 >= 1 && n1 <= total_samples
                        channel_busy[n0:n1] .= true
                        
                        for i in 1:M
                            idx = n0 + i - 1
                            if idx >= 1 && idx <= length(rx_signal) && i <= length(wave)
                                rx_signal[idx] += wave[i]
                            end
                        end
                        
                        if s == 1
                            first_start_sample = n0
                            # 送信シンボルとパケット開始サンプルを記録
                            push!(device_payloads[d], copy(payload))
                            push!(device_packet_starts[d], n0)
                        end
                        
                        push!(scheduled_windows, (n0, n1, d))
                    end
                    
                    tx_start_time += Ts
                end
                
                # パケット送信完了
                device_packet_count[d] += 1
                
                if tx_success
                    device_success_count[d] += 1
                    
                    if dc > 0
                        T_on = payload_len * Ts
                        T_off = T_on * (1 - dc) / dc
                        device_next_tx_time[d] = current_time + T_on + T_off
                    else
                        device_next_tx_time[d] = current_time + packet_interval
                    end
                else
                    device_next_tx_time[d] = current_time + packet_interval * 0.1
                end
            end
        end
        
        current_time += time_step
    end
    
    # 各端末のSNRに基づくAWGN付加
    for d in 1:num_devices
        if device_packet_count[d] > 0
            snr = device_snrs[d]
            for (start_idx, end_idx, tx_device) in scheduled_windows
                if tx_device == d
                    if start_idx >= 1 && end_idx <= length(rx_signal)
                        signal_segment = @view rx_signal[start_idx:end_idx]
                        add_awgn!(signal_segment, snr)
                    end
                end
            end
        end
    end
    
    # 復調＋PER判定
    errors = falses(num_devices)
    device_packet_errors = zeros(Int, num_devices)
    
    for d in 1:num_devices
        if device_packet_count[d] > 0
            num_payloads = length(device_payloads[d])
            num_starts = length(device_packet_starts[d])
            num_packets = min(num_payloads, num_starts)
            
            for packet_idx in 1:num_packets
                payload = device_payloads[d][packet_idx]
                packet_start = device_packet_starts[d][packet_idx]
                packet_has_error = false
                
                # パケット内の各シンボルを復調
                for s in 1:length(payload)
                    start_idx = packet_start + (s-1) * M
                    end_idx = start_idx + M - 1
                    
                    if end_idx <= total_samples && start_idx >= 1 && end_idx <= length(rx_signal)
                        y = @view rx_signal[start_idx:end_idx]
                        m_hat = demod_lora(sf, bw, y)
                        
                        if m_hat != payload[s]
                            packet_has_error = true
                            break
                        end
                    else
                        packet_has_error = true
                        break
                    end
                end
                
                if packet_has_error
                    device_packet_errors[d] += 1
                end
            end
            
            if device_packet_errors[d] > 0
                errors[d] = true
            end
        end
    end
    
    # 統計計算
    total_packets = sum(device_packet_count)
    total_errors = sum(device_packet_errors)
    per = total_packets > 0 ? total_errors / total_packets : 0.0
    
    return per, node_positions, device_snrs, errors, shadowing_values, device_packet_count, device_success_count
end

# ===== 継続的送信版SNRスイープ機能 =====
function run_continuous_snr_sweep_per(sf::Int, bw::Float64, num_devices::Int,
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=10.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    packet_interval::Float64=2.0,
    iter::Int=100,
    seed::Int=1234,
    save_path::String="results_LoRa_simple_model/continuous_snr_sweep_per.csv",
    cs_threshold_dBm::Union{Float64,Nothing}=nothing,
    use_shadowing::Bool=true)

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    per_vals = zeros(Float64, length(snrs))
    
    println("継続的送信SNRスイープ実行中...")
    println("SNR範囲: $(snr_min) ~ $(snr_max) dB, ステップ: $(snr_step) dB")
    println("反復回数: $iter")
    println("シミュレーション時間: $sim_time 秒")
    println("パケット送信間隔: $packet_interval 秒")
    
    for (i, snr) in enumerate(snrs)
        acc_per = 0.0
        
        for it in 1:iter
            local_rng = MersenneTwister(seed + it + i * 1000)
            
            per, positions, snrs_actual, errors, shadowing_vals, packet_counts, success_counts = continuous_csma_dc_per(sf, bw, num_devices;
                                               payload_range=payload_range,
                                               sim_time=sim_time,
                                               backoff_max=backoff_max,
                                               dc=dc,
                                               packet_interval=packet_interval,
                                               fixed_snr=snr,
                                               cs_threshold_dBm=cs_threshold_dBm,
                                               use_shadowing=use_shadowing,
                                               rng=local_rng)
            
            acc_per += per
        end
        
        per_vals[i] = acc_per / iter
        
        @printf("SNR = %5.1f dB, PER = %.6f\n", snr, per_vals[i])
    end

    # CSV保存
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    
    open(save_path, "w") do io
        println(io, "SNR_dB,PER")
        for i in 1:length(snrs)
            println(io, "$(snrs[i]),$(per_vals[i])")
        end
    end
    
    @info "継続的送信SNRスイープ結果をCSVに保存しました: $save_path"
    
    return snrs, per_vals
end

# ===== メイン実行部分 =====
println("継続的送信SNRスイープ専用プログラム")
println("設定:")
println("  SF: $(SF)")
println("  エリアサイズ: $(area_size) km")
println("  送信電力: $(Tx_dB) dBm")
println("  搬送波周波数: $(f_c) MHz")
println("  パスロス係数: α=$(α), β=$(β), γ=$(γ)")
println("  シャドウイング標準偏差: $(shadowing_std) dB")

# パラメータ設定
sf = SF
bw = 125e3
num_devices = 10

# SNRスイープ実行
println("\n=== 継続的送信SNRスイープ実行 ===")
snr_min, snr_max, snr_step = -20.0, 0.0, 0.5
iter_sweep = 1000

snrs_cont, per_vals_cont = run_continuous_snr_sweep_per(sf, bw, num_devices,
                               snr_min, snr_max, snr_step;
                               payload_range=(16,32),
                               sim_time=10.0,
                               backoff_max=0.01,
                               dc=0.01,
                               packet_interval=2.0,
                               iter=iter_sweep,
                               cs_threshold_dBm=-110.0,
                               use_shadowing=true,
                               save_path="LoRa_PER/LoRa_PL_DC_CS/results_CSMA/continuous_snr_sweep_per_sf$(sf)_dev$(num_devices)_iter$(iter_sweep).csv")

# 結果の表示
println("\n=== 継続的送信SNRスイープ結果 ===")
println("SNR範囲: $(snr_min) ~ $(snr_max) dB")
println("最小PER: $(minimum(per_vals_cont))")
println("最大PER: $(maximum(per_vals_cont))")

println("\nシミュレーション完了！")
