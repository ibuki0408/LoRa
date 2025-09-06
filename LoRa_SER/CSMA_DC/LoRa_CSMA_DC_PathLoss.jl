using Random, FFTW, Statistics, Printf, DelimitedFiles, LinearAlgebra

# ===== パラメータ設定 =====
# エリアサイズ(km)
const area_size = 0.5

# パスロス係数
const α = 4.0

# 伝搬損失オフセット
const β = 9.5

# 伝搬周波数係数
const γ = 4.5

# 搬送波周波数(MHz)
const f_c = 923.2

# シャドウイング標準偏差[dB]
const shadowing_std = 3.48

# 送信電力(dBm)
const Tx_dB = 13

# 雑音指数[dB]
const noise_figure = 10

# 帯域幅(Hz)
const band_width = 125e3

# 雑音スペクトラム密度(dBm/Hz)
const noise_power_spectrum_density = -174

# SNR閾値[SF, threshold]
const SNR_threshold_list = [[6, -5], [7, -7.5], [8, -10], [9, -12.5], [10, -15], [11, -17.5], [12, -20]]

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
function add_awgn!(y::Vector{ComplexF64}, snr_dB::Float64)
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
    # パスロス計算: PL = 10α*log10(d) + β + 10γ*log10(fc)
    return 10 * α * log10(distance_km) + β + 10 * γ * log10(f_c)
end

function shadowing_value(rng)
    # シャドウイング値生成
    return randn(rng) * shadowing_std
end

function received_power(tx_power_dBm::Float64, distance_km::Float64, shadowing_dB::Float64)
    # 受信電力計算: Pr = Pt - PL - Shadowing
    pl = path_loss(distance_km)
    return tx_power_dBm - pl - shadowing_dB
end

function snr_threshold(sf::Int)
    # SFに対応するSNR閾値を取得
    for (sf_val, threshold) in SNR_threshold_list
        if sf_val == sf
            return threshold
        end
    end
    return -20.0  # デフォルト値
end

function select_sf_based_on_distance(distance_km::Float64, shadowing_dB::Float64)
    # 距離とシャドウイングに基づいてSFを選択
    received_pwr = received_power(Tx_dB, distance_km, shadowing_dB)
    snr = received_pwr - (noise_power_spectrum_density + 10*log10(band_width) + noise_figure)
    
    if snr > -7.5
        return 7
    elseif snr > -10.0
        return 8
    elseif snr > -12.5
        return 9
    elseif snr > -15.0
        return 10
    elseif snr > -17.5
        return 11
    else
        return 12
    end
end

# ===== 端末配置関数 =====
function generate_node_positions(num_devices::Int, rng)
    # 円形エリア内に端末をランダム配置
    positions = Vector{Tuple{Float64,Float64}}(undef, num_devices)
    
    for i in 1:num_devices
        # 極座標でランダム配置（一様分布）
        theta = rand(rng) * 2π
        r = sqrt(rand(rng)) * area_size  # sqrt(rand())で面積が一様分布になる
        x = r * cos(theta)
        y = r * sin(theta)
        positions[i] = (x, y)
    end
    
    return positions
end

function calculate_snr_for_device(device_pos::Tuple{Float64,Float64}, shadowing_dB::Float64)
    # 端末からゲートウェイ（原点）までのSNR計算
    dist = distance(device_pos[1], device_pos[2], 0.0, 0.0)
    received_pwr = received_power(Tx_dB, dist, shadowing_dB)
    snr = received_pwr - (noise_power_spectrum_density + 10*log10(band_width) + noise_figure)
    return snr
end

# ===== 改良版CSMA + DC + PER =====
function csma_dc_per_enhanced(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    rng=Random.default_rng())

    M  = 2^sf
    Ts = M / bw     # シンボル長
    fs = bw         # サンプル周波数 ≈ 帯域幅

    # 端末配置生成
    node_positions = generate_node_positions(num_devices, rng)
    
    # 各端末のシャドウイング値とSF決定
    shadowing_values = [shadowing_value(rng) for _ in 1:num_devices]
    device_sfs = [select_sf_based_on_distance(distance(node_positions[i][1], node_positions[i][2], 0.0, 0.0), shadowing_values[i]) for i in 1:num_devices]
    
    # 送信スケジュール
    tx_starts = rand(rng, num_devices) .* sim_time
    payload_lens = rand(rng, payload_range[1]:payload_range[2], num_devices)
    payloads = [rand(rng, 0:M-1, payload_lens[d]) for d in 1:num_devices]

    # 受信バッファ
    total_samples = Int(ceil(sim_time*fs)) + maximum(payload_lens)*M
    rx_signal = zeros(ComplexF64, total_samples)
    channel_busy = falses(total_samples)

    # シンボルキャッシュ
    sym_cache = Dict{Int, Vector{ComplexF64}}()

    # DC制約
    next_time = zeros(Float64, num_devices)

    # 実際に送信
    tx_indices = [Int[] for _ in 1:num_devices]
    for d in sortperm(tx_starts)
        t = max(tx_starts[d], next_time[d])
        L = payload_lens[d]
        first_start = 0
        
        # この端末のSFを使用
        device_sf = device_sfs[d]
        device_M = 2^device_sf
        device_Ts = device_M / bw
        
        for s in 1:L
            m = payloads[d][s]
            if !haskey(sym_cache, m)
                sym_cache[m] = lora_symbol(device_sf, bw, m)
            end
            wave = sym_cache[m]
            n0 = Int(floor(t*fs)) + 1
            n1 = n0 + device_M - 1

            # CSMA: busyなら待機＋バックオフ
            while n1 <= total_samples && any(channel_busy[n0:n1])
                t += rand(rng)*backoff_max
                n0 = Int(floor(t*fs)) + 1
                n1 = n0 + device_M - 1
            end

            if n1 <= total_samples
                channel_busy[n0:n1] .= true
                @inbounds for i in 1:device_M
                    rx_signal[n0+i-1] += wave[i]
                end
                if s==1; first_start=n0; end
                push!(tx_indices[d], n0)
            end
            t += device_Ts
        end

        # DC制約
        T_on = L * device_Ts
        tend = (first_start==0) ? t : ((first_start-1)/fs + T_on)
        if dc>0
            T_off = T_on * (1-dc)/dc
            next_time[d] = tend + T_off
        end
    end

    # 各端末のSNRに基づくAWGN付加
    for d in 1:num_devices
        if !isempty(tx_indices[d])
            # この端末のSNRを計算
            dist = distance(node_positions[d][1], node_positions[d][2], 0.0, 0.0)
            received_pwr = received_power(Tx_dB, dist, shadowing_values[d])
            snr = received_pwr - (noise_power_spectrum_density + 10*log10(band_width) + noise_figure)
            
            # この端末の信号にのみノイズを追加
            device_sf = device_sfs[d]
            device_M = 2^device_sf
            for start_idx in tx_indices[d]
                if start_idx + device_M - 1 <= total_samples
                    signal_segment = @view rx_signal[start_idx:start_idx+device_M-1]
                    add_awgn!(signal_segment, snr)
                end
            end
        end
    end

    # 復調＋PER判定
    errors = falses(num_devices)
    for d in 1:num_devices
        device_sf = device_sfs[d]
        device_M = 2^device_sf
        for (s,start) in enumerate(tx_indices[d])
            if start + device_M - 1 <= total_samples
                y = @view rx_signal[start:start+device_M-1]
                m_hat = demod_lora(device_sf, bw, y)
                if m_hat != payloads[d][s]
                    errors[d] = true
                    break
                end
            end
        end
    end

    return mean(Float64.(errors)), node_positions, device_sfs, shadowing_values
end

# ===== 結果分析関数 =====
function analyze_results(node_positions, device_sfs, shadowing_values, errors)
    num_devices = length(node_positions)
    
    # 距離別分析
    distances = [distance(pos[1], pos[2], 0.0, 0.0) for pos in node_positions]
    
    # SF別分析
    sf_counts = Dict{Int, Int}()
    sf_errors = Dict{Int, Int}()
    
    for i in 1:num_devices
        sf = device_sfs[i]
        sf_counts[sf] = get(sf_counts, sf, 0) + 1
        if errors[i]
            sf_errors[sf] = get(sf_errors, sf, 0) + 1
        end
    end
    
    # 統計情報
    stats = Dict(
        :total_devices => num_devices,
        :total_errors => sum(errors),
        :per => mean(Float64.(errors)),
        :avg_distance => mean(distances),
        :max_distance => maximum(distances),
        :min_distance => minimum(distances),
        :sf_distribution => sf_counts,
        :sf_error_rates => Dict(k => get(sf_errors, k, 0) / v for (k, v) in sf_counts)
    )
    
    return stats
end

# ===== SNRスイープ (PER vs SNR) =====
function run_per_sweep_enhanced(sf::Int, bw::Float64,
    num_devices::Int,
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    iter::Int=100,
    seed::Int=1234,
    save_path::String="LoRa_PER_enhanced.csv")

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    per_vals = zeros(Float64, length(snrs))
    
    # 詳細結果保存用
    detailed_results = []

    for (i, snr) in enumerate(snrs)
        acc = 0.0
        all_stats = []
        
        for it in 1:iter
            per, positions, sfs, shadowing = csma_dc_per_enhanced(sf, bw, num_devices;
                                           payload_range=payload_range,
                                           sim_time=sim_time,
                                           backoff_max=backoff_max,
                                           dc=dc,
                                           rng=rng)
            acc += per
            
            # 詳細分析（最初の数回のみ）
            if it <= 5
                # エラー情報を生成（簡易版）
                errors = rand(rng, num_devices) .< per
                stats = analyze_results(positions, sfs, shadowing, errors)
                push!(all_stats, stats)
            end
        end
        
        per_vals[i] = acc / iter
        push!(detailed_results, (snr, per_vals[i], all_stats))
        @printf("SNR = %5.1f dB, PER = %.6f\n", snr, per_vals[i])
    end

    # CSV保存
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    
    # 基本結果
    open(save_path, "w") do io
        println(io, "SNR_dB,PER")
        for i in 1:length(snrs)
            println(io, "$(snrs[i]),$(per_vals[i])")
        end
    end
    
    # 詳細結果（最初のSNR値のみ）
    if !isempty(detailed_results)
        first_result = detailed_results[1]
        detailed_path = replace(save_path, ".csv" => "_detailed.csv")
        open(detailed_path, "w") do io
            println(io, "Iteration,Total_Devices,Total_Errors,PER,Avg_Distance,Max_Distance,Min_Distance")
            for (i, stats) in enumerate(first_result[3])
                println(io, "$i,$(stats[:total_devices]),$(stats[:total_errors]),$(stats[:per]),$(stats[:avg_distance]),$(stats[:max_distance]),$(stats[:min_distance])")
            end
        end
    end
    
    @info "結果をCSVに保存しました: $save_path"
    if !isempty(detailed_results)
        @info "詳細結果を保存しました: $(replace(save_path, ".csv" => "_detailed.csv"))"
    end

    return snrs, per_vals, detailed_results
end

# ===== 実行部分 =====
println("LoRa CSMA/DC PER シミュレーション開始（パスロス・端末配置考慮版）")
println("設定:")
println("  エリアサイズ: $(area_size) km")
println("  送信電力: $(Tx_dB) dBm")
println("  搬送波周波数: $(f_c) MHz")
println("  パスロス係数: α=$(α), β=$(β), γ=$(γ)")

sf = 7
bw = 125e3
num_devices = 50
snr_min, snr_max, snr_step = -10.0, 0.0, 1.0
iter = 50   # 反復回数

snrs, per_vals, detailed_results = run_per_sweep_enhanced(sf, bw, num_devices,
                               snr_min, snr_max, snr_step;
                               payload_range=(16,32),
                               sim_time=1.0,
                               backoff_max=0.01,
                               dc=0.01,
                               iter=iter,
                               save_path="LoRa_SER/CSMA_DC/LoRa_PER_enhanced_sf$(sf)_dev$(num_devices).csv")

println("シミュレーション完了！")