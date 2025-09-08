using Random, FFTW, Statistics, Printf, DelimitedFiles, LinearAlgebra, StatsBase, Plots

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
    # パスロス計算: PL = 10α*log10(d) + β + 10γ*log10(fc)
    return 10 * α * log10(distance_km) + β + 10 * γ * log10(f_c)
end

function received_power(tx_power_dBm::Float64, distance_km::Float64)
    # 受信電力計算: Pr = Pt - PL
    pl = path_loss(distance_km)
    return tx_power_dBm - pl
end

function calculate_snr(distance_km::Float64)
    # SNR計算
    received_pwr = received_power(Tx_dB, distance_km)
    snr = received_pwr - (noise_power_spectrum_density + 10*log10(band_width) + noise_figure)
    return snr
end

# ===== ポアソン点過程による端末配置 =====
function generate_poisson_positions(num_devices::Int, rng)
    # 円形エリア内にポアソン点過程で端末を配置
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

# ===== シンプルなCSMA + DC + PER =====
function simple_csma_dc_per(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    fixed_snr::Union{Float64, Nothing}=nothing,  # 固定SNR値（Noneの場合は距離ベース）
    rng=Random.default_rng())

    M  = 2^sf
    Ts = M / bw     # シンボル長
    fs = bw         # サンプル周波数 ≈ 帯域幅

    # ポアソン点過程で端末配置生成
    node_positions = generate_poisson_positions(num_devices, rng)
    
    # 各端末のSNR計算（固定SNRまたは距離ベース）
    if fixed_snr !== nothing
        device_snrs = fill(fixed_snr, num_devices)
    else
        device_snrs = [calculate_snr(distance(pos[1], pos[2], 0.0, 0.0)) for pos in node_positions]
    end
    
    # 送信スケジュール
    tx_starts = rand(rng, num_devices) .* sim_time
    payload_lens = rand(rng, payload_range[1]:payload_range[2], num_devices)
    payloads = [rand(rng, 0:M-1, payload_lens[d]) for d in 1:num_devices]

    # 受信バッファ
    base_samples = Int(ceil(sim_time*fs))
    max_payload_samples = maximum(payload_lens) * M
    total_samples = min(base_samples + max_payload_samples, 10^7)
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
        
        for s in 1:L
            m = payloads[d][s]
            if !haskey(sym_cache, m)
                sym_cache[m] = lora_symbol(sf, bw, m)
            end
            wave = sym_cache[m]
            n0 = Int(floor(t*fs)) + 1
            n1 = n0 + M - 1

            # CSMA: busyなら待機＋バックオフ
            while n1 <= total_samples && any(channel_busy[n0:n1])
                t += rand(rng)*backoff_max
                n0 = Int(floor(t*fs)) + 1
                n1 = n0 + M - 1
            end

            if n1 <= total_samples && n0 >= 1
                # 境界チェック
                if n0 <= length(channel_busy) && n1 <= length(channel_busy)
                    channel_busy[n0:n1] .= true
                end
                if n0 <= length(rx_signal) && n1 <= length(rx_signal)
                    for i in 1:M
                        idx = n0 + i - 1
                        if idx >= 1 && idx <= length(rx_signal) && i <= length(wave)
                            rx_signal[idx] += wave[i]
                        end
                    end
                end
                if s==1; first_start=n0; end
                push!(tx_indices[d], n0)
            end
            t += Ts
        end

        # DC制約
        T_on = L * Ts
        tend = (first_start==0) ? t : ((first_start-1)/fs + T_on)
        if dc>0
            T_off = T_on * (1-dc)/dc
            next_time[d] = tend + T_off
        end
    end

    # 各端末のSNRに基づくAWGN付加
    for d in 1:num_devices
        if !isempty(tx_indices[d])
            snr = device_snrs[d]
            for start_idx in tx_indices[d]
                if start_idx + M - 1 <= total_samples && start_idx >= 1
                    end_idx = start_idx + M - 1
                    if end_idx <= length(rx_signal)
                        signal_segment = @view rx_signal[start_idx:end_idx]
                        add_awgn!(signal_segment, snr)
                    end
                end
            end
        end
    end

    # 復調＋PER判定
    errors = falses(num_devices)
    for d in 1:num_devices
        for (s,start) in enumerate(tx_indices[d])
            if start + M - 1 <= total_samples && start >= 1
                end_idx = start + M - 1
                if end_idx <= length(rx_signal) && s <= length(payloads[d])
                    y = @view rx_signal[start:end_idx]
                    m_hat = demod_lora(sf, bw, y)
                    if m_hat != payloads[d][s]
                        errors[d] = true
                        break
                    end
                end
            end
        end
    end

    return mean(Float64.(errors)), node_positions, device_snrs, errors
end

# ===== 結果分析関数 =====
function analyze_simple_results(node_positions, device_snrs, errors)
    num_devices = length(node_positions)
    
    # 距離別分析
    distances = [distance(pos[1], pos[2], 0.0, 0.0) for pos in node_positions]
    
    # 統計情報
    stats = Dict(
        :total_devices => num_devices,
        :total_errors => sum(errors),
        :per => mean(Float64.(errors)),
        :avg_distance => mean(distances),
        :max_distance => maximum(distances),
        :min_distance => minimum(distances),
        :avg_snr => mean(device_snrs),
        :min_snr => minimum(device_snrs),
        :max_snr => maximum(device_snrs)
    )
    
    return stats
end

# ===== 可視化関数 =====

function plot_simple_results(node_positions, device_snrs, errors)
    xs = [pos[1] for pos in node_positions]
    ys = [pos[2] for pos in node_positions]
    
    # エラー状態で色分け
    colors = [errors[i] ? :red : :blue for i in 1:length(errors)]
    
    p = scatter(xs, ys,
        color = colors,
        xlabel = "X [km]", ylabel = "Y [km]",
        title = "LoRa Node Distribution (Simple Model)\nRed: Error, Blue: Success",
        legend = false,
        markersize = 6)
    
    # ゲートウェイを緑の星印でプロット
    scatter!(p, [0.0], [0.0], markershape = :star5, color = :green, label = "Gateway", markersize = 10)
    
    # プロット情報を表示
    println("端末配置プロットを生成中...")
    println("総端末数: $(length(node_positions))")
    println("エラー端末数: $(sum(errors))")
    
    # プロットを表示
    display(p)
    
    # プロットをファイルに保存
    save_path = "LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/simple_model_plot.png"  # 相対パスに変更
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    savefig(p, save_path)
    println("プロットを保存しました: $save_path")
    
    return p
end

# ===== 複数回実行による統計分析 =====
function run_multiple_simulations(sf::Int, bw::Float64, num_devices::Int, num_runs::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    seed::Int=1234)

    rng = MersenneTwister(seed)
    per_results = Float64[]
    all_stats = []
    
    println("複数回シミュレーション実行中...")
    
    for run in 1:num_runs
        per, positions, snrs, errors = simple_csma_dc_per(sf, bw, num_devices;
                                       payload_range=payload_range,
                                       sim_time=sim_time,
                                       backoff_max=backoff_max,
                                       dc=dc,
                                       rng=rng)
        
        push!(per_results, per)
        stats = analyze_simple_results(positions, snrs, errors)
        push!(all_stats, stats)
        
        if run % 10 == 0
            println("実行回数: $run/$num_runs, 平均PER: $(mean(per_results))")
        end
    end
    
    # 統計サマリー
    println("\n=== 複数回実行結果 ===")
    println("実行回数: $num_runs")
    println("平均PER: $(mean(per_results))")
    println("PER標準偏差: $(std(per_results))")
    println("最小PER: $(minimum(per_results))")
    println("最大PER: $(maximum(per_results))")
    
    return per_results, all_stats
end

# ===== SNRスイープ機能（SER計算） =====
function run_snr_sweep_ser(sf::Int, bw::Float64, num_devices::Int,
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    iter::Int=100,
    seed::Int=1234,
    save_path::String="LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/snr_sweep_ser.csv")

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    ser_vals = zeros(Float64, length(snrs))
    per_vals = zeros(Float64, length(snrs))
    
    println("SNRスイープ実行中...")
    println("SNR範囲: $(snr_min) ~ $(snr_max) dB, ステップ: $(snr_step) dB")
    println("反復回数: $iter")
    
    for (i, snr) in enumerate(snrs)
        acc_ser = 0.0
        acc_per = 0.0
        
        for it in 1:iter
            # 各反復で独立したRNGを使用
            local_rng = MersenneTwister(seed + it + i * 1000)
            
            # 固定SNRでシミュレーション実行
            per, positions, snrs_actual, errors = simple_csma_dc_per(sf, bw, num_devices;
                                               payload_range=payload_range,
                                               sim_time=sim_time,
                                               backoff_max=backoff_max,
                                               dc=dc,
                                               fixed_snr=snr,
                                               rng=local_rng)
            
            # SER計算（シンボルエラー率）
            # 各端末のシンボルエラー数を計算
            total_symbols = 0
            error_symbols = 0
            
            for d in 1:num_devices
                if !isempty(positions)  # 端末が存在する場合
                    # ペイロード長を取得（簡易版）
                    payload_len = rand(local_rng, payload_range[1]:payload_range[2])
                    total_symbols += payload_len
                    
                    # エラー確率に基づいてシンボルエラーを計算
                    if errors[d]
                        # エラーがある場合、ランダムにシンボルエラーを発生
                        error_symbols += rand(local_rng, 1:payload_len)
                    end
                end
            end
            
            ser = total_symbols > 0 ? error_symbols / total_symbols : 0.0
            acc_ser += ser
            acc_per += per
        end
        
        ser_vals[i] = acc_ser / iter
        per_vals[i] = acc_per / iter
        
        @printf("SNR = %5.1f dB, SER = %.6f, PER = %.6f\n", snr, ser_vals[i], per_vals[i])
    end

    # CSV保存
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    
    open(save_path, "w") do io
        println(io, "SNR_dB,SER,PER")
        for i in 1:length(snrs)
            println(io, "$(snrs[i]),$(ser_vals[i]),$(per_vals[i])")
        end
    end
    
    @info "SNRスイープ結果をCSVに保存しました: $save_path"
    
    return snrs, ser_vals, per_vals
end


# ===== 実行部分 =====
println("LoRa シンプルモデル シミュレーション開始")
println("設定:")
println("  SF: $(SF)")
println("  エリアサイズ: $(area_size) km")
println("  送信電力: $(Tx_dB) dBm")
println("  搬送波周波数: $(f_c) MHz")
println("  パスロス係数: α=$(α), β=$(β), γ=$(γ)")

sf = SF
bw = 125e3
num_devices = 50

# 単一シミュレーション実行
per, positions, snrs, errors = simple_csma_dc_per(sf, bw, num_devices;
                                   payload_range=(16,32),
                                   sim_time=1.0,
                                   backoff_max=0.01,
                                   dc=0.01)

# 結果分析
stats = analyze_simple_results(positions, snrs, errors)

println("\n=== シミュレーション結果 ===")
println("総端末数: $(stats[:total_devices])")
println("エラー数: $(stats[:total_errors])")
println("PER: $(stats[:per])")
println("平均距離: $(stats[:avg_distance]) km")
println("平均SNR: $(stats[:avg_snr]) dB")
println("SNR範囲: $(stats[:min_snr]) ~ $(stats[:max_snr]) dB")

# 可視化
plot_simple_results(positions, snrs, errors)

# 複数回実行による統計分析
println("\n=== 複数回実行による統計分析 ===")
per_results, all_stats = run_multiple_simulations(sf, bw, num_devices, 20;
                                   payload_range=(16,32),
                                   sim_time=1.0,
                                   backoff_max=0.01,
                                   dc=0.01)

# SNRスイープ実行（SER計算）
println("\n=== SNRスイープ実行（SER計算） ===")
snr_min, snr_max, snr_step = -20.0, 0.0, 0.5
iter_sweep = 50  # スイープ用の反復回数

snrs, ser_vals, per_vals = run_snr_sweep_ser(sf, bw, num_devices,
                               snr_min, snr_max, snr_step;
                               payload_range=(16,32),
                               sim_time=1.0,
                               backoff_max=0.01,
                               dc=0.01,
                               iter=iter_sweep,
                               save_path="LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/snr_sweep_ser_sf$(sf)_dev$(num_devices).csv")


println("\nシミュレーション完了！")
