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

# シャドウイング標準偏差[dB] - 追加
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
    # パスロス計算: PL = 10α*log10(d) + β + 10γ*log10(fc)
    return 10 * α * log10(distance_km) + β + 10 * γ * log10(f_c)
end

# シャドウイング値生成関数 - 追加
function shadowing_value(rng)
    # シャドウイング値生成（ログ正規分布）
    return randn(rng) * shadowing_std
end

# 受信電力計算関数を修正（シャドウイング追加）
function received_power(tx_power_dBm::Float64, distance_km::Float64, shadowing_dB::Float64)
    # 受信電力計算: Pr = Pt - PL - Shadowing
    pl = path_loss(distance_km)
    return tx_power_dBm - pl - shadowing_dB
end

# 既存の関数も保持（後方互換性のため）
function received_power(tx_power_dBm::Float64, distance_km::Float64)
    # シャドウイングなしの受信電力計算
    pl = path_loss(distance_km)
    return tx_power_dBm - pl
end

# SNR計算関数を修正（シャドウイング対応）
function calculate_snr(distance_km::Float64, shadowing_dB::Float64)
    # SNR計算（シャドウイング考慮）
    received_pwr = received_power(Tx_dB, distance_km, shadowing_dB)
    snr = received_pwr - (noise_power_spectrum_density + 10*log10(band_width) + noise_figure)
    return snr
end

# 既存の関数も保持
function calculate_snr(distance_km::Float64)
    # シャドウイングなしのSNR計算
    received_pwr = received_power(Tx_dB, distance_km)
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

# ===== 改良版CSMA + DC + PER（シャドウイング対応） =====
function simple_csma_dc_per(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    fixed_snr::Union{Float64, Nothing}=nothing,  # 固定SNR値（Noneの場合は距離ベース）
    cs_threshold_dBm::Union{Float64, Nothing}=nothing, # 受信電力しきい値（Noneで従来のbusy配列）
    use_shadowing::Bool=true,  # シャドウイング使用フラグ - 追加
    rng=Random.default_rng())

    M  = 2^sf
    Ts = M / bw     # シンボル長
    fs = bw         # サンプル周波数 ≈ 帯域幅

    # ポアソン点過程で端末配置生成
    node_positions = generate_poisson_positions(num_devices, rng)
    
    # シャドウイング値生成 - 追加
    shadowing_values = use_shadowing ? [shadowing_value(rng) for _ in 1:num_devices] : zeros(num_devices)
    
    # 各端末のSNR計算（固定SNRまたは距離ベース）
    if fixed_snr !== nothing
        device_snrs = fill(fixed_snr, num_devices)
    else
        if use_shadowing
            device_snrs = [calculate_snr(distance(pos[1], pos[2], 0.0, 0.0), shadowing_values[i]) for (i, pos) in enumerate(node_positions)]
        else
            device_snrs = [calculate_snr(distance(pos[1], pos[2], 0.0, 0.0)) for pos in node_positions]
        end
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

    # 既にスケジュール済みのシンボルウィンドウ (start_sample, end_sample, tx_device)
    scheduled_windows = Tuple{Int,Int,Int}[]

    # 受信電力ベースのキャリアセンス関数（デバイスd視点）
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
            while n1 <= total_samples && (
                (cs_threshold_dBm === nothing && any(channel_busy[n0:n1])) ||
                (cs_threshold_dBm !== nothing && sensed_power_at_device(d, n0, n1) >= cs_threshold_dBm)
            )
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
                # 他端末のキャリアセンス用にスケジュールを記録
                push!(scheduled_windows, (n0, n1, d))
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

    # シャドウイング値を返り値に追加
    return mean(Float64.(errors)), node_positions, device_snrs, errors, shadowing_values
end

# ===== 結果分析関数 =====
function analyze_simple_results(node_positions, device_snrs, errors, shadowing_values=nothing)
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
    
    # シャドウイング統計を追加
    if shadowing_values !== nothing
        stats[:avg_shadowing] = mean(shadowing_values)
        stats[:min_shadowing] = minimum(shadowing_values)
        stats[:max_shadowing] = maximum(shadowing_values)
        stats[:shadowing_std] = std(shadowing_values)
    end
    
    return stats
end

# ===== 可視化関数 =====

# シンプルな端末配置可視化関数
function plot_node_positions(node_positions; 
    sf::Int=7, num_devices::Int=100, save_prefix::String="node_positions")
    xs = [pos[1] for pos in node_positions]
    ys = [pos[2] for pos in node_positions]
    
    p = scatter(xs, ys,
        color = :blue,
        xlabel = "X [km]", ylabel = "Y [km]",
        title = "LoRa Node Distribution\nSF=$(sf), Devices=$(num_devices)",
        legend = false,
        markersize = 4,
        aspect_ratio = :equal,
        alpha = 0.7)
    
    # ゲートウェイを緑の星印でプロット
    scatter!(p, [0.0], [0.0], markershape = :star5, color = :green, label = "Gateway", markersize = 12)
    
    # エリア境界を円で表示
    theta = 0:0.01:2π
    circle_x = area_size .* cos.(theta)
    circle_y = area_size .* sin.(theta)
    plot!(p, circle_x, circle_y, color = :black, linestyle = :dash, linewidth = 2, label = "Area Boundary")
    
    # プロット情報を表示
    println("端末配置プロットを生成中...")
    println("総端末数: $(length(node_positions))")
    
    # プロットを表示
    display(p)
    
    # プロットをファイルに保存
    save_path = "LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/$(save_prefix)_sf$(sf)_dev$(num_devices).png"
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    savefig(p, save_path)
    println("プロットを保存しました: $save_path")
    
    return p
end

# エラー状態の可視化関数
function plot_node_positions_with_errors(node_positions, device_snrs, errors; 
    sf::Int=7, num_devices::Int=100, save_prefix::String="node_positions_errors")
    xs = [pos[1] for pos in node_positions]
    ys = [pos[2] for pos in node_positions]
    
    # エラー状態で色分け
    colors = [errors[i] ? :red : :blue for i in 1:length(errors)]
    
    p = scatter(xs, ys,
        color = colors,
        xlabel = "X [km]", ylabel = "Y [km]",
        title = "LoRa Node Distribution (Error Status)\nRed: Error, Blue: Success",
        legend = false,
        markersize = 6,
        aspect_ratio = :equal)
    
    # ゲートウェイを緑の星印でプロット
    scatter!(p, [0.0], [0.0], markershape = :star5, color = :green, label = "Gateway", markersize = 12)
    
    # エリア境界を円で表示
    theta = 0:0.01:2π
    circle_x = area_size .* cos.(theta)
    circle_y = area_size .* sin.(theta)
    plot!(p, circle_x, circle_y, color = :black, linestyle = :dash, linewidth = 2, label = "Area Boundary")
    
    # プロット情報を表示
    println("エラー状態プロットを生成中...")
    println("総端末数: $(length(node_positions))")
    println("エラー端末数: $(sum(errors))")
    println("成功率: $(1.0 - mean(Float64.(errors)))")
    
    # プロットを表示
    display(p)
    
    # プロットをファイルに保存
    save_path = "LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/$(save_prefix)_sf$(sf)_dev$(num_devices).png"
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    savefig(p, save_path)
    println("プロットを保存しました: $save_path")
    
    return p
end

# SNR分布の可視化
function plot_snr_distribution(node_positions, device_snrs; 
    sf::Int=7, num_devices::Int=100, save_prefix::String="snr_distribution")
    xs = [pos[1] for pos in node_positions]
    ys = [pos[2] for pos in node_positions]
    
    p = scatter(xs, ys,
        zcolor = device_snrs,  # zcolorを使用
        xlabel = "X [km]", ylabel = "Y [km]",
        title = "LoRa Node SNR Distribution\nSF=$(sf), Devices=$(num_devices)",
        legend = false,
        markersize = 6,
        aspect_ratio = :equal,
        colorbar_title = "SNR [dB]",
        colorbar_titlefontsize = 10)
    
    # ゲートウェイを緑の星印でプロット
    scatter!(p, [0.0], [0.0], markershape = :star5, color = :green, label = "Gateway", markersize = 12)
    
    # エリア境界を円で表示
    theta = 0:0.01:2π
    circle_x = area_size .* cos.(theta)
    circle_y = area_size .* sin.(theta)
    plot!(p, circle_x, circle_y, color = :black, linestyle = :dash, linewidth = 2, label = "Area Boundary")
    
    # プロット情報を表示
    println("SNR分布プロットを生成中...")
    println("総端末数: $(length(node_positions))")
    println("平均SNR: $(mean(device_snrs)) dB")
    println("SNR範囲: $(minimum(device_snrs)) ~ $(maximum(device_snrs)) dB")
    
    # プロットを表示
    display(p)
    
    # プロットをファイルに保存
    save_path = "LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/$(save_prefix)_sf$(sf)_dev$(num_devices).png"
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    savefig(p, save_path)
    println("プロットを保存しました: $save_path")
    
    return p
end

# 統合可視化関数
function visualize_simulation_results(node_positions, device_snrs, errors; 
    sf::Int=7, num_devices::Int=100, save_prefix::String="simulation_results")
    
    println("\n=== シミュレーション結果可視化 ===")
    
    # 1. 基本的な端末配置
    p1 = plot_node_positions(node_positions; sf=sf, num_devices=num_devices, save_prefix=save_prefix)
    
    # 2. SNR分布
    p3 = plot_snr_distribution(node_positions, device_snrs; sf=sf, num_devices=num_devices, save_prefix=save_prefix)
    
    # 3. エラー状態
    p4 = plot_node_positions_with_errors(node_positions, device_snrs, errors; sf=sf, num_devices=num_devices, save_prefix=save_prefix)
    
    # 4. 統計情報の表示
    distances = [distance(pos[1], pos[2], 0.0, 0.0) for pos in node_positions]
    println("\n=== 統計情報 ===")
    println("総端末数: $(length(node_positions))")
    println("エラー端末数: $(sum(errors))")
    println("成功率: $(1.0 - mean(Float64.(errors)))")
    println("平均距離: $(mean(distances)) km")
    println("平均SNR: $(mean(device_snrs)) dB")
    println("SNR範囲: $(minimum(device_snrs)) ~ $(maximum(device_snrs)) dB")
    
    return p1, p3, p4
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

# ===== SNRスイープ機能（PER計算のみ） =====
function run_snr_sweep_per(sf::Int, bw::Float64, num_devices::Int,
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    iter::Int=100,
    seed::Int=1234,
    save_path::String="LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/snr_sweep_per.csv",
    cs_threshold_dBm::Union{Float64,Nothing}=nothing,   # 追加
    use_shadowing::Bool=true)                           # 追加

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    per_vals = zeros(Float64, length(snrs))
    
    println("SNRスイープ実行中...")
    println("SNR範囲: $(snr_min) ~ $(snr_max) dB, ステップ: $(snr_step) dB")
    println("反復回数: $iter")
    
    for (i, snr) in enumerate(snrs)
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
                                               cs_threshold_dBm=cs_threshold_dBm,   # 追加
                                               use_shadowing=use_shadowing,         # 追加
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
    
    @info "SNRスイープ結果をCSVに保存しました: $save_path"
    
    return snrs, per_vals
end

# 理論的なSER計算関数を削除（不要になったため）

# ===== 実行部分 =====
println("LoRa シンプルモデル シミュレーション開始")
println("設定:")
println("  SF: $(SF)")
println("  エリアサイズ: $(area_size) km")
println("  送信電力: $(Tx_dB) dBm")
println("  搬送波周波数: $(f_c) MHz")
println("  パスロス係数: α=$(α), β=$(β), γ=$(γ)")
println("  シャドウイング標準偏差: $(shadowing_std) dB")

sf = SF
bw = 125e3
num_devices = 2

# 実行フラグ（SNRスイープ以外はデフォルトで無効化）
RUN_SINGLE = false  # 可視化のため有効化
RUN_MULTIPLE = false
RUN_SNR_SWEEP = true  # PER計算のみに変更
RUN_VISUALIZATION_ONLY = false  # 可視化のみ実行

if RUN_VISUALIZATION_ONLY
    # 可視化のみ実行（テスト用）
    println("\n=== 可視化テスト実行 ===")
    test_positions = generate_poisson_positions(num_devices, Random.default_rng())
    plot_node_positions(test_positions; sf=sf, num_devices=num_devices, save_prefix="nodes_only_test")
end

if RUN_SINGLE
    per, positions, snrs, errors = simple_csma_dc_per(sf, bw, num_devices;
                                       payload_range=(16,32),
                                       sim_time=1.0,
                                       backoff_max=0.01,
                                       dc=0.01)

    stats = analyze_simple_results(positions, snrs, errors)
    println("\n=== シミュレーション結果 ===")
    println("総端末数: $(stats[:total_devices])")
    println("エラー数: $(stats[:total_errors])")
    println("PER: $(stats[:per])")
    println("平均距離: $(stats[:avg_distance]) km")
    println("平均SNR: $(stats[:avg_snr]) dB")
    println("SNR範囲: $(stats[:min_snr]) ~ $(stats[:max_snr]) dB")
    
    # 統合可視化実行 → 端末配置のみ
    plot_node_positions(positions; sf=sf, num_devices=num_devices, save_prefix="nodes_only")
end

if RUN_MULTIPLE
    println("\n=== 複数回実行による統計分析 ===")
    per_results, all_stats = run_multiple_simulations(sf, bw, num_devices, 20;
                                       payload_range=(16,32),
                                       sim_time=1.0,
                                       backoff_max=0.01,
                                       dc=0.01)
end

if RUN_SNR_SWEEP
    println("\n=== SNRスイープ実行（PER計算のみ） ===")
    snr_min, snr_max, snr_step = -20.0, 0.0, 0.5
    iter_sweep = 1000 # スイープ用の反復回数
    snrs, per_vals = run_snr_sweep_per(sf, bw, num_devices,
                                   snr_min, snr_max, snr_step;
                                   payload_range=(16,32),
                                   sim_time=1.0,
                                   backoff_max=0.01,
                                   dc=0.01,
                                   iter=iter_sweep,
                                   cs_threshold_dBm=-110.0,
                                   use_shadowing=true,
                                   save_path="LoRa_SER/LoRa_Pathloss_CSMA_DC/results_simple_model/snr_sweep_per_sf$(sf)_dev$(num_devices)_iter$(iter_sweep).csv")
end

println("\nシミュレーション完了！")
