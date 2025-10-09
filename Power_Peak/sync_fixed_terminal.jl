using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots, Dates

# ===== 信号パラメータ =====
struct SignalParameters
    duration_s::Float64       # 信号持続時間（s）
    center_frequency_hz::Float64 # 中心周波数（Hz）
    bandwidth_hz::Float64     # 帯域幅（Hz）
    sampling_rate_hz::Float64 # サンプリングレート（Hz）
    tx_power_dbm::Float64     # 送信電力（dBm）
end

# ===== ノイズパラメータ =====
struct NoiseParameters
    noise_power_dbm::Float64  # ノイズ電力（dBm）
    snr_db::Float64          # 信号対雑音比（dB）
end

# ===== パスロスパラメータ =====
struct PathLossParameters
    distance_m::Float64       # 送信機と受信機の距離（m）
    frequency_hz::Float64     # 周波数（Hz）
    path_loss_exponent::Float64  # パスロス指数（通常2.0-4.0）
    reference_distance_m::Float64  # 参照距離（m）
    reference_path_loss_db::Float64  # 参照距離でのパスロス（dB）
end

# ===== ポアソン点過程端末配置パラメータ =====
struct TerminalDeploymentParameters
    deployment_mode::String   # 配置モード: "poisson" または "fixed"
    lambda::Float64           # ポアソン点過程の密度（点/m²）
    num_terminals::Int        # 固定端末数（fixedモード時）
    area_size_m::Float64      # エリアサイズ（m）
    min_distance_m::Float64   # 最小距離（m）
    max_distance_m::Float64  # 最大距離（m）
    frequency_hz::Float64    # 周波数（Hz）
    path_loss_exponent::Float64  # パスロス指数
    reference_distance_m::Float64  # 参照距離（m）
    reference_path_loss_db::Float64  # 参照距離でのパスロス（dB）
end

# ===== シャドウイングパラメータ =====
struct ShadowingParameters
    enabled::Bool            # シャドウイング有効/無効
    std_db::Float64         # シャドウイング標準偏差（dB）
    correlation_distance_m::Float64  # 相関距離（m）
    correlation_coefficient::Float64  # 相関係数
end

# ===== 端末情報 =====
struct TerminalInfo
    x_m::Float64             # X座標（m）
    y_m::Float64             # Y座標（m）
    distance_m::Float64      # 基地局からの距離（m）
    path_loss_db::Float64    # パスロス（dB）
    shadowing_db::Float64    # シャドウイング（dB）
    total_loss_db::Float64   # 総損失（パスロス+シャドウイング）（dB）
    rx_power_dbm::Float64    # 受信電力（dBm）
end

# ===== シミュレーションパラメータ（一括設定） =====
struct SimulationParameters
    # 信号パラメータ
    signal_duration_us::Float64      # 信号持続時間（μs）
    center_frequency_ghz::Float64     # 中心周波数（GHz）
    signal_bandwidth_mhz::Float64    # 同期信号帯域幅（MHz）
    terminal_bandwidth_mhz::Float64  # 端末受信帯域幅（MHz）
    sampling_rate_mhz::Float64       # サンプリングレート（MHz）
    tx_power_dbm::Float64            # 送信電力（dBm）
    
    # 受信環境パラメータ
    snr_db::Float64                  # 信号対雑音比（dB）
    shadowing_enabled::Bool          # シャドウイング有効/無効
    shadowing_std_db::Float64        # シャドウイング標準偏差（dB）
    
    # 端末配置パラメータ
    deployment_mode::String          # 配置モード: "poisson" または "fixed"
    num_terminals::Int               # 端末数（fixedモード時）
    area_size_m::Float64            # エリアサイズ（m）
    min_distance_m::Float64         # 最小距離（m）
    max_distance_m::Float64          # 最大距離（m）
    
    # パスロスパラメータ
    path_loss_exponent::Float64      # パスロス指数
    reference_distance_m::Float64   # 参照距離（m）
    reference_path_loss_db::Float64  # 参照距離でのパスロス（dB）
    
    # シミュレーション設定
    signal_interval_ms::Float64     # 信号間隔（ms）
    total_duration_ms::Float64       # 総持続時間（ms）
    power_threshold_dbm::Float64     # ピーク検出閾値（dBm）
end

# ===== QPSKシンボル生成関数 =====
"""
QPSKシンボル（4つの候補の中からランダムに1つ）を生成する。
平均電力が1になるように正規化済み。
"""
function generate_qpsk_symbol()
    # 実部と虚部がそれぞれ-1か+1のどちらかをランダムに選ぶ
    real_part = rand([-1, 1])
    imag_part = rand([-1, 1])
    
    # 複素数を生成し、平均電力が1になるようにsqrt(2)で割る
    return (real_part + im * imag_part) / sqrt(2)
end

# ===== 送信信号生成関数 =====
"""
周波数領域で、各サブキャリアにランダムなQPSKシンボルを配置することで
広帯域信号を生成する。
"""
function generate_ofdm_signal_with_qpsk(bandwidth_hz::Float64, duration_s::Float64, fs_hz::Float64)
    # FFTのサイズ（＝時間領域のサンプル数）を決定
    N = Int(ceil(duration_s * fs_hz))
    
    # 周波数領域の配列をゼロで初期化
    freq_domain_signal = zeros(ComplexF64, N)
    
    # 信号を配置する帯域幅を周波数ビンの数に変換
    bw_bins = Int(floor(bandwidth_hz / (fs_hz / N)))
    
    # 中心付近の指定帯域幅のビンにランダムなシンボルを設定
    start_bin = div(N, 2) - div(bw_bins, 2)
    end_bin = div(N, 2) + div(bw_bins, 2)
    
    for i in start_bin:end_bin
        # 各サブキャリアに、ランダムなQPSKシンボルを1つ配置する
        freq_domain_signal[i] = generate_qpsk_symbol()
    end
    
    # IFFTを適用して時間領域の信号に変換
    time_domain_signal_complex = ifft(ifftshift(freq_domain_signal))
    
    return time_domain_signal_complex
end

# ===== 送信側：周期的同期信号生成 =====
"""
周期的に同期信号を送信する送信機のシミュレーション
"""
function generate_periodic_sync_signals(params::SignalParameters, interval_ms::Float64, total_duration_ms::Float64)
    # パラメータ設定
    signal_duration_ms = params.duration_s * 1000
    sampling_rate_hz = params.sampling_rate_hz
    total_samples = Int(ceil(total_duration_ms * sampling_rate_hz / 1000))
    
    # 信号送信タイミングを計算
    signal_times = Float64[]
    current_time = 0.0
    while current_time + signal_duration_ms <= total_duration_ms
        push!(signal_times, current_time)
        current_time += interval_ms
    end
    
    println("送信側パラメータ:")
    println("• 信号間隔: $(interval_ms) ms")
    println("• 信号持続時間: $(signal_duration_ms) ms")
    println("• 総持続時間: $(total_duration_ms) ms")
    println("• 送信信号数: $(length(signal_times))")
    println("• 送信電力: $(params.tx_power_dbm) dBm")
    println()
    
    # 送信信号配列を初期化
    tx_signal = zeros(ComplexF64, total_samples)
    time_axis = collect((0:total_samples-1) / sampling_rate_hz * 1000)  # ms単位
    
    # 各送信タイミングで信号を生成・配置
    signal_count = 0
    for (i, signal_time) in enumerate(signal_times)
        # 個別のQPSK信号を生成
        individual_signal = generate_ofdm_signal_with_qpsk(
            params.bandwidth_hz,
            params.duration_s,
            sampling_rate_hz
        )
        
        # 送信電力で正規化
        tx_power_linear = 10^(params.tx_power_dbm / 10) * 1e-3  # mWからWに変換
        individual_signal = individual_signal * sqrt(tx_power_linear)
        
        # 信号を送信配列に配置
        signal_samples = length(individual_signal)
        start_sample = Int(ceil(signal_time * sampling_rate_hz / 1000)) + 1
        end_sample = min(start_sample + signal_samples - 1, total_samples)
        
        if start_sample <= total_samples
            actual_samples = end_sample - start_sample + 1
            tx_signal[start_sample:end_sample] = individual_signal[1:actual_samples]
            signal_count += 1
            println("送信信号 $(signal_count): 時刻 $(round(signal_time, digits=3)) ms ～ $(round(signal_time + signal_duration_ms, digits=3)) ms")
        end
    end
    
    return time_axis, tx_signal, signal_count
end

# ===== シャドウイング計算関数 =====
"""
シャドウイングを計算する（固定値またはランダム）
"""
function calculate_shadowing(shadowing_params::ShadowingParameters, distance_m::Float64, fixed_value::Bool=true)
    if !shadowing_params.enabled
        return 0.0
    end
    
    if fixed_value
        # 固定値のシャドウイング（全端末で同じ値）
        return shadowing_params.std_db
    else
        # 対数正規分布によるシャドウイング（ランダム）
        shadowing_db = randn() * shadowing_params.std_db
        return shadowing_db
    end
end

# ===== 端末配置関数 =====
"""
ポアソン点過程または固定端末数で端末を配置する
"""
function deploy_terminals(deployment_params::TerminalDeploymentParameters, shadowing_params::ShadowingParameters, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    area_size = deployment_params.area_size_m
    
    if deployment_params.deployment_mode == "poisson"
        # ポアソン点過程モード
        lambda = deployment_params.lambda
        num_points = rand(Poisson(lambda * area_size^2))
        println("ポアソン点過程モード: 密度$(lambda)点/m², 生成点数$(num_points)")
        
    elseif deployment_params.deployment_mode == "fixed"
        # 固定端末数モード
        num_points = deployment_params.num_terminals
        println("固定端末数モード: 端末数$(num_points)（固定位置配置）")
        
    else
        error("不正な配置モード: $(deployment_params.deployment_mode)")
    end
    
    # 固定端末位置の設定
    fixed_positions = [
        (250.0, 250.0),    # 東方向 50m
        (-50.0, 0.0),   # 西方向 50m
        (0.0, 50.0),    # 北方向 50m
        (0.0, -50.0),   # 南方向 50m
        (35.4, 35.4),   # 北東方向 50m
        (-35.4, 35.4),  # 北西方向 50m
        (35.4, -35.4),  # 南東方向 50m
        (-35.4, -35.4)  # 南西方向 50m
    ]
    
    # 固定位置に端末を配置
    for i in 1:min(num_points, length(fixed_positions))
        # 固定位置を取得
        x, y = fixed_positions[i]
        
        # 基地局からの距離を計算（確認用）
        distance = sqrt(x^2 + y^2)
        
        # 固定位置の距離制限をチェック
        if distance >= deployment_params.min_distance_m && distance <= deployment_params.max_distance_m
            # パスロスを計算
            distance_ratio = distance / deployment_params.reference_distance_m
            path_loss_db = deployment_params.reference_path_loss_db + 
                          10 * deployment_params.path_loss_exponent * log10(distance_ratio)
            
            # シャドウイングを計算（端末ごとに異なる値）
            shadowing_db = calculate_shadowing(shadowing_params, distance, false)
            
            # 総損失を計算
            total_loss_db = path_loss_db + shadowing_db
            
            # 受信電力を計算（パラメータから取得した送信電力）
            rx_power_dbm = tx_power_dbm - total_loss_db
            
            # 端末情報を作成
            terminal = TerminalInfo(x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm)
            push!(terminals, terminal)
        end
    end
    
    return terminals
end

# ===== パスロス計算関数 =====
"""
距離に基づくパスロスを計算する
"""
function calculate_path_loss(path_loss_params::PathLossParameters)
    # フリースペースパスロス + 距離による減衰
    # PL(d) = PL(d0) + 10*n*log10(d/d0)
    # ここで n はパスロス指数、d0 は参照距離
    
    distance_ratio = path_loss_params.distance_m / path_loss_params.reference_distance_m
    path_loss_db = path_loss_params.reference_path_loss_db + 
                  10 * path_loss_params.path_loss_exponent * log10(distance_ratio)
    
    return path_loss_db
end

# ===== 受信側：パスロス + AWGN付加 =====
"""
受信信号にパスロスとAWGN（Additive White Gaussian Noise）を付加する
"""
function add_path_loss_and_noise(tx_signal::Vector{ComplexF64}, path_loss_params::PathLossParameters, 
                                noise_params::NoiseParameters, shadowing_params::ShadowingParameters, sampling_rate_hz::Float64)
    N = length(tx_signal)
    
    # パスロス計算
    path_loss_db = calculate_path_loss(path_loss_params)
    path_loss_linear = 10^(-path_loss_db / 10)  # dBから線形値に変換
    
    # シャドウイング計算
    shadowing_db = calculate_shadowing(shadowing_params, path_loss_params.distance_m, false)
    shadowing_linear = 10^(-shadowing_db / 10)  # dBから線形値に変換
    
    # パスロスとシャドウイングを適用した受信信号
    total_loss_linear = path_loss_linear * shadowing_linear
    rx_signal_with_path_loss = tx_signal * sqrt(total_loss_linear)
    
    # ノイズなし信号の電力計算
    signal_power_linear = mean(abs2.(rx_signal_with_path_loss))
    signal_power_dbm = 10 * log10(signal_power_linear * 1000)  # WからdBmに変換
    
    # 固定ノイズ床（受信機の熱雑音）- 物理的に正しいモデル
    noise_power_dbm = noise_params.noise_power_dbm  # 固定値（-109.0 dBm）
    
    # ノイズ電力の線形値に変換
    noise_power_linear = 10^(noise_power_dbm / 10) * 1e-3  # dBmからWに変換
    
    # AWGN（複素ガウシアンノイズ）を生成
    # 実部と虚部が独立した正規分布（平均0、分散σ²/2）
    noise_std = sqrt(noise_power_linear / 2)  # 実部・虚部の標準偏差
    noise_real = randn(N) * noise_std
    noise_imag = randn(N) * noise_std
    noise = noise_real + im * noise_imag
    
    # 受信信号にノイズを付加
    rx_signal_with_noise = rx_signal_with_path_loss + noise
    
    # ノイズ付加後の受信電力計算
    rx_power_linear = mean(abs2.(rx_signal_with_noise))
    rx_power_dbm = 10 * log10(rx_power_linear * 1000)  # WからdBmに変換
    
    # 実際のSNRを計算
    actual_noise_power = mean(abs2.(noise))
    actual_snr_db = 10 * log10(signal_power_linear / actual_noise_power)
    
    println("受信側パスロス・シャドウイング・AWGNパラメータ:")
    println("• 距離: $(round(path_loss_params.distance_m, digits=1)) m")
    println("• パスロス: $(round(path_loss_db, digits=2)) dB")
    println("• シャドウイング: $(round(shadowing_db, digits=2)) dB")
    println("• 総損失: $(round(path_loss_db + shadowing_db, digits=2)) dB")
    println("• 受信信号電力（ノイズなし）: $(round(signal_power_dbm, digits=2)) dBm")
    println("• 受信信号電力（ノイズあり）: $(round(rx_power_dbm, digits=2)) dBm")
    println("• 固定ノイズ電力: $(round(noise_power_dbm, digits=2)) dBm")
    println("• 実際のSNR: $(round(actual_snr_db, digits=2)) dB")
    println()
    
    return rx_signal_with_noise, actual_snr_db, path_loss_db
end

# ===== 受信電力ピーク検出 =====
"""
受信電力から同期信号のピークを検出する
"""
function detect_sync_peaks(rx_power::Vector{Float64}, time_axis::Vector{Float64}, 
                          signal_count::Int, interval_ms::Float64, power_threshold_db::Float64)
    # 電力閾値を線形値に変換
    power_threshold_linear = 10^(power_threshold_db / 10) * 1e-3  # dBmからWに変換
    
    # ピーク検出用の閾値（平均電力の倍数）
    avg_power = mean(rx_power)
    detection_threshold = max(power_threshold_linear, 3 * avg_power)
    
    # ピーク検出
    detected_peaks = []
    peak_times = Float64[]
    peak_powers = Float64[]
    
    # 各信号送信タイミング周辺でピークを検索
    for i in 1:signal_count
        expected_start = (i-1) * interval_ms
        expected_end = expected_start + 1.0  # 1msの検索窓
        
        # 時間窓のインデックスを計算
        start_idx = max(1, Int(ceil(expected_start * length(rx_power) / time_axis[end])))
        end_idx = min(length(rx_power), Int(ceil(expected_end * length(rx_power) / time_axis[end])))
        
        if start_idx <= end_idx
            window_power = rx_power[start_idx:end_idx]
            window_time = time_axis[start_idx:end_idx]
            
            # この窓での最大電力とそのインデックス
            max_power = maximum(window_power)
            max_idx = argmax(window_power)
            
            if max_power > detection_threshold
                peak_time = window_time[max_idx]
                push!(detected_peaks, (peak_time, max_power))
                push!(peak_times, peak_time)
                push!(peak_powers, max_power)
            end
        end
    end
    
    detection_rate = length(detected_peaks) / signal_count * 100
    
    println("受信電力ピーク検出結果:")
    println("• 検出閾値: $(round(detection_threshold*1000, digits=3)) mW ($(round(10*log10(detection_threshold*1000), digits=1)) dBm)")
    println("• 検出ピーク数: $(length(detected_peaks))/$signal_count")
    println("• 検出率: $(round(detection_rate, digits=1))%")
    
    if !isempty(peak_powers)
        println("• ピーク電力範囲: $(round(minimum(peak_powers)*1000, digits=3)) - $(round(maximum(peak_powers)*1000, digits=3)) mW")
        println("• ピーク電力範囲: $(round(10*log10(minimum(peak_powers)*1000), digits=1)) - $(round(10*log10(maximum(peak_powers)*1000), digits=1)) dBm")
    end
    println()
    
    return detected_peaks, peak_times, peak_powers, detection_rate
end

# ===== 受信電力可視化関数 =====
"""
受信電力とピーク検出結果の可視化
"""
function plot_received_power(time_axis::Vector{Float64}, rx_power::Vector{Float64}, 
                           peak_times::Vector{Float64}, peak_powers::Vector{Float64},
                           title::String)
    # データを間引いてプロット（最大10000点）
    max_points = 10000
    step = max(1, div(length(time_axis), max_points))
    time_subset = time_axis[1:step:end]
    power_subset = rx_power[1:step:end]
    
    # 受信電力（線形スケール）
    p1 = plot(time_subset, power_subset * 1000, label="Received Power", linewidth=1, color=:blue,
              xlabel="Time (ms)", ylabel="Power (mW)", title="$(title) - Received Power (Linear)")
    
    # 検出されたピークをマーカーで表示
    if !isempty(peak_times)
        scatter!(p1, peak_times, peak_powers * 1000, label="Detected Peaks", 
                markersize=8, color=:red, marker=:circle)
    end
    
    # 受信電力（対数スケール）
    p2 = plot(time_subset, 10*log10.(power_subset * 1000), label="Received Power", linewidth=1, color=:blue,
              xlabel="Time (ms)", ylabel="Power (dBm)", title="$(title) - Received Power (Log)")
    
    # 検出されたピークをマーカーで表示
    if !isempty(peak_times)
        scatter!(p2, peak_times, 10*log10.(peak_powers * 1000), label="Detected Peaks", 
                markersize=8, color=:red, marker=:circle)
    end
    
    # 電力ヒストグラム
    p3 = histogram(power_subset * 1000, bins=50, label="Power Distribution", color=:green, alpha=0.7,
                   xlabel="Power (mW)", ylabel="Count", title="$(title) - Power Histogram")
    
    # ピーク電力の統計
    if !isempty(peak_powers)
        peak_stats = [
            ("Min Peak", minimum(peak_powers) * 1000),
            ("Max Peak", maximum(peak_powers) * 1000),
            ("Avg Peak", mean(peak_powers) * 1000),
            ("Peak Count", length(peak_powers))
        ]
        
        stats_text = join(["$name: $(round(val, digits=3))" for (name, val) in peak_stats], "\n")
        annotate!(p3, 0.7, 0.8, text(stats_text, 10, :left))
    end
    
    return plot(p1, p2, p3, layout=(2,2), size=(1200, 800))
end

# ===== 受信電力CSV出力関数 =====
function save_received_power_to_csv(time_axis::Vector{Float64}, rx_power::Vector{Float64}, 
                                  peak_times::Vector{Float64}, peak_powers::Vector{Float64},
                                  params::SignalParameters, noise_params::NoiseParameters, 
                                  terminals::Vector{TerminalInfo}, representative_terminal::TerminalInfo,
                                  signal_count::Int, detection_rate::Float64, output_dir::String, shadowing_params::ShadowingParameters)
    # 実行時刻を取得
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    # 受信電力データ（高サンプリングレート）
    power_data = DataFrame(
        time_ms = time_axis,
        power_mw = rx_power * 1000,  # WからmWに変換
        power_dbm = 10*log10.(rx_power * 1000)  # dBmに変換
    )
    
    # 検出されたピークデータ
    if !isempty(peak_times)
        peak_data = DataFrame(
            peak_time_ms = peak_times,
            peak_power_mw = peak_powers * 1000,
            peak_power_dbm = 10*log10.(peak_powers * 1000),
            peak_index = [findfirst(x -> abs(x - t) < 1e-6, time_axis) for t in peak_times]
        )
    else
        peak_data = DataFrame(
            peak_time_ms = Float64[],
            peak_power_mw = Float64[],
            peak_power_dbm = Float64[],
            peak_index = Int[]
        )
    end
    
    # 端末配置データ
    terminal_data = DataFrame(
        terminal_id = collect(1:length(terminals)),
        x_m = [t.x_m for t in terminals],
        y_m = [t.y_m for t in terminals],
        distance_m = [t.distance_m for t in terminals],
        path_loss_db = [t.path_loss_db for t in terminals],
        shadowing_db = [t.shadowing_db for t in terminals],
        total_loss_db = [t.total_loss_db for t in terminals],
        rx_power_dbm = [t.rx_power_dbm for t in terminals],
        is_representative = [t == representative_terminal for t in terminals]
    )
    
    # パラメータ情報
    params_data = DataFrame(
        parameter = ["Duration (ms)", "Bandwidth (MHz)", "Center Frequency (GHz)", "Sampling Rate (MHz)", 
                    "TX Power (dBm)", "Noise Power (dBm)", "SNR (dB)", "Num Terminals", "Min Distance (m)", 
                    "Max Distance (m)", "Path Loss Exponent", "Shadowing Enabled", "Shadowing Std (dB)", 
                    "Representative Distance (m)", "Representative Path Loss (dB)", "Representative Shadowing (dB)", 
                    "Representative Total Loss (dB)", "Representative RX Power (dBm)", "Signal Interval (ms)", 
                    "Total Signals", "Detected Peaks", "Detection Rate (%)", "Execution Time"],
        value = [params.duration_s*1000, params.bandwidth_hz/1e6, params.center_frequency_hz/1e9, 
                params.sampling_rate_hz/1e6, params.tx_power_dbm, noise_params.noise_power_dbm, 
                noise_params.snr_db, length(terminals), minimum([t.distance_m for t in terminals]),
                maximum([t.distance_m for t in terminals]), terminals[1].path_loss_db / (10 * log10(terminals[1].distance_m / 1.0)),
                shadowing_params.enabled, shadowing_params.std_db, representative_terminal.distance_m, 
                round(representative_terminal.path_loss_db, digits=2), round(representative_terminal.shadowing_db, digits=2),
                round(representative_terminal.total_loss_db, digits=2), round(representative_terminal.rx_power_dbm, digits=2), 
                20.0, signal_count, length(peak_times), detection_rate, timestamp]
    )
    
    # CSVファイルに保存（received_powerとsync_simulationのみ）
    CSV.write("$(output_dir)/received_power_data_$(timestamp).csv", power_data)
    CSV.write("$(output_dir)/sync_simulation_parameters_$(timestamp).csv", params_data)
    
    return power_data, peak_data, params_data
end

# ===== パラメータ設定関数 =====
"""
シミュレーションパラメータを一括設定する
"""
function create_simulation_parameters()
    # ===== シミュレーションパラメータ（ここで一括変更可能） =====
    return SimulationParameters(
        # 信号パラメータ
        142.8,              # 信号持続時間（μs）
        4.7,                # 中心周波数（GHz）
        1.0,                # 同期信号帯域幅（MHz）- 広帯域信号
        0.01,                # 端末受信帯域幅（MHz）- 狭帯域受信
        2.0,                # サンプリングレート（MHz）- 同期信号帯域幅の2倍
        20.0,               # 送信電力（dBm）
        
        # 受信環境パラメータ
        0.0,                # SNR（dB）- 非常に低いSNR（悪い受信環境） -15dB
        true,               # シャドウイング有効/無効
        8.0,                # シャドウイング標準偏差（dB）- 都市環境
        
        # 端末配置パラメータ
        "fixed",            # 配置モード: "poisson" または "fixed"
        1,                  # 端末数（fixedモード時）
        100.0,              # エリアサイズ（m）
        10.0,               # 最小距離（m）
        500.0,              # 最大距離（m）
        
        # パスロスパラメータ
        3.0,                # パスロス指数（フリースペース）
        1.0,                # 参照距離（m）
        32.4,               # 参照距離でのパスロス（dB）- 1m@4.7GHz
        
        # シミュレーション設定
        20.0,               # 信号間隔（ms）
        110.0,              # 総持続時間（ms）
        -130.0              # ピーク検出閾値（dBm）- 狭帯域用に調整
    )
end

# ===== パラメータ設定例 =====
"""
環境別パラメータ設定例（コメントアウトを外して使用）
"""
function create_parameter_examples()
    # 郊外環境（軽微なシャドウイング）
    # return SimulationParameters(142.8, 4.7, 1.0, 2.0, 20.0, 15.0, true, 4.0, "fixed", 1, 100.0, 10.0, 200.0, 2.0, 1.0, 32.4, 20.0, 100.0, -60.0)
    
    # 都市環境（標準的なシャドウイング）
    # return SimulationParameters(142.8, 4.7, 1.0, 2.0, 20.0, 10.0, true, 8.0, "fixed", 1, 100.0, 10.0, 200.0, 2.0, 1.0, 32.4, 20.0, 100.0, -60.0)
    
    # 密集都市環境（強いシャドウイング）
    # return SimulationParameters(142.8, 4.7, 1.0, 2.0, 20.0, 5.0, true, 12.0, "fixed", 1, 100.0, 10.0, 200.0, 3.0, 1.0, 32.4, 20.0, 100.0, -60.0)
    
    # 理想環境（シャドウイング無効）
    # return SimulationParameters(142.8, 4.7, 1.0, 2.0, 20.0, 25.0, false, 0.0, "fixed", 1, 100.0, 10.0, 200.0, 2.0, 1.0, 32.4, 20.0, 100.0, -60.0)
    
    # 複数端末（10台）
    # return SimulationParameters(142.8, 4.7, 1.0, 2.0, 20.0, 10.0, true, 8.0, "fixed", 10, 150.0, 20.0, 300.0, 2.0, 1.0, 32.4, 20.0, 100.0, -60.0)
    
    # ポアソン点過程（低密度）
    # return SimulationParameters(142.8, 4.7, 1.0, 2.0, 20.0, 10.0, true, 8.0, "poisson", 0, 100.0, 10.0, 200.0, 2.0, 1.0, 32.4, 20.0, 100.0, -60.0)
    
    # デフォルト設定
    return create_simulation_parameters()
end

# ===== メイン実行関数 =====
function main()
    # 乱数シードを固定（再現性のため）- 最初に設定
    Random.seed!(1234)
    
    println("="^60)
    title_str = "同期信号受信シミュレーション（受信電力出力）"
    println(title_str)
    println("="^60)

    # ===== パラメータ設定（一括） =====
    sim_params = create_parameter_examples()
    
    # パラメータを表示
    println("シミュレーションパラメータ:")
    println("• 信号持続時間: $(sim_params.signal_duration_us) μs")
    println("• 中心周波数: $(sim_params.center_frequency_ghz) GHz")
    println("• 同期信号帯域幅: $(sim_params.signal_bandwidth_mhz) MHz")
    println("• 端末受信帯域幅: $(sim_params.terminal_bandwidth_mhz) MHz")
    println("• サンプリングレート: $(sim_params.sampling_rate_mhz) MHz")
    println("• 送信電力: $(sim_params.tx_power_dbm) dBm")
    println("• SNR: $(sim_params.snr_db) dB")
    println("• シャドウイング: $(sim_params.shadowing_enabled ? "有効" : "無効")")
    if sim_params.shadowing_enabled
        println("• シャドウイング標準偏差: $(sim_params.shadowing_std_db) dB")
    end
    println("• 配置モード: $(sim_params.deployment_mode)")
    println("• 端末数: $(sim_params.num_terminals)")
    println("• エリアサイズ: $(sim_params.area_size_m)m×$(sim_params.area_size_m)m")
    println("• パスロス指数: $(sim_params.path_loss_exponent)")
    println("• 信号間隔: $(sim_params.signal_interval_ms) ms")
    println("• 総持続時間: $(sim_params.total_duration_ms) ms")
    println("• ピーク検出閾値: $(sim_params.power_threshold_dbm) dBm")
    println()
    
    # パラメータを既存の構造体に変換
    params = SignalParameters(
        sim_params.signal_duration_us * 1e-6,  # μs → s
        sim_params.center_frequency_ghz * 1e9,  # GHz → Hz
        sim_params.signal_bandwidth_mhz * 1e6,  # MHz → Hz (同期信号帯域幅)
        sim_params.sampling_rate_mhz * 1e6,    # MHz → Hz
        sim_params.tx_power_dbm
    )
    
    # ノイズパラメータ設定（固定ノイズ床: 熱雑音 + ノイズフィギュア）
    # 熱雑音密度: -174 dBm/Hz（端末受信帯域幅を使用）
    noise_figure_db = 5.0
    terminal_bandwidth_hz = sim_params.terminal_bandwidth_mhz * 1e6  # MHz → Hz
    fixed_noise_power_dbm = -174 + 10 * log10(terminal_bandwidth_hz) + noise_figure_db
    noise_params = NoiseParameters(fixed_noise_power_dbm, sim_params.snr_db)
    
    # 端末配置パラメータ設定
    deployment_params = TerminalDeploymentParameters(
        sim_params.deployment_mode,
        0.0001,  # 密度（poissonモード時）
        sim_params.num_terminals,
        sim_params.area_size_m,
        sim_params.min_distance_m,
        sim_params.max_distance_m,
        sim_params.center_frequency_ghz * 1e9,  # GHz → Hz
        sim_params.path_loss_exponent,
        sim_params.reference_distance_m,
        sim_params.reference_path_loss_db
    )
    
    # シャドウイングパラメータ設定
    shadowing_params = ShadowingParameters(
        sim_params.shadowing_enabled,
        sim_params.shadowing_std_db,
        50.0,  # 相関距離
        0.5    # 相関係数
    )
    
    # 端末配置
    println("端末配置を生成中...")
    terminals = deploy_terminals(deployment_params, shadowing_params, sim_params.tx_power_dbm)
    
    # 端末配置情報を表示
    println("端末配置情報:")
    println("• 配置モード: $(deployment_params.deployment_mode)")
    println("• 生成された端末数: $(length(terminals))")
    if deployment_params.deployment_mode == "poisson"
        println("• 密度: $(deployment_params.lambda) 点/m²")
    else
        println("• 指定端末数: $(deployment_params.num_terminals)")
    end
    println("• エリアサイズ: $(deployment_params.area_size_m)m×$(deployment_params.area_size_m)m")
    println()
    
    for (i, terminal) in enumerate(terminals)
        println("• 端末$(i): 位置($(round(terminal.x_m, digits=1)), $(round(terminal.y_m, digits=1))) m, 距離$(round(terminal.distance_m, digits=1)) m")
        println("  - パスロス: $(round(terminal.path_loss_db, digits=1)) dB")
        println("  - シャドウイング: $(round(terminal.shadowing_db, digits=1)) dB")
        println("  - 総損失: $(round(terminal.total_loss_db, digits=1)) dB")
        println("  - 受信電力: $(round(terminal.rx_power_dbm, digits=1)) dBm")
    end
    println()
    
    # 代表端末を選択（最も近い端末）
    if !isempty(terminals)
        representative_terminal = terminals[argmin([t.distance_m for t in terminals])]
        println("代表端末（最も近い端末）:")
        println("• 位置: ($(round(representative_terminal.x_m, digits=1)), $(round(representative_terminal.y_m, digits=1))) m")
        println("• 距離: $(round(representative_terminal.distance_m, digits=1)) m")
        println("• パスロス: $(round(representative_terminal.path_loss_db, digits=1)) dB")
        println("• シャドウイング: $(round(representative_terminal.shadowing_db, digits=1)) dB")
        println("• 総損失: $(round(representative_terminal.total_loss_db, digits=1)) dB")
        println("• 受信電力: $(round(representative_terminal.rx_power_dbm, digits=1)) dBm")
    else
        println("警告: 端末が生成されませんでした。密度を上げるか、エリアサイズを大きくしてください。")
        return
    end
    println()
    
    # 送信側：周期的同期信号生成
    println("送信側：周期的同期信号を生成中...")
    time_axis, tx_signal, signal_count = generate_periodic_sync_signals(
        params, sim_params.signal_interval_ms, sim_params.total_duration_ms
    )
    
    # 代表端末での受信シミュレーション
    println("代表端末での受信シミュレーション中...")
    
    # 代表端末のパスロスパラメータを作成
    representative_path_loss_params = PathLossParameters(
        representative_terminal.distance_m,
        deployment_params.frequency_hz,
        deployment_params.path_loss_exponent,
        deployment_params.reference_distance_m,
        deployment_params.reference_path_loss_db
    )
    
    # 受信側：パスロス + シャドウイング + ノイズ付加
    println("受信側：パスロス、シャドウイング、ノイズを付加中...")
    rx_signal, actual_snr_db, path_loss_db = add_path_loss_and_noise(tx_signal, representative_path_loss_params, noise_params, shadowing_params, params.sampling_rate_hz)
    
    # 受信電力の計算
    rx_power = abs2.(rx_signal)
    
    # 受信電力ピーク検出
    println("受信電力ピークを検出中...")
    detected_peaks, peak_times, peak_powers, detection_rate = detect_sync_peaks(
        rx_power, time_axis, signal_count, sim_params.signal_interval_ms, sim_params.power_threshold_dbm
    )
    
    # 出力ディレクトリを作成（固定端末用）
    output_dir = "results_sync_simulation_fixed"
    mkpath(output_dir)
    
    # 実行時刻を取得
    execution_timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    println("実行時刻: $(execution_timestamp)")
    println()
    
    # 出力処理（受信電力CSVとパラメータCSVのみ）
    println("受信電力CSVファイルを生成中...")
    save_received_power_to_csv(time_axis, rx_power, peak_times, peak_powers, 
                              params, noise_params, terminals, representative_terminal,
                              signal_count, detection_rate, output_dir, shadowing_params)
    println("CSVファイルを '$(output_dir)' に保存しました。")
    println("• received_power_data_$(execution_timestamp).csv - 受信電力データ（mW, dBm）")
    println("• sync_simulation_parameters_$(execution_timestamp).csv - パラメータ情報")
    println()
    
    println("\n" * "="^60)
    println("プログラム完了")
    println("="^60)
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    # コード内でSNRを設定（上記のSNR設定部分を変更してください）
    main()  # デフォルトSNR=20dBで実行
end
