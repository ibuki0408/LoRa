# ===== メインシミュレーション実行ファイル =====

using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots, Dates

# 各機能モジュールを読み込み
include("signal_generation.jl")
include("path_loss.jl")
include("shadowing.jl")
include("noise_generation.jl")
include("terminal_deployment.jl")
# include("signal_visualization.jl")  # コメントアウト

# ===== シミュレーションパラメータ =====
struct SimulationParameters
    # 信号パラメータ
    signal_duration_us::Float64      # 信号持続時間（μs）
    center_frequency_ghz::Float64     # 中心周波数（GHz）
    signal_bandwidth_mhz::Float64    # 同期信号帯域幅（MHz）
    terminal_bandwidth_mhz::Float64  # 端末受信帯域幅（MHz）
    tx_sampling_rate_mhz::Float64    # 送信側サンプリングレート（MHz）
    rx_sampling_rate_mhz::Float64    # 受信側サンプリングレート（MHz）
    tx_power_dbm::Float64            # 送信電力（dBm）
    
    # 受信環境パラメータ
    shadowing_enabled::Bool          # シャドウイング有効/無効
    shadowing_std_db::Float64        # シャドウイング標準偏差（dB）
    
    # 端末配置パラメータ
    deployment_mode::String          # 配置モード: "poisson" または "fixed"
    num_terminals::Int               # 端末数
    area_size_m::Float64            # エリアサイズ（m）
    min_distance_m::Float64         # 最小距離（m）
    max_distance_m::Float64         # 最大距離（m）
    path_loss_exponent::Float64      # パスロス指数
    reference_distance_m::Float64   # 参照距離（m）
    reference_path_loss_db::Float64 # 参照距離でのパスロス（dB）
    
    # シミュレーション制御パラメータ
    signal_interval_ms::Float64      # 信号間隔（ms）
    total_duration_ms::Float64       # 総持続時間（ms）
    power_threshold_dbm::Float64     # ピーク検出閾値（dBm）
end

# ===== 端末情報 =====
struct TerminalInfo
    x_m::Float64            # X座標（m）
    y_m::Float64            # Y座標（m）
    distance_m::Float64     # 基地局からの距離（m）
    path_loss_db::Float64   # パスロス（dB）
    shadowing_db::Float64   # シャドウイング（dB）
    total_loss_db::Float64  # 総損失（パスロス+シャドウイング）（dB）
    rx_power_dbm::Float64   # 受信電力（dBm）
end

# ===== パスロス・シャドウイング・ノイズ付加 =====
function add_path_loss_and_noise(tx_signal::Vector{ComplexF64}, path_loss_params::PathLossParameters, 
                                noise_params::NoiseParameters, shadowing_params::ShadowingParameters, sampling_rate_hz::Float64, 
                                terminal_shadowing_db::Float64)
    N = length(tx_signal)
    
    # パスロス計算
    path_loss_db = calculate_path_loss(path_loss_params)
    path_loss_linear = 10^(-path_loss_db / 10)
    
    # シャドウイング計算（端末の固定値を使用）
    shadowing_db = terminal_shadowing_db
    shadowing_linear = 10^(-shadowing_db / 10)
    
    # 総損失を適用
    total_loss_linear = path_loss_linear * shadowing_linear
    rx_signal_with_path_loss = tx_signal * sqrt(total_loss_linear)
    
    # 信号電力計算（パスロス後）
    signal_power_linear = mean(abs2.(rx_signal_with_path_loss))
    signal_power_dbm = 10 * log10(signal_power_linear * 1000)
    
    # ノイズ生成（固定ノイズ床）
    noise_power_dbm = noise_params.noise_power_dbm  # 固定値
    noise = generate_awgn_noise(noise_power_dbm, N)
    
    # ノイズ付加後の受信信号
    rx_signal_with_noise = rx_signal_with_path_loss + noise
    
    # 受信電力計算（ノイズ後）
    rx_power_linear = mean(abs2.(rx_signal_with_noise))
    rx_power_dbm = 10 * log10(rx_power_linear * 1000)
    
    # 実際のSNRを計算
    actual_noise_power = mean(abs2.(noise))
    actual_snr_db = 10 * log10(signal_power_linear / actual_noise_power)
    
    println("受信側パスロス・シャドウイング・AWGNパラメータ:")
    println("• 距離: $(round(path_loss_params.distance_m, digits=2)) m")
    println("• パスロス: $(round(path_loss_db, digits=2)) dB")
    println("• シャドウイング: $(round(shadowing_db, digits=2)) dB（端末固定値）")
    println("• 総損失: $(round(path_loss_db + shadowing_db, digits=2)) dB")
    println("• 受信信号電力（ノイズなし）: $(round(signal_power_dbm, digits=2)) dBm")
    println("• 受信信号電力（ノイズあり）: $(round(rx_power_dbm, digits=2)) dBm")
    println("• 固定ノイズ電力: $(round(noise_power_dbm, digits=2)) dBm")
    println("• 実際のSNR: $(round(actual_snr_db, digits=2)) dB")
    println()
    
    return rx_signal_with_noise, actual_snr_db, path_loss_db
end

# ===== 端末配置 =====
function deploy_terminals(deployment_mode::String, num_terminals::Int, area_size_m::Float64, 
                         min_distance_m::Float64, max_distance_m::Float64, frequency_hz::Float64,
                         path_loss_exponent::Float64, reference_distance_m::Float64, reference_path_loss_db::Float64,
                         shadowing_params::ShadowingParameters, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    
    if deployment_mode == "fixed"
        # 固定端末数モード
        println("固定端末数モード: 端末数$(num_terminals)（固定位置配置）")
        
        # 固定位置を定義
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
        
        for i in 1:min(num_terminals, length(fixed_positions))
            x, y = fixed_positions[i]
            distance = sqrt(x^2 + y^2)
            
            # パスロス計算
            path_loss_params = PathLossParameters(
                distance, frequency_hz, path_loss_exponent,
                reference_distance_m, reference_path_loss_db
            )
            path_loss_db = calculate_path_loss(path_loss_params)
            
            # シャドウイング計算
            shadowing_db = calculate_shadowing(shadowing_params, distance, false)
            
            # 総損失
            total_loss_db = path_loss_db + shadowing_db
            
            # 受信電力計算
            rx_power_dbm = tx_power_dbm - total_loss_db
            
            terminal = TerminalInfo(x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm)
            push!(terminals, terminal)
            
            println("• 端末$(i): 位置($(x), $(y)) m, 距離$(round(distance, digits=1)) m")
            println("  - パスロス: $(round(path_loss_db, digits=2)) dB")
            println("  - シャドウイング: $(round(shadowing_db, digits=2)) dB")
            println("  - 総損失: $(round(total_loss_db, digits=2)) dB")
            println("  - 受信電力: $(round(rx_power_dbm, digits=1)) dBm")
        end
    end
    
    return terminals
end

# ===== 同期ピーク検出 =====
function detect_sync_peaks(rx_power::Vector{Float64}, time_axis::Vector{Float64}, 
                          signal_count::Int, interval_ms::Float64, power_threshold_db::Float64)
    detected_peaks = Int[]
    peak_times = Float64[]
    peak_powers = Float64[]
    
    # 動的閾値計算
    avg_power = mean(rx_power)
    dynamic_threshold = max(3 * avg_power, 10^(power_threshold_db / 10) * 1e-3)
    
    # 各信号タイミングでピーク検出
    for i in 1:signal_count
        expected_time = 10.0 + (i - 1) * interval_ms  # 開始遅延 + (i-1) * 間隔
        time_tolerance = interval_ms * 0.1  # 10%の許容誤差
        
        # 期待時刻周辺のサンプルを検索
        start_idx = findfirst(t -> t >= expected_time - time_tolerance, time_axis)
        end_idx = findlast(t -> t <= expected_time + time_tolerance, time_axis)
        
        if start_idx !== nothing && end_idx !== nothing
            search_range = start_idx:end_idx
            max_power_idx = argmax(rx_power[search_range])
            actual_idx = start_idx + max_power_idx - 1
            
            if rx_power[actual_idx] > dynamic_threshold
                push!(detected_peaks, actual_idx)
                push!(peak_times, time_axis[actual_idx])
                push!(peak_powers, rx_power[actual_idx])
            end
        end
    end
    
    detection_rate = length(detected_peaks) / signal_count * 100
    
    println("受信電力ピーク検出結果:")
    println("• 検出閾値: $(round(dynamic_threshold * 1000, digits=1)) mW ($(round(10*log10(dynamic_threshold*1000), digits=1)) dBm)")
    println("• 検出ピーク数: $(length(detected_peaks))/$(signal_count)")
    println("• 検出率: $(round(detection_rate, digits=1))%")
    if !isempty(peak_powers)
        println("• ピーク電力範囲: $(round(minimum(peak_powers)*1000, digits=1)) - $(round(maximum(peak_powers)*1000, digits=1)) mW")
        println("• ピーク電力範囲: $(round(10*log10(minimum(peak_powers)*1000), digits=1)) - $(round(10*log10(maximum(peak_powers)*1000), digits=1)) dBm")
    end
    println()
    
    return detected_peaks, peak_times, peak_powers, detection_rate
end

# ===== CSV保存 =====
function save_received_power_to_csv(time_axis::Vector{Float64}, rx_power::Vector{Float64}, 
                                   peak_times::Vector{Float64}, peak_powers::Vector{Float64},
                                   params::SignalParameters, noise_params::NoiseParameters, 
                                   terminals::Vector{TerminalInfo}, representative_terminal::TerminalInfo,
                                   signal_count::Int, detection_rate::Float64, output_dir::String, 
                                   shadowing_params::ShadowingParameters)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    # 受信電力データ
    min_length = min(length(time_axis), length(rx_power))
    power_data = DataFrame(
        time_ms = time_axis[1:min_length],
        power_mw = rx_power[1:min_length] * 1000,
        power_dbm = 10*log10.(rx_power[1:min_length] * 1000)
    )
    CSV.write("$(output_dir)/received_power_data_$(timestamp).csv", power_data)
    
    # パラメータデータ
    params_data = DataFrame(
        parameter = [
            "signal_duration_us", "center_frequency_ghz", "signal_bandwidth_mhz", "terminal_bandwidth_mhz",
            "tx_sampling_rate_mhz", "rx_sampling_rate_mhz", "tx_power_dbm", "snr_db",
            "shadowing_enabled", "shadowing_std_db", "deployment_mode", "num_terminals",
            "area_size_m", "min_distance_m", "max_distance_m", "path_loss_exponent",
            "reference_distance_m", "reference_path_loss_db", "signal_interval_ms",
            "total_duration_ms", "power_threshold_dbm", "detection_rate_percent"
        ],
        value = [
            params.duration_s * 1e6, params.center_frequency_hz / 1e9, params.bandwidth_hz / 1e6, 0.125,
            61.44, 0.25, params.tx_power_dbm, noise_params.snr_db,
            shadowing_params.enabled, shadowing_params.std_db, "fixed", length(terminals),
            100.0, 10.0, 500.0, 3.0, 1.0, 0.0, 20.0, 110.0, -100.0, detection_rate
        ]
    )
    CSV.write("$(output_dir)/QPSK_simulation_parameters_$(timestamp).csv", params_data)
end

# ===== シミュレーションパラメータ作成 =====
function create_simulation_parameters()
    return SimulationParameters(
        # 信号パラメータ
        66.67,              # 信号持続時間（μs）
        4.7,                # 中心周波数（GHz）
        3.6,                # 同期信号帯域幅（MHz）- 広帯域信号
        0.125,                # 端末受信帯域幅（MHz）- より現実的な帯域
        3.84,                # 送信側サンプリングレート（MHz）- 同期信号帯域幅の16倍
        0.25,                # 受信側サンプリングレート（MHz）- 端末受信帯域幅の2倍
        43.0,               # 送信電力（dBm）
        
        # 受信環境パラメータ
        true,               # シャドウイング有効/無効
        8.0,                # シャドウイング標準偏差（dB）
        
        # 端末配置パラメータ
        "fixed",          # 配置モード（ポアソン点過程）
        1,                  # 端末数（固定）
        100.0,              # エリアサイズ（m）
        10.0,               # 最小距離（m）
        500.0,              # 最大距離（m）
        3.0,                # パスロス指数
        1.0,                # 参照距離（m）
        0.0,                # 参照距離でのパスロス（dB）
        
        # シミュレーション制御パラメータ
        20.0,               # 信号間隔（ms）
        110.0,              # 総持続時間（ms）
        -100.0              # ピーク検出閾値（dBm）
    )
end

# ===== メイン実行関数 =====
function main_simulation()
    Random.seed!(1234)  # ランダムシードを固定
    
    println("="^60)
    title_str = "分割されたシミュレーション実行"
    println(title_str)
    println("="^60)

    # ===== 1. パラメータ設定 =====
    sim_params = create_simulation_parameters()
    
    println("シミュレーションパラメータ:")
    println("• 同期信号帯域幅: $(sim_params.signal_bandwidth_mhz) MHz")
    println("• 端末受信帯域幅: $(sim_params.terminal_bandwidth_mhz) MHz")
    println("• 送信側サンプリングレート: $(sim_params.tx_sampling_rate_mhz) MHz (高レート)")
    println("• 受信側サンプリングレート: $(sim_params.rx_sampling_rate_mhz) MHz (低レート)")
    println()
    
    # ===== 2. ノイズパラメータ設定 =====
    noise_figure_db = 5.0
    terminal_bandwidth_hz = sim_params.terminal_bandwidth_mhz * 1e6  # MHz → Hz
    fixed_noise_power_dbm = -174 + 10 * log10(terminal_bandwidth_hz) + noise_figure_db
    noise_params = NoiseParameters(fixed_noise_power_dbm, 100.0)  # 固定SNR値
    
    # ===== 3. シャドウイングパラメータ設定 =====
    shadowing_params = ShadowingParameters(
        sim_params.shadowing_enabled,
        sim_params.shadowing_std_db,
        50.0,  # correlation_distance_m
        0.5    # correlation_coefficient
    )
    
    # ===== 4. 端末配置 =====
    println("端末配置を生成中...")
    deployment_params = TerminalDeploymentParameters(
        sim_params.deployment_mode,
        1.0,  # lambda (ポアソン点過程の密度)
        sim_params.num_terminals,  # 端末数
        sim_params.area_size_m,
        sim_params.min_distance_m,
        sim_params.max_distance_m,
        sim_params.center_frequency_ghz * 1e9,
        sim_params.path_loss_exponent,
        sim_params.reference_distance_m,
        sim_params.reference_path_loss_db
    )
    
    terminals = deploy_terminals(sim_params.deployment_mode, sim_params.num_terminals, sim_params.area_size_m, 
                               sim_params.min_distance_m, sim_params.max_distance_m, 
                               sim_params.center_frequency_ghz * 1e9, sim_params.path_loss_exponent,
                               sim_params.reference_distance_m, sim_params.reference_path_loss_db,
                               shadowing_params, sim_params.tx_power_dbm)
    
    if isempty(terminals)
        println("警告: 端末が生成されませんでした。密度を上げるか、エリアサイズを大きくしてください。")
        return
    end
    
    # 端末配置統計分析
    analyze_terminal_deployment(terminals)
    println()

    # ===== 5. 高いレートで「理想的な送信信号」を生成 =====
    params_tx = SignalParameters(
        sim_params.signal_duration_us * 1e-6,
        sim_params.center_frequency_ghz * 1e9,
        sim_params.signal_bandwidth_mhz * 1e6,
        sim_params.tx_sampling_rate_mhz * 1e6,
        sim_params.tx_power_dbm
    )
    
    println("送信側：高レートで理想信号を生成中...")
    time_axis_tx, tx_signal_high_rate, signal_count = generate_periodic_sync_signals(
        params_tx, sim_params.signal_interval_ms, sim_params.total_duration_ms
    )
    println("生成された高レート信号のサンプル数: $(length(tx_signal_high_rate))")
    
    # ===== 送信信号の詳細情報出力 =====
    println("=== 送信信号の詳細情報 ===")
    signal_amplitude = abs.(tx_signal_high_rate)
    signal_power = abs2.(tx_signal_high_rate)
    
    println("• 時間軸情報:")
    println("  - サンプル数: $(length(tx_signal_high_rate))")
    println("  - 時間範囲: $(round(minimum(time_axis_tx), digits=3)) - $(round(maximum(time_axis_tx), digits=3)) ms")
    println("  - サンプリング間隔: $(round((time_axis_tx[2] - time_axis_tx[1]), digits=6)) ms")
    println("  - サンプリングレート: $(sim_params.tx_sampling_rate_mhz) MHz")
    
    println("• 振幅情報:")
    println("  - 最大振幅: $(round(maximum(signal_amplitude), digits=6))")
    println("  - 最小振幅: $(round(minimum(signal_amplitude), digits=6))")
    println("  - 平均振幅: $(round(mean(signal_amplitude), digits=6))")
    println("  - 振幅標準偏差: $(round(std(signal_amplitude), digits=6))")
    
    println("• 電力情報:")
    signal_power_mw = signal_power * 1000
    signal_power_dbm = 10 * log10.(signal_power_mw)
    println("  - 最大電力: $(round(maximum(signal_power_mw), digits=6)) mW ($(round(maximum(signal_power_dbm), digits=2)) dBm)")
    println("  - 最小電力: $(round(minimum(signal_power_mw), digits=6)) mW ($(round(minimum(signal_power_dbm), digits=2)) dBm)")
    println("  - 平均電力: $(round(mean(signal_power_mw), digits=6)) mW ($(round(mean(signal_power_dbm), digits=2)) dBm)")
    
    # 信号の実部・虚部情報
    signal_real = real.(tx_signal_high_rate)
    signal_imag = imag.(tx_signal_high_rate)
    println("• 複素信号情報:")
    println("  - 実部範囲: $(round(minimum(signal_real), digits=6)) - $(round(maximum(signal_real), digits=6))")
    println("  - 虚部範囲: $(round(minimum(signal_imag), digits=6)) - $(round(maximum(signal_imag), digits=6))")
    println("  - 実部平均: $(round(mean(signal_real), digits=6))")
    println("  - 虚部平均: $(round(mean(signal_imag), digits=6))")
    
    println()

    # ===== 6. ダウンサンプリング =====
    println("受信側：ダウンサンプリングを実行中...")
    rate_ratio = sim_params.rx_sampling_rate_mhz / sim_params.tx_sampling_rate_mhz
    tx_signal_low_rate = resample(tx_signal_high_rate, rate_ratio)
    println("ダウンサンプリング後の信号のサンプル数: $(length(tx_signal_low_rate))")
    
    # ===== ダウンサンプリング後信号の詳細情報出力 =====
    println("=== ダウンサンプリング後信号の詳細情報 ===")
    low_rate_amplitude = abs.(tx_signal_low_rate)
    low_rate_power = abs2.(tx_signal_low_rate)
    
    println("• 時間軸情報:")
    println("  - サンプル数: $(length(tx_signal_low_rate))")
    println("  - サンプリングレート: $(sim_params.rx_sampling_rate_mhz) MHz")
    println("  - ダウンサンプリング比: $(round(rate_ratio, digits=6))")
    
    println("• 振幅情報:")
    println("  - 最大振幅: $(round(maximum(low_rate_amplitude), digits=6))")
    println("  - 最小振幅: $(round(minimum(low_rate_amplitude), digits=6))")
    println("  - 平均振幅: $(round(mean(low_rate_amplitude), digits=6))")
    println("  - 振幅標準偏差: $(round(std(low_rate_amplitude), digits=6))")
    
    println("• 電力情報:")
    low_rate_power_mw = low_rate_power * 1000
    low_rate_power_dbm = 10 * log10.(low_rate_power_mw)
    println("  - 最大電力: $(round(maximum(low_rate_power_mw), digits=6)) mW ($(round(maximum(low_rate_power_dbm), digits=2)) dBm)")
    println("  - 最小電力: $(round(minimum(low_rate_power_mw), digits=6)) mW ($(round(minimum(low_rate_power_mw), digits=2)) dBm)")
    println("  - 平均電力: $(round(mean(low_rate_power_mw), digits=6)) mW ($(round(mean(low_rate_power_dbm), digits=2)) dBm)")
    
    # ダウンサンプリング前後の比較
    println("• ダウンサンプリング前後の比較:")
    println("  - 振幅比: $(round(mean(low_rate_amplitude) / mean(signal_amplitude), digits=6))")
    println("  - 電力比: $(round(mean(low_rate_power_mw) / mean(signal_power_mw), digits=6))")
    println("  - サンプル数比: $(round(length(tx_signal_low_rate) / length(tx_signal_high_rate), digits=6))")
    
    println()

    # ===== 7. 受信シミュレーション =====
    time_axis_rx = collect((0:length(tx_signal_low_rate)-1) / (sim_params.rx_sampling_rate_mhz * 1e6) * 1000)
    
    representative_terminal = terminals[argmin([t.distance_m for t in terminals])]
    representative_path_loss_params = PathLossParameters(
        representative_terminal.distance_m,
        sim_params.center_frequency_ghz * 1e9,
        sim_params.path_loss_exponent,
        sim_params.reference_distance_m,
        sim_params.reference_path_loss_db
    )
    
    println("受信側：パスロス、シャドウイング、ノイズを付加中...")
    rx_signal, actual_snr_db, path_loss_db = add_path_loss_and_noise(
        tx_signal_low_rate, representative_path_loss_params, noise_params, shadowing_params, 
        sim_params.rx_sampling_rate_mhz * 1e6, representative_terminal.shadowing_db
    )
    
    rx_power = abs2.(rx_signal)
    
    # ===== 8. ピーク検出 =====
    println("受信電力ピークを検出中...")
    detected_peaks, peak_times, peak_powers, detection_rate = detect_sync_peaks(
        rx_power, time_axis_rx, signal_count, sim_params.signal_interval_ms, sim_params.power_threshold_dbm
    )
    
    # ===== 9. 信号可視化（無効化） =====
    # println("信号可視化を実行中...")
    # visualize_signal_analysis(tx_signal_high_rate, time_axis_tx, 
    #                          tx_signal_low_rate, time_axis_rx,
    #                          rx_power, peak_times, peak_powers,
    #                          "results_modular_simulation")
    
    # ===== 10. 結果保存 =====
    output_dir = "results_modular_simulation"
    mkpath(output_dir)
    
    execution_timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    println("受信電力CSVファイルを生成中...")
    save_received_power_to_csv(time_axis_rx, rx_power, peak_times, peak_powers, 
                              params_tx, noise_params, terminals, representative_terminal,
                              signal_count, detection_rate, output_dir, shadowing_params)
    
    # ===== 送信信号データのCSV保存 =====
    println("送信信号データをCSVファイルに保存中...")
    
    # 高レート送信信号データ
    tx_high_rate_data = DataFrame(
        time_ms = time_axis_tx,
        real_part = real.(tx_signal_high_rate),
        imag_part = imag.(tx_signal_high_rate),
        amplitude = abs.(tx_signal_high_rate),
        power_mw = abs2.(tx_signal_high_rate) * 1000,
        power_dbm = 10*log10.(abs2.(tx_signal_high_rate) * 1000)
    )
    CSV.write("$(output_dir)/tx_high_rate_signal_$(execution_timestamp).csv", tx_high_rate_data)
    
    # 低レート送信信号データ
    time_axis_low_rate = collect((0:length(tx_signal_low_rate)-1) / (sim_params.rx_sampling_rate_mhz * 1e6) * 1000)
    tx_low_rate_data = DataFrame(
        time_ms = time_axis_low_rate,
        real_part = real.(tx_signal_low_rate),
        imag_part = imag.(tx_signal_low_rate),
        amplitude = abs.(tx_signal_low_rate),
        power_mw = abs2.(tx_signal_low_rate) * 1000,
        power_dbm = 10*log10.(abs2.(tx_signal_low_rate) * 1000)
    )
    CSV.write("$(output_dir)/tx_low_rate_signal_$(execution_timestamp).csv", tx_low_rate_data)
    
    println("CSVファイルを '$(output_dir)' に保存しました。")
    println("• received_power_data_$(execution_timestamp).csv - 受信電力データ")
    println("• tx_high_rate_signal_$(execution_timestamp).csv - 高レート送信信号データ")
    println("• tx_low_rate_signal_$(execution_timestamp).csv - 低レート送信信号データ")
    println()
    
    println("="^60)
    println("プログラム完了")
    println("="^60)
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    main_simulation()
end
