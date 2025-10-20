using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots, Dates

# ===== 信号パラメータ =====
struct SignalParameters
    duration_s::Float64       # 信号持続時間（s）
    center_frequency_hz::Float64  # 中心周波数（Hz）
    bandwidth_hz::Float64    # 帯域幅（Hz）
    sampling_rate_hz::Float64    # サンプリングレート（Hz）
    tx_power_dbm::Float64    # 送信電力（dBm）
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
    x_m::Float64            # X座標（m）
    y_m::Float64            # Y座標（m）
    distance_m::Float64     # 基地局からの距離（m）
    path_loss_db::Float64   # パスロス（dB）
    shadowing_db::Float64   # シャドウイング（dB）
    total_loss_db::Float64  # 総損失（パスロス+シャドウイング）（dB）
    rx_power_dbm::Float64   # 受信電力（dBm）
end

# ===== ノイズパラメータ =====
struct NoiseParameters
    noise_power_dbm::Float64  # ノイズ電力（dBm）
end

# ===== シミュレーションパラメータ（一括設定） =====
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

# ===== パスロス計算 =====
function calculate_path_loss(path_loss_params::PathLossParameters)
    distance_ratio = path_loss_params.distance_m / path_loss_params.reference_distance_m
    path_loss_db = path_loss_params.reference_path_loss_db + 
                   10 * path_loss_params.path_loss_exponent * log10(distance_ratio)
    return path_loss_db
end

# ===== シャドウイング計算 =====
function calculate_shadowing(shadowing_params::ShadowingParameters, distance_m::Float64, fixed_value::Bool=false)
    if !shadowing_params.enabled
        return 0.0
    end
    if fixed_value
        return shadowing_params.std_db
    else
        shadowing_db = randn() * shadowing_params.std_db
        return shadowing_db
    end
end

# ===== OFDM信号生成（QPSK変調） =====
function generate_ofdm_signal_with_qpsk(bandwidth_hz::Float64, duration_s::Float64, fs_hz::Float64, tx_power_dbm::Float64)
    # サンプル数とFFTサイズを計算
    N = Int(ceil(duration_s * fs_hz))
    N = max(N, 1)  # 最低1サンプル
    
    # 帯域幅に対応する周波数ビン数を計算
    bw_bins = Int(ceil(bandwidth_hz / fs_hz * N))
    bw_bins = max(bw_bins, 1)  # 最低1ビン
    
    # 周波数ドメインでQPSKシンボルを生成
    freq_domain = zeros(ComplexF64, N)
    
    # 正の周波数成分（中央から帯域幅分）
    start_bin = N ÷ 2 - bw_bins ÷ 2
    end_bin = N ÷ 2 + bw_bins ÷ 2
    
    for i in start_bin:end_bin
        if i > 0 && i <= N
            # QPSKシンボル生成（ランダム）
            real_part = (rand() > 0.5 ? 1.0 : -1.0) / sqrt(2)
            imag_part = (rand() > 0.5 ? 1.0 : -1.0) / sqrt(2)
            freq_domain[i] = real_part + 1im * imag_part
        end
    end
    
    # 負の周波数成分（共役対称）
    for i in 1:(N÷2)
        if N - i + 1 <= N && N - i + 1 > N ÷ 2
            freq_domain[N - i + 1] = conj(freq_domain[i])
        end
    end
    
    # IFFTで時間ドメイン信号を生成
    time_domain = ifft(freq_domain)
    
    # 送信電力に合わせてスケーリング（正規化の代わり）
    linear_power = 10^(tx_power_dbm / 10) * 1e-3  # dBm → W
    time_domain = time_domain * sqrt(linear_power)
    
    return time_domain
end

# ===== 周期的同期信号生成 =====
function generate_periodic_sync_signals(params::SignalParameters, interval_ms::Float64, total_duration_ms::Float64)
    signal_duration_ms = params.duration_s * 1000
    sampling_rate_hz = params.sampling_rate_hz
    total_samples = Int(ceil(total_duration_ms * sampling_rate_hz / 1000))
    
    # 信号送信タイミングを計算
    signal_times = Float64[]
    start_delay_ms = 10.0  # 開始遅延（ms）
    current_time = start_delay_ms
    while current_time + signal_duration_ms <= total_duration_ms
        push!(signal_times, current_time)
        current_time += interval_ms
    end
    
    # 送信信号配列を初期化
    tx_signal = zeros(ComplexF64, total_samples)
    time_axis = collect((0:total_samples-1) / sampling_rate_hz * 1000)  # ms単位
    
    # 各送信タイミングで信号を生成・配置
    signal_count = 0
    for signal_time in signal_times
        # 信号生成
        signal_samples = generate_ofdm_signal_with_qpsk(params.bandwidth_hz, params.duration_s, params.sampling_rate_hz, params.tx_power_dbm)
        
        # 信号配置
        start_sample = Int(ceil(signal_time * sampling_rate_hz / 1000)) + 1
        end_sample = min(start_sample + length(signal_samples) - 1, total_samples)
        
        if start_sample <= total_samples
            signal_length = end_sample - start_sample + 1
            tx_signal[start_sample:end_sample] = signal_samples[1:signal_length]
            signal_count += 1
            println("送信信号 $(signal_count): 時刻 $(signal_time) ms ～ $(signal_time + signal_duration_ms) ms")
        end
    end
    
    println("送信側パラメータ:")
    println("• 開始遅延: $(start_delay_ms) ms")
    println("• 信号間隔: $(interval_ms) ms")
    println("• 信号持続時間: $(signal_duration_ms) ms")
    println("• 総持続時間: $(total_duration_ms) ms")
    println("• 送信信号数: $(signal_count)")
    println("• 送信電力: $(params.tx_power_dbm) dBm")
    println()
    
    return time_axis, tx_signal, signal_count
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
    noise_power_linear = 10^(noise_power_dbm / 10) * 1e-3
    noise = sqrt(noise_power_linear) * (randn(N) + 1im * randn(N)) / sqrt(2)
    
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
function deploy_terminals(deployment_params::TerminalDeploymentParameters, shadowing_params::ShadowingParameters, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    
    if deployment_params.deployment_mode == "poisson"
        # ポアソン点過程モード（ランダム配置）
        lambda = deployment_params.lambda
        area_size = deployment_params.area_size_m
        num_points = rand(Poisson(lambda * area_size^2))
        println("ポアソン点過程モード: 密度$(lambda)点/m², 生成点数$(num_points)")
        
        for i in 1:num_points
            # 円内のランダムな位置を生成（面積均等分布）
            r_max = deployment_params.max_distance_m
            r_min = deployment_params.min_distance_m
            
            # 面積均等分布のための距離生成
            r = sqrt(r_min^2 + (r_max^2 - r_min^2) * rand())
            theta = 2 * π * rand()  # 0から2πのランダムな角度
            
            # 直交座標に変換
            x = r * cos(theta)
            y = r * sin(theta)
            distance = sqrt(x^2 + y^2)
            
            # 距離制限をチェック
            if distance >= deployment_params.min_distance_m && distance <= deployment_params.max_distance_m
                # パスロス計算
                path_loss_params = PathLossParameters(
                    distance, deployment_params.frequency_hz, deployment_params.path_loss_exponent,
                    deployment_params.reference_distance_m, deployment_params.reference_path_loss_db
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
                
                println("• 端末$(i): 位置($(round(x, digits=1)), $(round(y, digits=1))) m, 距離$(round(distance, digits=1)) m")
                println("  - パスロス: $(round(path_loss_db, digits=2)) dB")
                println("  - シャドウイング: $(round(shadowing_db, digits=2)) dB")
                println("  - 総損失: $(round(total_loss_db, digits=2)) dB")
                println("  - 受信電力: $(round(rx_power_dbm, digits=1)) dBm")
            end
        end
        
    elseif deployment_params.deployment_mode == "fixed"
        # 固定端末数モード
        println("固定端末数モード: 端末数$(deployment_params.num_terminals)（固定位置配置）")
        
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
        
        for i in 1:min(deployment_params.num_terminals, length(fixed_positions))
            x, y = fixed_positions[i]
            distance = sqrt(x^2 + y^2)
            
            # パスロス計算
            path_loss_params = PathLossParameters(
                distance, deployment_params.frequency_hz, deployment_params.path_loss_exponent,
                deployment_params.reference_distance_m, deployment_params.reference_path_loss_db
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
        expected_time = 5.0 + (i - 1) * interval_ms  # 開始遅延 + (i-1) * 間隔
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
            "tx_sampling_rate_mhz", "rx_sampling_rate_mhz", "tx_power_dbm", "noise_power_dbm",
            "shadowing_enabled", "shadowing_std_db", "deployment_mode", "num_terminals",
            "area_size_m", "min_distance_m", "max_distance_m", "path_loss_exponent",
            "reference_distance_m", "reference_path_loss_db", "signal_interval_ms",
            "total_duration_ms", "power_threshold_dbm", "detection_rate_percent"
        ],
        value = [
            params.duration_s * 1e6, params.center_frequency_hz / 1e9, params.bandwidth_hz / 1e6, 0.01,
            16.0, 0.02, params.tx_power_dbm, noise_params.noise_power_dbm,
            shadowing_params.enabled, shadowing_params.std_db, "fixed", length(terminals),
            100.0, 10.0, 500.0, 3.0, 1.0, 0.0, 20.0, 110.0, -130.0, detection_rate
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
        7.68,                # 送信側サンプリングレート（MHz）- 同期信号帯域幅の16倍
        0.25,                # 受信側サンプリングレート（MHz）- 端末受信帯域幅の2倍
        43.0,               # 送信電力（dBm）
        
        # 受信環境パラメータ
        true,               # シャドウイング有効/無効
        8.0,                # シャドウイング標準偏差（dB）
        
        # 端末配置パラメータ
        "fixed",          # 配置モード（ランダム配置）
        1,                  # 端末数
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

# ===== リサンプリング対応版メイン実行関数 =====
function main_resampling()
    #Random.seed!(1234)
    
    println("="^60)
    title_str = "同期信号受信シミュレーション（リサンプリング対応版）"
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
    
    # ===== サブキャリア間隔の算出 =====
    println("=== サブキャリア間隔の算出 ===")
    
    # 送信側のサブキャリア間隔
    tx_sampling_rate_hz = sim_params.tx_sampling_rate_mhz * 1e6
    signal_duration_s = sim_params.signal_duration_us * 1e-6
    tx_fft_size = Int(ceil(signal_duration_s * tx_sampling_rate_hz))
    tx_scs_hz = tx_sampling_rate_hz / tx_fft_size
    tx_scs_khz = tx_scs_hz / 1e3
    
    # 受信側のサブキャリア間隔
    rx_sampling_rate_hz = sim_params.rx_sampling_rate_mhz * 1e6
    rx_fft_size = Int(ceil(signal_duration_s * rx_sampling_rate_hz))
    rx_scs_hz = rx_sampling_rate_hz / rx_fft_size
    rx_scs_khz = rx_scs_hz / 1e3
    
    println("送信側:")
    println("• FFTサイズ: $(tx_fft_size)")
    println("• サブキャリア間隔: $(round(tx_scs_hz, digits=2)) Hz")
    println("• サブキャリア間隔: $(round(tx_scs_khz, digits=2)) kHz")
    println()
    
    println("受信側:")
    println("• FFTサイズ: $(rx_fft_size)")
    println("• サブキャリア間隔: $(round(rx_scs_hz, digits=2)) Hz")
    println("• サブキャリア間隔: $(round(rx_scs_khz, digits=2)) kHz")
    println()
    
    # 5G NR標準との比較
    println("5G NR標準との比較:")
    println("• 5G NR 120kHz SCS: 120.0 kHz")
    println("• 送信側SCS: $(round(tx_scs_khz, digits=2)) kHz (差: $(round(tx_scs_khz - 120, digits=2)) kHz)")
    println("• 受信側SCS: $(round(rx_scs_khz, digits=2)) kHz (差: $(round(rx_scs_khz - 120, digits=2)) kHz)")
    println()

    # ノイズパラメータ設定（端末受信帯域幅を使用）
    noise_figure_db = 5.0
    terminal_bandwidth_hz = sim_params.terminal_bandwidth_mhz * 1e6  # MHz → Hz
    fixed_noise_power_dbm = -174 + 10 * log10(terminal_bandwidth_hz) + noise_figure_db
    noise_params = NoiseParameters(fixed_noise_power_dbm)
    
    # 端末配置パラメータ設定
    deployment_params = TerminalDeploymentParameters(
        sim_params.deployment_mode,
        0.001,  # lambda (ポアソン点過程の密度)
        sim_params.num_terminals,
        sim_params.area_size_m,
        sim_params.min_distance_m,
        sim_params.max_distance_m,
        sim_params.center_frequency_ghz * 1e9,  # frequency_hz
        sim_params.path_loss_exponent,
        sim_params.reference_distance_m,
        sim_params.reference_path_loss_db
    )
    
    # シャドウイングパラメータ設定
    shadowing_params = ShadowingParameters(
        sim_params.shadowing_enabled,
        sim_params.shadowing_std_db,
        50.0,  # correlation_distance_m
        0.5    # correlation_coefficient
    )
    
    # 端末配置
    println("端末配置を生成中...")
    terminals = deploy_terminals(deployment_params, shadowing_params, sim_params.tx_power_dbm)
    
    if isempty(terminals)
        println("警告: 端末が生成されませんでした。密度を上げるか、エリアサイズを大きくしてください。")
        return
    end
    println()

    # ===== 3. 高いレートで「理想的な送信信号」を生成 =====
    # 送信側の高レート設定でSignalParametersを作成
    params_tx = SignalParameters(
        sim_params.signal_duration_us * 1e-6,
        sim_params.center_frequency_ghz * 1e9,
        sim_params.signal_bandwidth_mhz * 1e6,
        sim_params.tx_sampling_rate_mhz * 1e6, # ★ 送信側の高レートを使用
        sim_params.tx_power_dbm
    )
    
    println("送信側：高レートで理想信号を生成中...")
    time_axis_tx, tx_signal_high_rate, signal_count = generate_periodic_sync_signals(
        params_tx, sim_params.signal_interval_ms, sim_params.total_duration_ms
    )
    println("生成された高レート信号のサンプル数: $(length(tx_signal_high_rate))")
    println()

    # ===== 4. `resample`関数でダウンサンプリング（端末のADCを模擬）=====
    println("受信側：ダウンサンプリングを実行中...")
    rate_ratio = sim_params.rx_sampling_rate_mhz / sim_params.tx_sampling_rate_mhz # レート比を計算 (例: 0.02 / 16.0)
    
    # ダウンサンプリング実行（resample関数が自動的にアンチエイリアシングフィルタを適用）
    tx_signal_low_rate = resample(tx_signal_high_rate, rate_ratio)
    
    println("ダウンサンプリング後の信号のサンプル数: $(length(tx_signal_low_rate))")
    println("端末受信帯域幅: $(sim_params.terminal_bandwidth_mhz) MHz")
    println("注意: resample関数が自動的に帯域制限フィルタを適用済み")
    println()

    # ===== 5. ダウンサンプリング後の信号に対してノイズ付加と各種処理を行う =====
    # 時間軸も受信側のレートに合わせて再計算
    time_axis_rx = collect((0:length(tx_signal_low_rate)-1) / (sim_params.rx_sampling_rate_mhz * 1e6) * 1000) # ms単位

    # 代表端末での受信シミュレーション
    println("代表端末での受信シミュレーション中...")
    
    # 代表端末のパスロスパラメータを作成（最も近い端末を選択）
    representative_terminal = terminals[argmin([t.distance_m for t in terminals])]
    representative_path_loss_params = PathLossParameters(
        representative_terminal.distance_m,
        sim_params.center_frequency_ghz * 1e9,  # frequency_hz
        sim_params.path_loss_exponent,
        sim_params.reference_distance_m,
        sim_params.reference_path_loss_db
    )
    
    println("受信側：パスロス、シャドウイング、ノイズを付加中...")
    # ★ ダウンサンプリング後の信号に対してノイズを付加（端末の固定シャドウイング値を使用）
    rx_signal, actual_snr_db, path_loss_db = add_path_loss_and_noise(
        tx_signal_low_rate, representative_path_loss_params, noise_params, shadowing_params, sim_params.rx_sampling_rate_mhz * 1e6, 
        representative_terminal.shadowing_db
    )
    
    rx_power = abs2.(rx_signal)
    
    println("受信電力ピークを検出中...")
    detected_peaks, peak_times, peak_powers, detection_rate = detect_sync_peaks(
        rx_power, time_axis_rx, signal_count, sim_params.signal_interval_ms, sim_params.power_threshold_dbm
    )
    
    # 出力ディレクトリを作成
    output_dir = "results_QPSK_simulation_resampling"
    mkpath(output_dir)
    
    # 実行時刻を取得
    execution_timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    # 出力処理
    println("実行時刻: $(execution_timestamp)")
    println()
    
    println("受信電力CSVファイルを生成中...")
    save_received_power_to_csv(time_axis_rx, rx_power, peak_times, peak_powers, 
                              params_tx, noise_params, terminals, representative_terminal,
                              signal_count, detection_rate, output_dir, shadowing_params)
    println("CSVファイルを '$(output_dir)' に保存しました。")
    println("• received_power_data_$(execution_timestamp).csv - 受信電力データ（mW, dBm）")
    println("• QPSK_simulation_parameters_$(execution_timestamp).csv - パラメータ情報")
    println()
    
    println("="^60)
    println("プログラム完了")
    println("="^60)
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    main_resampling()
end
