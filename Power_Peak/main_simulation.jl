# ===== メインシミュレーション実行ファイル =====

using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots, Dates

# 各機能モジュールを読み込み
include("modules/signal_generation.jl")
include("modules/path_loss.jl")
include("modules/shadowing.jl")
include("modules/noise_generation.jl")
include("modules/terminal_deployment.jl")
include("archive/adaptive_threshold.jl")
include("modules/local_clock.jl")
include("archive/sync_observation.jl")
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
    area_size_m::Float64             # エリアサイズ（m）
    min_distance_m::Float64         # 最小距離（m）
    max_distance_m::Float64         # 最大距離（m）
    path_loss_exponent::Float64      # パスロス指数
    reference_distance_m::Float64   # 参照距離（m）
    reference_path_loss_db::Float64 # 参照距離でのパスロス（dB）
    
    # シミュレーション制御パラメータ
    signal_interval_ms::Float64      # 信号間隔（ms）
    total_duration_ms::Float64       # 総持続時間（ms）
    power_threshold_dbm::Float64     # ピーク検出閾値（dBm）
    
    # クロックパラメータ
    clock_precision::String          # クロック精度: "high", "standard", "low"
    temperature_variation_c::Float64 # 温度変動（°C）
    sync_accuracy_ns::Float64       # 同期精度（ns）
    representative_selection::String # 代表端末選択: "nearest" | "median" | "farthest"
end

# ===== ノイズフロア推定（事前窓） =====
function estimate_noise_floor_dbm(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}; pre_signal_end_ms::Float64=9.0)
    idxs = findall(t -> t < pre_signal_end_ms, time_axis_ms)
    if isempty(idxs)
        # フォールバック: 全区間の中央値
        noise_power_w = median(rx_power)
        return 10 * log10(noise_power_w * 1000)
    end
    segment = rx_power[idxs]
    # ロバストに中央値を使用
    noise_power_w = median(segment)
    return 10 * log10(noise_power_w * 1000)
end

# ===== シンプルな閾値超過検出（期待時刻不要） =====
"""
受信電力が閾値を超えている全ての時間区間を記録
"""
function detect_threshold_crossings(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}, threshold_w::Float64)
    crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[],      # 閾値上抜け時刻
        :end_times => Float64[],        # 閾値下抜け時刻
        :peak_times => Float64[],       # ピーク時刻
        :peak_powers_w => Float64[],    # ピーク電力（W）
        :peak_powers_dbm => Float64[]   # ピーク電力（dBm）
    )
    
    is_above_threshold = false
    segment_start_idx = 0
    segment_peak_idx = 0
    segment_peak_power = 0.0
    
    for i in 2:length(rx_power)
        p_prev = rx_power[i-1]
        p_curr = rx_power[i]
        t_prev = time_axis_ms[i-1]
        t_curr = time_axis_ms[i]
        
        # 上抜け検出
        if !is_above_threshold && p_prev <= threshold_w && p_curr > threshold_w
            # 線形補間で正確な上抜け時刻を計算
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            
            push!(crossings[:start_times], crossing_time)
            is_above_threshold = true
            segment_start_idx = i
            segment_peak_idx = i
            segment_peak_power = p_curr
        end
        
        # 閾値超過区間内でピークを追跡
        if is_above_threshold
            if p_curr > segment_peak_power
                segment_peak_idx = i
                segment_peak_power = p_curr
            end
        end
        
        # 下抜け検出
        if is_above_threshold && p_prev > threshold_w && p_curr <= threshold_w
            # 線形補間で正確な下抜け時刻を計算
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            
            push!(crossings[:end_times], crossing_time)
            
            # この区間のピークを記録
            peak_time = time_axis_ms[segment_peak_idx]
            push!(crossings[:peak_times], peak_time)
            push!(crossings[:peak_powers_w], segment_peak_power)
            push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
            
            is_above_threshold = false
        end
    end
    
    # 最後まで閾値を超えていた場合の処理
    if is_above_threshold
        # 最後の時刻を下抜け時刻として記録（実際には信号終了）
        if !isempty(time_axis_ms)
            push!(crossings[:end_times], time_axis_ms[end])
            peak_time = time_axis_ms[segment_peak_idx]
            push!(crossings[:peak_times], peak_time)
            push!(crossings[:peak_powers_w], segment_peak_power)
            push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
        end
    end
    
    return crossings
end

# ===== デバウンス処理フィルタ（★ここが新しい関数） =====
"""
検出された区間が、直前の有効な区間から近すぎる場合に除外する。
ピーク時刻（peak_times）を基準に判断する。
"""
function filter_by_debounce(crossings::Dict{Symbol, Vector{Float64}}, debounce_time_ms::Float64)
    
    # 処理対象がない場合はそのまま返す
    if isempty(crossings[:peak_times])
        return crossings
    end

    # フィルタリング後の結果を格納する新しい辞書
    filtered_crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[],
        :end_times => Float64[],
        :peak_times => Float64[],
        :peak_powers_w => Float64[],
        :peak_powers_dbm => Float64[]
    )

    # 最後に「有効」として採用された信号のピーク時刻
    # -Inf で初期化し、最初の信号は必ず採用されるようにする
    last_accepted_peak_time_ms = -Inf

    # 検出されたすべての区間を順にチェック
    for i in 1:length(crossings[:peak_times])
        current_peak_time = crossings[:peak_times][i]

        # 現在のピーク時刻が、最後に採用したピーク時刻 + デバウンス時間 よりも後か？
        if current_peak_time >= last_accepted_peak_time_ms + debounce_time_ms
            # 【有効な検出】
            # 採用リストに追加
            push!(filtered_crossings[:start_times], crossings[:start_times][i])
            push!(filtered_crossings[:end_times], crossings[:end_times][i])
            push!(filtered_crossings[:peak_times], current_peak_time)
            push!(filtered_crossings[:peak_powers_w], crossings[:peak_powers_w][i])
            push!(filtered_crossings[:peak_powers_dbm], crossings[:peak_powers_dbm][i])

            # 「最後に採用した時刻」を更新する
            last_accepted_peak_time_ms = current_peak_time
        else
            # 【無効な検出】
            # デバウンス期間内のため、この区間は無視する
        end
    end

    return filtered_crossings
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

# ===== 端末配置（terminal_deployment.jlを使用） =====
# 端末配置は terminal_deployment.jl の deploy_terminals() 関数を使用

# ===== 同期ピーク検出（固定閾値） =====
function detect_sync_peaks(rx_power::Vector{Float64}, time_axis::Vector{Float64}, 
                          signal_count::Int, interval_ms::Float64, power_threshold_db::Float64)
    detected_peaks = Int[]
    detection_times = Float64[]  # 閾値上抜け時刻（ms）
    peak_powers = Float64[]
    
    # 動的閾値計算
    avg_power = mean(rx_power)
    dynamic_threshold = max(3 * avg_power, 10^(power_threshold_db / 10) * 1e-3)
    
    # 各信号タイミングでピーク検出
    for i in 1:signal_count
        expected_time = 10.0 + (i - 1) * interval_ms  # 開始遅延 + (i-1) * 間隔
        time_tolerance = interval_ms * 0.3  # 30%の許容誤差（広げて検出率向上）
        
        # 期待時刻周辺のサンプルを検索
        start_idx = findfirst(t -> t >= expected_time - time_tolerance, time_axis)
        end_idx = findlast(t -> t <= expected_time + time_tolerance, time_axis)
        
        if start_idx !== nothing && end_idx !== nothing
            search_range = start_idx:end_idx
            max_power_idx = argmax(rx_power[search_range])
            actual_idx = start_idx + max_power_idx - 1
            
            if rx_power[actual_idx] > dynamic_threshold
                # 上抜け点を探索（線形補間）
                crossing_time_ms = time_axis[actual_idx]
                for i in max(start_idx+1, actual_idx-100):actual_idx
                    p0 = rx_power[i-1]; p1 = rx_power[i]
                    if p0 <= dynamic_threshold && p1 > dynamic_threshold
                        t0 = time_axis[i-1]; t1 = time_axis[i]
                        frac = (dynamic_threshold - p0) / (p1 - p0 + 1e-20)
                        crossing_time_ms = t0 + frac * (t1 - t0)
                        break
                    end
                end
                push!(detected_peaks, actual_idx)
                push!(detection_times, crossing_time_ms)
                push!(peak_powers, rx_power[actual_idx])
            end
        end
    end
    
    detection_rate = length(detected_peaks) / signal_count * 100
    
    println("受信電力ピーク検出結果:")
    println("• 検出閾値: $(round(dynamic_threshold * 1000, digits=1)) mW ($(round(10*log10(dynamic_threshold*1000), digits=1)) dBm)")
    println("• 検出ピーク数: $(length(detected_peaks))/$(signal_count)")
    println("• 検出率: $(round(detection_rate, digits=1))%")
    
    return detected_peaks, detection_times, peak_powers, detection_rate
end

# ===== 同期ピーク検出（端末別閾値） =====
function detect_sync_peaks_adaptive(rx_power::Vector{Float64}, time_axis::Vector{Float64}, 
                                    signal_count::Int, interval_ms::Float64, terminal_threshold::Float64)
    detected_peaks = Int[]
    detection_times = Float64[]  # 閾値上抜け時刻（ms）
    peak_powers = Float64[]
    
    # 端末別閾値を使用
    threshold = terminal_threshold
    
    # 各信号タイミングでピーク検出
    for i in 1:signal_count
        expected_time = 10.0 + (i - 1) * interval_ms  # 開始遅延 + (i-1) * 間隔
        time_tolerance = interval_ms * 0.3  # 30%の許容誤差（広げて検出率向上）
        
        # 期待時刻周辺のサンプルを検索
        start_idx = findfirst(t -> t >= expected_time - time_tolerance, time_axis)
        end_idx = findlast(t -> t <= expected_time + time_tolerance, time_axis)
        
        if start_idx !== nothing && end_idx !== nothing
            search_range = start_idx:end_idx
            max_power_idx = argmax(rx_power[search_range])
            actual_idx = start_idx + max_power_idx - 1
            
            if rx_power[actual_idx] > threshold
                # 上抜け点を探索（線形補間）
                crossing_time_ms = time_axis[actual_idx]
                for i in max(start_idx+1, actual_idx-100):actual_idx
                    p0 = rx_power[i-1]; p1 = rx_power[i]
                    if p0 <= threshold && p1 > threshold
                        t0 = time_axis[i-1]; t1 = time_axis[i]
                        frac = (threshold - p0) / (p1 - p0 + 1e-20)
                        crossing_time_ms = t0 + frac * (t1 - t0)
                        break
                    end
                end
                push!(detected_peaks, actual_idx)
                push!(detection_times, crossing_time_ms)
                push!(peak_powers, rx_power[actual_idx])
            end
        end
    end
    
    detection_rate = length(detected_peaks) / signal_count * 100
    
    println("受信電力ピーク検出結果（適応的閾値）:")
    println("• 検出閾値: $(round(threshold * 1000, digits=1)) mW ($(round(10*log10(threshold*1000), digits=1)) dBm)")
    println("• 検出ピーク数: $(length(detected_peaks))/$(signal_count)")
    println("• 検出率: $(round(detection_rate, digits=1))%")
    
    return detected_peaks, detection_times, peak_powers, detection_rate
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
        4.7,                 # 中心周波数（GHz）
        3.6,                # 同期信号帯域幅（MHz）- 広帯域信号
        0.125,                # 端末受信帯域幅（MHz）- より現実的な帯域
        7.68,                # 送信側サンプリングレート（MHz）- 同期信号帯域幅の16倍
        0.25,                # 受信側サンプリングレート（MHz）- 端末受信帯域幅の2倍
        43.0,               # 送信電力（dBm）
        
        # 受信環境パラメータ
        true,               # シャドウイング有効/無効
        0.0,                 # シャドウイング標準偏差（dB）
        
        # 端末配置パラメータ
        "random_fixed",        # 配置モード（固定数ランダム配置）
        1,                  # 端末数（1端末）
        100.0,              # エリアサイズ（m）
        10.0,                # 最小距離（m）
        500.0,              # 最大距離（m）
        3.0,                # パスロス指数
        1.0,                # 参照距離（m）
        0.0,                 # 参照距離でのパスロス（dB）
        
        # シミュレーション制御パラメータ
        20.0,               # 信号間隔（ms）
        220.0,              # 総持続時間（ms）
        -100.0,             # ピーク検出閾値（dBm）
        
        # クロックパラメータ
        "none",         # クロック精度
        5.0,                # 温度変動（°C）
        100.0,              # 同期精度（ns）
        "median"           # 代表端末選択方法（"nearest"/"median"/"farthest"）
    )
end

# ===== メイン実行関数 =====
function main_simulation()
    Random.seed!(time_ns())  # 実行毎に異なるシードでランダム化（現在時刻ns）
    
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
    
    # 端末配置パラメータを作成
    deployment_params = TerminalDeploymentParameters(
        sim_params.deployment_mode,
        0.001,  # lambda (ポアソン点過程の密度) - 現実的な密度に調整
        sim_params.num_terminals,
        sim_params.area_size_m,
        sim_params.min_distance_m,
        sim_params.max_distance_m,
        sim_params.center_frequency_ghz * 1e9,
        sim_params.path_loss_exponent,
        sim_params.reference_distance_m,
        sim_params.reference_path_loss_db
    )
    
    # terminal_deployment.jlの関数を使用
    terminals = deploy_terminals(deployment_params, shadowing_params.std_db, shadowing_params.enabled, sim_params.tx_power_dbm)
    
    if isempty(terminals)
        println("警告: 端末が生成されませんでした。密度を上げるか、エリアサイズを大きくしてください。")
        return
    end
    
    # 端末配置統計分析
    analyze_terminal_deployment(terminals)
    println()
    
    println("端末数: $(length(terminals)) (シンプル化のため1端末のみ使用)")
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

    # ===== 7. 受信シミュレーション（1端末のみ） =====
    time_axis_rx = collect((0:length(tx_signal_low_rate)-1) / (sim_params.rx_sampling_rate_mhz * 1e6) * 1000)
    
    # 端末が1つなので、それを取得
    terminal = terminals[1]
    terminal_path_loss_params = PathLossParameters(
        terminal.distance_m,
        sim_params.center_frequency_ghz * 1e9,
        sim_params.path_loss_exponent,
        sim_params.reference_distance_m,
        sim_params.reference_path_loss_db
    )
    
    println("受信側：パスロス、シャドウイング、ノイズを付加中...")
    println("端末情報:")
    println("  - 位置: ($(round(terminal.x_m, digits=1)), $(round(terminal.y_m, digits=1))) m")
    println("  - 距離: $(round(terminal.distance_m, digits=1)) m")
    println("  - シャドウイング: $(round(terminal.shadowing_db, digits=2)) dB")
    
    rx_signal, actual_snr_db, path_loss_db = add_path_loss_and_noise(
        tx_signal_low_rate, terminal_path_loss_params, noise_params, shadowing_params, 
        sim_params.rx_sampling_rate_mhz * 1e6, terminal.shadowing_db
    )
    
    rx_power = abs2.(rx_signal)
    
    # ===== 8. ノイズフロア推定と閾値設定 =====
    println("\n=== ノイズフロア推定 ===")
    noise_floor_dbm = estimate_noise_floor_dbm(rx_power, time_axis_rx; pre_signal_end_ms=9.0)
    margin_db = 15.0
    threshold_dbm = noise_floor_dbm + margin_db
    threshold_w = 10^(threshold_dbm / 10) * 1e-3  # dBm -> W
    
    println("• ノイズフロア: $(round(noise_floor_dbm, digits=2)) dBm")
    println("• マージン: $(round(margin_db, digits=1)) dB")
    println("• 検出閾値: $(round(threshold_dbm, digits=2)) dBm ($(round(threshold_w * 1000, digits=6)) mW)")
    println()
    
    # ===== 9. シンプルな閾値超過検出（期待時刻不要） (★ここが変更点) =====
    println("=== 閾値超過検出 ===")
    println("受信電力が閾値を超えている全ての時間区間を検出中...")
    
    # 1. 生の閾値超過をすべて検出
    crossings_raw = detect_threshold_crossings(rx_power, time_axis_rx, threshold_w)
    
    num_crossings_raw = length(crossings_raw[:start_times])
    println("• 検出された閾値超過区間数（フィルタ前）: $num_crossings_raw")
    
    # ===== 2. サンプル数フィルタ =====
    sampling_interval_ms = 1.0 / (sim_params.rx_sampling_rate_mhz * 1e6) * 1000  # サンプリング間隔（ms）
    min_samples = 2  # 最低2サンプル以上にまたがる区間のみ採用
    
    filtered_indices = Int[]
    for i in 1:num_crossings_raw
        duration_ms = crossings_raw[:end_times][i] - crossings_raw[:start_times][i]
        num_samples = ceil(duration_ms / sampling_interval_ms)
        
        if num_samples >= min_samples
            push!(filtered_indices, i)
        end
    end
    
    # フィルタリング後の結果を作成 (crossings_filtered_samples)
    local crossings_filtered_samples::Dict{Symbol, Vector{Float64}}
    if !isempty(filtered_indices)
        crossings_filtered_samples = Dict{Symbol, Vector{Float64}}(
            :start_times => crossings_raw[:start_times][filtered_indices],
            :end_times => crossings_raw[:end_times][filtered_indices],
            :peak_times => crossings_raw[:peak_times][filtered_indices],
            :peak_powers_w => crossings_raw[:peak_powers_w][filtered_indices],
            :peak_powers_dbm => crossings_raw[:peak_powers_dbm][filtered_indices]
        )
    else
        crossings_filtered_samples = Dict{Symbol, Vector{Float64}}(
            :start_times => Float64[],
            :end_times => Float64[],
            :peak_times => Float64[],
            :peak_powers_w => Float64[],
            :peak_powers_dbm => Float64[]
        )
    end
    
    num_crossings_samples = length(crossings_filtered_samples[:start_times])
    println("• フィルタリング後（最低$(min_samples)サンプル以上）: $num_crossings_samples")

    # ===== 3. デバウンス処理フィルタ（ここが新しい） =====
    # 信号持続時間(約0.067ms)より長く、信号間隔(20ms)より十分短い値を設定
    debounce_time_ms = 1.0  # 1.0ms 以内に次の有効なピークが来たら無視する
    
    # 最終的な検出結果
    crossings = filter_by_debounce(crossings_filtered_samples, debounce_time_ms)
    
    num_crossings = length(crossings[:start_times])
    println("• フィルタリング後（デバウンス処理 $(debounce_time_ms) ms）: $num_crossings")
    
    
    if num_crossings > 0
        println("\n=== 検出された同期信号区間（フィルタリング後） ===")
        for i in 1:num_crossings
            start_time = crossings[:start_times][i]
            end_time = crossings[:end_times][i]
            peak_time = crossings[:peak_times][i]
            peak_power_dbm = crossings[:peak_powers_dbm][i]
            duration = end_time - start_time
            
            num_samples = ceil(duration / sampling_interval_ms)
            println("区間 $i:")
            println("  - 開始時刻: $(round(start_time, digits=3)) ms (閾値上抜け)")
            println("  - 終了時刻: $(round(end_time, digits=3)) ms (閾値下抜け)")
            println("  - 持続時間: $(round(duration, digits=3)) ms (約$(Int(num_samples))サンプル)")
            println("  - ピーク時刻: $(round(peak_time, digits=3)) ms")
            println("  - ピーク電力: $(round(peak_power_dbm, digits=2)) dBm")
            println()
        end
    else
        println("警告: 閾値を超える区間が検出されませんでした")
        println("      閾値が高すぎる可能性があります")
    end
    
    # ===== 10. 結果保存（シンプル版） =====
    # println("信号可視化を実行中...")
    # visualize_signal_analysis(tx_signal_high_rate, time_axis_tx, 
    #                          tx_signal_low_rate, time_axis_rx,
    #                          rx_power, peak_times, peak_powers,
    #                          "results_modular_simulation")
    
    output_dir = "results_modular_simulation"
    mkpath(output_dir)
    
    execution_timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    println("\n=== CSVファイル保存 ===")
    
    # 1. 受信電力データ（時系列）
    println("受信電力データを保存中...")
    min_length = min(length(time_axis_rx), length(rx_power))
    power_data = DataFrame(
        time_ms = time_axis_rx[1:min_length],
        power_mw = rx_power[1:min_length] * 1000,
        power_dbm = 10*log10.(rx_power[1:min_length] * 1000)
    )
    CSV.write("$(output_dir)/received_power_data_$(execution_timestamp).csv", power_data)
    
    # 2. 検出された閾値超過区間の情報
    if num_crossings > 0
        println("検出された閾値超過区間の情報を保存中...")
        durations_ms = crossings[:end_times] .- crossings[:start_times]
        num_samples_vec = [Int(ceil(d / sampling_interval_ms)) for d in durations_ms]
        
        crossings_data = DataFrame(
            interval_id = 1:num_crossings,
            start_time_ms = crossings[:start_times],
            end_time_ms = crossings[:end_times],
            duration_ms = durations_ms,
            num_samples = num_samples_vec,
            peak_time_ms = crossings[:peak_times],
            peak_power_mw = crossings[:peak_powers_w] .* 1000,
            peak_power_dbm = crossings[:peak_powers_dbm]
        )
        CSV.write("$(output_dir)/threshold_crossings_$(execution_timestamp).csv", crossings_data)
        println("  - 検出区間数: $num_crossings")
    end
    
    # 3. 端末情報と閾値設定
    println("端末情報と閾値設定を保存中...")
    terminal_info = DataFrame(
        terminal_id = [1],
        x_m = [terminal.x_m],
        y_m = [terminal.y_m],
        distance_m = [terminal.distance_m],
        path_loss_db = [terminal.path_loss_db],
        shadowing_db = [terminal.shadowing_db],
        rx_power_dbm = [terminal.rx_power_dbm],
        noise_floor_dbm = [noise_floor_dbm],
        threshold_margin_db = [margin_db],
        threshold_dbm = [threshold_dbm]
    )
    CSV.write("$(output_dir)/terminal_info_$(execution_timestamp).csv", terminal_info)
    
    println("\nCSVファイルを '$(output_dir)' に保存しました:")
    println("  • received_power_data_$(execution_timestamp).csv - 受信電力時系列データ")
    if num_crossings > 0
        println("  • threshold_crossings_$(execution_timestamp).csv - 検出された閾値超過区間")
    end
    println("  • terminal_info_$(execution_timestamp).csv - 端末情報と閾値設定")
    println()
    
    println("="^60)
    println("プログラム完了")
    println("="^60)
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    main_simulation()
end