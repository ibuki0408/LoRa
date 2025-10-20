using Random, Statistics, Printf, FFTW, LinearAlgebra, Plots, DSP, Distributions

# ===== 同期ビーコン信号シミュレーション =====

# ビーコン信号の種類
@enum BeaconType GaussianPulse RectangularPulse ChirpPulse

# ビーコン信号のパラメータ
struct BeaconParameters
    signal_type::BeaconType
    duration_ms::Float64      # 信号持続時間（ms）
    peak_power_dBm::Float64   # ピーク電力（dBm）
    center_frequency::Float64 # 中心周波数（Hz）
    bandwidth::Float64         # 帯域幅（Hz）
end

# 端末の位置と受信特性
mutable struct Terminal
    id::Int
    x::Float64              # X座標（km）
    y::Float64              # Y座標（km）
    clock_drift_ppm::Float64 # クロックドリフト（ppm）
    local_time::Float64      # ローカル時刻
    detected_peaks::Vector{Tuple{Float64, Float64, Float64}} # (検出時刻, 電力, 信頼度)
    peak_detection_count::Int # ピーク検出回数
end

# ビーコン送信機
mutable struct BeaconTransmitter
    x::Float64              # X座標（km）
    y::Float64              # Y座標（km）
    parameters::BeaconParameters
    transmission_times::Vector{Float64} # 送信時刻履歴
end

# シミュレーション環境
mutable struct SimulationEnvironment
    terminals::Vector{Terminal}
    beacon::BeaconTransmitter
    path_loss_alpha::Float64    # パスロス係数α
    path_loss_beta::Float64     # パスロス係数β
    path_loss_gamma::Float64    # パスロス係数γ
    shadowing_std::Float64      # シャドウイング標準偏差
    noise_floor_dBm::Float64    # ノイズフロア
    sampling_rate::Float64      # サンプリングレート（Hz）
    enable_bandwidth_loss::Bool # 帯域幅損失を有効にするかどうか
end

# ===== ポアソン点過程による端末配置 =====

function generate_poisson_points(rng::AbstractRNG, radius::Float64, density::Float64, target_count::Int)
    """
    ポアソン点過程で端末位置を生成
    
    Args:
        rng: 乱数生成器
        radius: カバレッジ半径（km）
        density: 端末密度（台/km²）
        target_count: 目標端末数
    
    Returns:
        Vector{Tuple{Float64, Float64}}: 端末位置のリスト
    """
    positions = Tuple{Float64, Float64}[]
    
    # ポアソン点過程の実装
    # 円形エリア内でランダムに点を配置
    for i in 1:target_count
        # 極座標でランダムに位置を生成
        theta = rand(rng) * 2π
        r = radius * sqrt(rand(rng))  # 半径方向の一様分布
        
        # 直交座標に変換
        x = r * cos(theta)
        y = r * sin(theta)
        
        push!(positions, (x, y))
    end
    
    return positions
end

function generate_poisson_points_advanced(rng::AbstractRNG, radius::Float64, density::Float64, target_count::Int)
    """
    高度なポアソン点過程（実際のポアソン分布に基づく）
    """
    positions = Tuple{Float64, Float64}[]
    
    # エリア内の期待端末数を計算
    area = π * radius^2
    expected_count = density * area
    
    # ポアソン分布で実際の端末数を決定
    actual_count = rand(rng, Poisson(expected_count))
    
    # 目標端末数に調整
    actual_count = min(actual_count, target_count)
    
    for i in 1:actual_count
        # 極座標でランダムに位置を生成
        theta = rand(rng) * 2π
        r = radius * sqrt(rand(rng))
        
        x = r * cos(theta)
        y = r * sin(theta)
        
        push!(positions, (x, y))
    end
    
    return positions
end

# ===== ビーコン信号生成関数 =====

function generate_gaussian_pulse(parameters::BeaconParameters, duration_ms::Float64, 
                                sampling_rate::Float64, peak_power_dBm::Float64)
    num_samples = Int(ceil(duration_ms * sampling_rate / 1000.0))
    samples = Float64[]
    
    center_time = duration_ms / 2.0  # 信号の中心時刻
    sigma = duration_ms / 6.0      # 標準偏差（3σで99.7%をカバー）
    
    for i in 1:num_samples
        time_ms = (i - 1) * 1000.0 / sampling_rate
        # ガウシアンパルス
        amplitude = peak_power_dBm * exp(-((time_ms - center_time)^2) / (2 * sigma^2))
        push!(samples, amplitude)
    end
    
    return samples
end

function generate_rectangular_pulse(parameters::BeaconParameters, duration_ms::Float64, 
                                   sampling_rate::Float64, peak_power_dBm::Float64)
    num_samples = Int(ceil(duration_ms * sampling_rate / 1000.0))
    samples = Float64[]
    
    for i in 1:num_samples
        time_ms = (i - 1) * 1000.0 / sampling_rate
        # 矩形パルス
        if time_ms <= duration_ms
            amplitude = peak_power_dBm
        else
            amplitude = -200.0  # 信号外は非常に低い電力
        end
        push!(samples, amplitude)
    end
    
    return samples
end

function generate_chirp_pulse(parameters::BeaconParameters, duration_ms::Float64, 
                             sampling_rate::Float64, peak_power_dBm::Float64)
    num_samples = Int(ceil(duration_ms * sampling_rate / 1000.0))
    samples = Float64[]
    
    for i in 1:num_samples
        time_ms = (i - 1) * 1000.0 / sampling_rate
        # チャープ信号（周波数が時間とともに変化）
        if time_ms <= duration_ms
            # 線形チャープ
            freq_slope = parameters.bandwidth / (duration_ms / 1000.0)
            instantaneous_freq = parameters.center_frequency + freq_slope * (time_ms / 1000.0)
            phase = 2π * instantaneous_freq * (time_ms / 1000.0)
            amplitude = peak_power_dBm * sin(phase)
        else
            amplitude = -200.0
        end
        push!(samples, amplitude)
    end
    
    return samples
end

# ビーコン信号生成の統合関数
function generate_beacon_signal(parameters::BeaconParameters, duration_ms::Float64, 
                               sampling_rate::Float64, peak_power_dBm::Float64)
    if parameters.signal_type == GaussianPulse
        return generate_gaussian_pulse(parameters, duration_ms, sampling_rate, peak_power_dBm)
    elseif parameters.signal_type == RectangularPulse
        return generate_rectangular_pulse(parameters, duration_ms, sampling_rate, peak_power_dBm)
    elseif parameters.signal_type == ChirpPulse
        return generate_chirp_pulse(parameters, duration_ms, sampling_rate, peak_power_dBm)
    else
        error("Unknown beacon signal type")
    end
end

# ===== 伝搬モデル =====

function calculate_path_loss(distance_km::Float64, alpha::Float64, beta::Float64, gamma::Float64, 
                           frequency_hz::Float64)
    return 10 * alpha * log10(distance_km) + beta + 10 * gamma * log10(frequency_hz / 1e6)
end

function calculate_position_based_shadowing(x::Float64, y::Float64, shadowing_std::Float64)
    """
    端末の位置に基づいて決定論的なシャドウイング値を計算
    
    Args:
        x, y: 端末の座標（km）
        shadowing_std: シャドウイング標準偏差（dB）
    
    Returns:
        Float64: 位置に基づくシャドウイング値（dB）
    """
    # 位置をハッシュ化してシードとして使用
    # 座標を整数に変換してハッシュ値を生成
    x_int = Int(round(x * 1000))  # メートル単位に変換して整数化
    y_int = Int(round(y * 1000))
    
    # ハッシュ値からシードを生成
    seed = abs(hash((x_int, y_int))) % 2^31
    
    # シードを設定してランダム値を生成
    rng = MersenneTwister(seed)
    shadowing_value = randn(rng) * shadowing_std
    
    return shadowing_value
end

function calculate_received_power(tx_power_dBm::Float64, distance_km::Float64, 
                                path_loss_dB::Float64, shadowing_dB::Float64, noise_dB::Float64)
    return tx_power_dBm - path_loss_dB - shadowing_dB + noise_dB
end

# ===== ピーク検出アルゴリズム =====

function detect_peak_in_signal(signal_power::Vector{Float64}, sampling_rate::Float64, 
                              noise_floor_dBm::Float64, threshold_dB::Float64)
    if isempty(signal_power)
        return nothing, 0.0
    end
    
    # 最大電力とそのインデックスを検索
    max_power, max_index = findmax(signal_power)
    
    # 閾値チェック
    if max_power < noise_floor_dBm + threshold_dB
        return nothing, 0.0
    end
    
    # ピーク時刻を計算
    peak_time = (max_index - 1) / sampling_rate
    
    return peak_time, max_power
end

# 高度なピーク検出（複数候補から最適なものを選択）
function detect_peak_advanced(signal_power::Vector{Float64}, sampling_rate::Float64, 
                             noise_floor_dBm::Float64, threshold_dB::Float64)
    if isempty(signal_power)
        return nothing, 0.0, 0.0
    end
    
    # 閾値以上の電力を持つサンプルを検索
    valid_indices = findall(p -> p >= noise_floor_dBm + threshold_dB, signal_power)
    
    if isempty(valid_indices)
        return nothing, 0.0, 0.0
    end
    
    # 最大電力のサンプルを選択
    max_power = maximum(signal_power[valid_indices])
    max_index = valid_indices[argmax(signal_power[valid_indices])]
    
    # ピーク時刻を計算
    peak_time = (max_index - 1) / sampling_rate
    
    # 信頼度を計算（ノイズフロアからの距離に基づく）
    confidence = min(1.0, (max_power - noise_floor_dBm) / (threshold_dB * 2))
    
    return peak_time, max_power, confidence
end

# ===== ブロードキャスト同期信号のピーク検出シミュレーション =====

function simulate_broadcast_peak_detection(env::SimulationEnvironment, transmission_time::Float64, 
                                          rng::AbstractRNG)
    detection_results = Tuple{Int, Float64, Float64, Float64}[]  # (端末ID, 検出時刻, 電力, 信頼度)
    
    println("  ブロードキャスト送信時刻: $(round(transmission_time, digits=6))s")
    println("  送信機位置: ($(env.beacon.x)km, $(env.beacon.y)km)")
    
    # ブロードキャスト信号の生成（全端末共通の基準信号）
    # 実際のブロードキャストでは、送信機から全方向に同じ信号が送信される
    base_beacon_signal = generate_beacon_signal(
        env.beacon.parameters, 
        env.beacon.parameters.duration_ms, 
        env.sampling_rate, 
        env.beacon.parameters.peak_power_dBm  # 基準電力
    )
    
    println("  ブロードキャスト信号: $(env.beacon.parameters.signal_type), 持続時間=$(env.beacon.parameters.duration_ms)ms")
    println("  ブロードキャスト帯域幅: $(env.beacon.parameters.bandwidth/1e3)kHz")
    println("  帯域幅損失: $(env.enable_bandwidth_loss ? "有効" : "無効")")
    
    if env.enable_bandwidth_loss
        println("  端末受信帯域幅: 125kHz (狭帯域)")
        println("  帯域幅比: $(round(125e3/env.beacon.parameters.bandwidth, digits=3))")
        println("  帯域幅損失: $(round(10*log10(125e3/env.beacon.parameters.bandwidth), digits=1))dB")
    else
        println("  端末受信帯域幅: $(env.beacon.parameters.bandwidth/1e3)kHz (広帯域)")
        println("  帯域幅損失: 0dB")
    end
    
    # 各端末での受信処理（ブロードキャスト信号を受信）
    for terminal in env.terminals
        # 距離計算
        distance_km = sqrt((terminal.x - env.beacon.x)^2 + (terminal.y - env.beacon.y)^2)
        
        # 伝搬損失計算
        path_loss_dB = calculate_path_loss(distance_km, env.path_loss_alpha, 
                                         env.path_loss_beta, env.path_loss_gamma, 
                                         env.beacon.parameters.center_frequency)
        
        # シャドウイング（端末の位置に基づいて固定値）
        shadowing_dB = calculate_position_based_shadowing(terminal.x, terminal.y, env.shadowing_std)
        
        # 受信電力計算（ブロードキャスト信号の受信電力）
        received_power_dBm = calculate_received_power(
            env.beacon.parameters.peak_power_dBm, distance_km, 
            path_loss_dB, shadowing_dB, 0.0
        )
        
        # 帯域幅損失の計算（有効な場合のみ）
        if env.enable_bandwidth_loss
            # 帯域幅損失ありバージョン（狭帯域端末での受信）
            # 端末の受信帯域幅（例：125kHz）とブロードキャスト信号の帯域幅（例：1MHz）の比
            terminal_bandwidth = 125e3  # 端末の受信帯域幅（125kHz）
            broadcast_bandwidth = env.beacon.parameters.bandwidth  # ブロードキャスト信号の帯域幅
            
            # 帯域幅比による電力損失
            bandwidth_ratio = terminal_bandwidth / broadcast_bandwidth
            bandwidth_loss_dB = 10 * log10(bandwidth_ratio)  # 帯域幅損失
            
            # 実際の受信電力（広帯域信号の一部のみ）
            actual_received_power_dBm = received_power_dBm + bandwidth_loss_dB
        else
            # 帯域幅損失なしバージョン（広帯域端末での受信）
            actual_received_power_dBm = received_power_dBm
        end
        
        # 受信信号の生成（帯域幅損失の有無に応じて）
        power_scaling = 10^((actual_received_power_dBm - env.beacon.parameters.peak_power_dBm) / 20.0)
        terminal_beacon_signal = base_beacon_signal .* power_scaling
        
        # 各端末固有のノイズを追加
        noise_signal = env.noise_floor_dBm .+ randn(rng, length(terminal_beacon_signal)) * 2.0
        total_signal = terminal_beacon_signal .+ noise_signal
        
        # ピーク検出
        peak_time, peak_power, confidence = detect_peak_advanced(
            total_signal, env.sampling_rate, 
            env.noise_floor_dBm, 3.0  # 3dB閾値
        )
        
        if peak_time !== nothing
            # 検出結果を記録
            detection_time = transmission_time + peak_time
            push!(detection_results, (terminal.id, detection_time, peak_power, confidence))
            
            # 端末の履歴を更新
            push!(terminal.detected_peaks, (detection_time, peak_power, confidence))
            terminal.peak_detection_count += 1
            
            println("    端末$(terminal.id): ピーク検出! 時刻=$(round(detection_time, digits=6))s, 電力=$(round(peak_power, digits=1))dBm, 信頼度=$(round(confidence, digits=2)) (距離=$(round(distance_km*1000, digits=1))m)")
        else
            println("    端末$(terminal.id): ピーク未検出 (距離=$(round(distance_km*1000, digits=1))m, 受信電力=$(round(received_power_dBm, digits=1))dBm)")
        end
    end
    
    return detection_results
end

# ===== シミュレーション実行 =====

function run_sync_simulation(num_terminals::Int=10, 
                            beacon_interval::Float64=0.04, simulation_duration::Float64=1.0,
                            enable_bandwidth_loss::Bool=true, duration_ms::Float64=5.0,
                            rng::AbstractRNG=Random.default_rng())
    
    println("=== 同期ビーコン信号シミュレーション開始 ===")
    
    # 環境設定（広帯域ブロードキャスト信号）
    env = SimulationEnvironment(
        Terminal[],  # 端末リスト
        BeaconTransmitter(0.0, 0.0, BeaconParameters(ChirpPulse, duration_ms, 20.0, 923.2e6, 1e6), Float64[]),  # 1MHz広帯域
        4.0, 9.5, 4.5,  # パスロス係数
        3.48,  # シャドウイング標準偏差
        -100.0,  # ノイズフロア
        10000.0,  # サンプリングレート（10kHz）
        enable_bandwidth_loss  # 帯域幅損失の有無
    )
    
    # ポアソン点過程で端末を配置
    println("ポアソン点過程で端末を配置中...")
    
    # ポアソン点過程のパラメータ
    coverage_radius = 1.0  # 1km半径のカバレッジエリア
    terminal_density = num_terminals / (π * coverage_radius^2)  # 端末密度（台/km²）
    
    # ポアソン点過程で端末位置を生成
    terminal_positions = generate_poisson_points(rng, coverage_radius, terminal_density, num_terminals)
    
    for (i, (x, y)) in enumerate(terminal_positions)
        terminal = Terminal(
            i, x, y, 
            rand(rng) * 100.0 - 50.0,  # ±50ppmドリフト
            0.0, Tuple{Float64, Float64, Float64}[], 0
        )
        push!(env.terminals, terminal)
    end
    
    println("端末配置完了: $(length(env.terminals))台")
    println("カバレッジ半径: $(coverage_radius)km")
    println("端末密度: $(round(terminal_density, digits=2))台/km²")
    
    # ビーコン送信時刻を生成
    beacon_times = collect(0.0:beacon_interval:simulation_duration)
    
    all_detection_results = Tuple{Int, Float64, Float64, Float64}[]
    
    println("端末数: $num_terminals")
    println("ブロードキャスト送信回数: $(length(beacon_times))")
    println("送信間隔: $(beacon_interval)s")
    println("信号タイプ: $(env.beacon.parameters.signal_type)")
    println("信号持続時間: $(env.beacon.parameters.duration_ms)ms")
    println("送信電力: $(env.beacon.parameters.peak_power_dBm)dBm")
    println("送信機位置: ($(env.beacon.x)km, $(env.beacon.y)km)")
    println("帯域幅損失: $(env.enable_bandwidth_loss ? "有効" : "無効")")
    println()
    
    # 各ブロードキャスト送信でのシミュレーション
    for (i, beacon_time) in enumerate(beacon_times)
        println("=== ブロードキャスト送信 $i/$(length(beacon_times)) ===")
        
        # ブロードキャストピーク検出シミュレーション
        detection_results = simulate_broadcast_peak_detection(env, beacon_time, rng)
        
        # 結果を蓄積
        append!(all_detection_results, detection_results)
        
        # ビーコン送信履歴を更新
        push!(env.beacon.transmission_times, beacon_time)
        println()
    end
    
    return env, all_detection_results
end

# ===== ファイル名生成関数 =====

function generate_filename_prefix(num_terminals::Int, beacon_interval::Float64, simulation_duration::Float64, 
                                signal_type::BeaconType, peak_power_dBm::Float64, duration_ms::Float64)
    """
    シミュレーションパラメータからファイル名プレフィックスを生成
    """
    # 信号タイプ
    signal_suffix = string(signal_type)
    
    # パラメータを組み合わせてファイル名プレフィックスを生成
    prefix = "sync_$(num_terminals)term_$(Int(beacon_interval*1000))ms_$(Int(simulation_duration))s_narrow_$(signal_suffix)_$(Int(peak_power_dBm))dBm_$(Int(duration_ms))ms"
    
    return prefix
end

# ===== 結果分析と可視化 =====

function analyze_peak_detection_results(env::SimulationEnvironment, all_detection_results::Vector{Tuple{Int, Float64, Float64, Float64}})
    
    println("\n=== ピーク検出結果分析 ===")
    
    if isempty(all_detection_results)
        println("検出されたピークがありません。")
        return
    end
    
    # 基本統計
    terminal_ids = [result[1] for result in all_detection_results]
    detection_times = [result[2] for result in all_detection_results]
    powers = [result[3] for result in all_detection_results]
    confidences = [result[4] for result in all_detection_results]
    
    println("ピーク検出統計:")
    println("  総検出数: $(length(all_detection_results))")
    println("  検出端末数: $(length(unique(terminal_ids)))")
    
    # 電力統計
    mean_power = mean(powers)
    std_power = std(powers)
    println("\n受信電力統計:")
    println("  平均電力: $(round(mean_power, digits=1)) dBm")
    println("  標準偏差: $(round(std_power, digits=1)) dBm")
    println("  最大電力: $(round(maximum(powers), digits=1)) dBm")
    println("  最小電力: $(round(minimum(powers), digits=1)) dBm")
    
    # 信頼度統計
    mean_confidence = mean(confidences)
    std_confidence = std(confidences)
    println("\n信頼度統計:")
    println("  平均信頼度: $(round(mean_confidence, digits=3))")
    println("  標準偏差: $(round(std_confidence, digits=3))")
    println("  最大信頼度: $(round(maximum(confidences), digits=3))")
    println("  最小信頼度: $(round(minimum(confidences), digits=3))")
    
    # 端末別統計
    println("\n端末別ピーク検出統計:")
    for terminal in env.terminals
        terminal_detections = [result for result in all_detection_results if result[1] == terminal.id]
        if !isempty(terminal_detections)
            terminal_powers = [det[3] for det in terminal_detections]
            terminal_confidences = [det[4] for det in terminal_detections]
            mean_terminal_power = mean(terminal_powers)
            mean_terminal_confidence = mean(terminal_confidences)
            println("  端末$(terminal.id): 検出数=$(length(terminal_detections)), 平均電力=$(round(mean_terminal_power, digits=1))dBm, 平均信頼度=$(round(mean_terminal_confidence, digits=3))")
        else
            println("  端末$(terminal.id): 検出数=0")
        end
    end
end

function export_power_data_to_csv(env::SimulationEnvironment, all_detection_results::Vector{Tuple{Int, Float64, Float64, Float64}}, 
                                  filename_prefix::String="")
    
    if isempty(all_detection_results)
        println("エクスポートするデータがありません。")
        return
    end
    
    # 各端末ごとにデータを分離
    terminal_data = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    
    for result in all_detection_results
        terminal_id = result[1]
        detection_time = result[2]
        power = result[3]
        
        if !haskey(terminal_data, terminal_id)
            terminal_data[terminal_id] = Tuple{Float64, Float64}[]
        end
        push!(terminal_data[terminal_id], (detection_time, power))
    end
    
    # 各端末ごとにCSVファイルを生成
    for terminal in env.terminals
        terminal_id = terminal.id
        
        if !haskey(terminal_data, terminal_id)
            println("端末$(terminal_id): 検出データがありません。")
            continue
        end
        
        # 連続的な電力レベルデータを生成
        time_range = 0.0:0.001:0.6  # 1ms間隔で0.6秒まで
        power_levels = Float64[]
        
        # ノイズフロアレベル
        noise_floor = -100.0
        
        for t in time_range
            # 基本ノイズレベル（より安定したノイズ）
            base_power = noise_floor + randn() * 0.3
            
            # この端末の検出時刻の周辺でピークを生成
            for (det_time, det_power) in terminal_data[terminal_id]
                time_diff = abs(t - det_time)
                if time_diff <= 0.015  # 15ms以内でピークを生成
                    # より鋭いガウシアン形状のピークを生成
                    sigma = 0.002  # 2msの標準偏差
                    peak_amplitude = det_power - noise_floor
                    peak_contribution = peak_amplitude * exp(-(time_diff^2) / (2 * sigma^2))
                    base_power += peak_contribution
                end
            end
            
            push!(power_levels, base_power)
        end
        
        # 端末別CSVファイルに出力（パラメータを含むファイル名）
        csv_filename = "results_charp/$(filename_prefix)_terminal_$(terminal_id)_power_data.csv"
        open(csv_filename, "w") do io
            # ヘッダー行
            println(io, "time_seconds,power_dbm")
            
            # 連続的な電力データを出力
            for (i, (t, power)) in enumerate(zip(time_range, power_levels))
                println(io, "$(t),$(power)")
            end
        end
        
        println("端末$(terminal_id)の受信電力データを '$(csv_filename)' に保存しました。")
    end
    
    println("データ形式: 時間(秒), 電力(dBm)")
    
    return length(terminal_data)
end

# 端末配置の可視化
function plot_terminal_layout(env::SimulationEnvironment, filename_prefix::String="")
    """
    端末配置の可視化
    """
    # 端末位置を取得
    terminal_x = [t.x for t in env.terminals]
    terminal_y = [t.y for t in env.terminals]
    
    # ゲートウェイ位置
    gw_x = env.beacon.x
    gw_y = env.beacon.y
    
    # プロット作成
    p = plot(terminal_x, terminal_y,
             xlabel="X Position (km)",
             ylabel="Y Position (km)",
             title="Terminal Layout (Poisson Point Process)",
             marker=:circle,
             markersize=6,
             color=:blue,
             label="Terminals",
             aspect_ratio=:equal)
    
    # ゲートウェイをマーク
    scatter!([gw_x], [gw_y],
             color=:red,
             markersize=10,
             marker=:star,
             label="Gateway (Beacon)")
    
    # カバレッジエリアを表示
    coverage_radius = 1.0
    theta_range = 0:0.1:2π
    coverage_x = coverage_radius .* cos.(theta_range)
    coverage_y = coverage_radius .* sin.(theta_range)
    
    plot!(coverage_x, coverage_y,
          color=:red,
          linestyle=:dash,
          linewidth=2,
          label="Coverage Area")
    
    # 各端末からゲートウェイまでの距離を表示
    for (i, terminal) in enumerate(env.terminals)
        distance = sqrt((terminal.x - gw_x)^2 + (terminal.y - gw_y)^2)
        annotate!(terminal.x, terminal.y, text("$(round(distance*1000, digits=0))m", 8, :black))
    end
    
    # パラメータを含むファイル名を生成
    layout_filename = isempty(filename_prefix) ? "results_charp/terminal_layout.png" : "results_charp/$(filename_prefix)_terminal_layout.png"
    savefig(p, layout_filename)
    println("端末配置プロットを '$(layout_filename)' に保存しました。")
    
    return p
end

# ===== メイン実行 =====

function main()
    println("ブロードキャスト同期信号ピーク検出シミュレーション")
    println("="^60)
    
    # シミュレーションパラメータ
    num_terminals = 1
    beacon_interval = 0.04
    simulation_duration = 1.0
    signal_type = ChirpPulse
    peak_power_dBm = 20.0
    bandwidth = 1e6
    duration_ms = 5.0  # 送信持続時間（ms）
    
    # シミュレーション実行
    println("\n=== シミュレーション実行 ===")
    env, all_detection_results = run_sync_simulation(
        1,  # num_terminals (1台の端末)
        0.04,  # beacon_interval (40ms間隔)
        0.6,   # simulation_duration (2.0秒間)
        true,  # enable_bandwidth_loss (帯域幅損失あり)
        duration_ms  # 送信持続時間
    )
    
    # ファイル名プレフィックスを生成
    filename_prefix = generate_filename_prefix(num_terminals, beacon_interval, simulation_duration, 
                                             signal_type, peak_power_dBm, duration_ms)
    
    # 端末配置の可視化
    plot_terminal_layout(env, filename_prefix)
    
    # 結果分析
    println("\n=== 結果分析 ===")
    analyze_peak_detection_results(env, all_detection_results)
    
    # CSV出力
    println("\n=== CSV出力 ===")
    export_power_data_to_csv(env, all_detection_results, filename_prefix)
    
    println("\nブロードキャストピーク検出シミュレーション完了！")
    
    return env, all_detection_results
end

# メイン実行
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
