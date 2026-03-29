# ===== 同期信号観測機能 =====

using Random, Statistics, LinearAlgebra

# ===== 同期信号観測記録 =====
struct SyncObservation
    terminal_id::Int                    # 端末ID
    global_detection_time_s::Float64    # グローバル時刻での検出時刻（s）
    local_detection_time_s::Float64     # ローカル時刻での検出時刻（s）
    clock_offset_ns::Float64           # クロックオフセット（ns）
    frequency_offset_ppm::Float64      # 周波数オフセット（ppm）
    signal_power_dbm::Float64          # 受信電力（dBm）
    detection_confidence::Float64      # 検出信頼度（0-1）
    temperature_c::Float64             # 温度（°C）
    observation_id::Int                # 観測ID
end

# ===== 同期信号観測器 =====
mutable struct SyncObserver
    terminal_id::Int                    # 端末ID
    local_clock::LocalClock            # ローカルクロック
    observations::Vector{SyncObservation}  # 観測記録
    next_observation_id::Int           # 次の観測ID
    detection_threshold_dbm::Float64   # 検出閾値（dBm）
    confidence_threshold::Float64      # 信頼度閾値
end

# ===== 同期信号観測器初期化 =====
function initialize_sync_observer(terminal_id::Int, local_clock::LocalClock, 
                                 detection_threshold_dbm::Float64, confidence_threshold::Float64 = 0.7)
    return SyncObserver(
        terminal_id,
        local_clock,
        SyncObservation[],
        1,
        detection_threshold_dbm,
        confidence_threshold
    )
end

# ===== 同期信号観測 =====
function observe_sync_signal(observer::SyncObserver, global_detection_time_s::Float64,
                            signal_power_dbm::Float64, detection_confidence::Float64 = 1.0)
    
    # ローカル時刻を計算
    local_detection_time_s = convert_to_local_time(observer.local_clock, global_detection_time_s)
    
    # クロックオフセットを計算
    clock_offset_ns = (local_detection_time_s - global_detection_time_s) * 1e9
    
    # 現在の周波数オフセットを取得
    current_frequency_offset_ppm = observer.local_clock.frequency_offset_ppm
    
    # 現在の温度を取得
    current_temperature_c = observer.local_clock.temperature_c
    
    # 観測記録を作成
    observation = SyncObservation(
        observer.terminal_id,
        global_detection_time_s,
        local_detection_time_s,
        clock_offset_ns,
        current_frequency_offset_ppm,
        signal_power_dbm,
        detection_confidence,
        current_temperature_c,
        observer.next_observation_id
    )
    
    # 観測記録を追加
    push!(observer.observations, observation)
    observer.next_observation_id += 1
    
    return observation
end

# ===== 複数同期信号の一括観測 =====
function observe_multiple_sync_signals(observer::SyncObserver, global_detection_times_s::Vector{Float64},
                                      signal_powers_dbm::Vector{Float64}, 
                                      detection_confidences::Vector{Float64} = Float64[])
    
    observations = SyncObservation[]
    
    # 信頼度が指定されていない場合は1.0で初期化
    if isempty(detection_confidences)
        detection_confidences = ones(Float64, length(global_detection_times_s))
    end
    
    for (i, global_time) in enumerate(global_detection_times_s)
        power = i <= length(signal_powers_dbm) ? signal_powers_dbm[i] : -100.0
        confidence = i <= length(detection_confidences) ? detection_confidences[i] : 1.0
        
        observation = observe_sync_signal(observer, global_time, power, confidence)
        push!(observations, observation)
    end
    
    return observations
end

# ===== 観測記録の分析 =====
function analyze_sync_observations(observer::SyncObserver)
    if isempty(observer.observations)
        println("端末 $(observer.terminal_id): 観測記録がありません")
        return Dict()
    end
    
    # 基本統計を計算
    global_times = [obs.global_detection_time_s for obs in observer.observations]
    local_times = [obs.local_detection_time_s for obs in observer.observations]
    clock_offsets = [obs.clock_offset_ns for obs in observer.observations]
    frequency_offsets = [obs.frequency_offset_ppm for obs in observer.observations]
    signal_powers = [obs.signal_power_dbm for obs in observer.observations]
    confidences = [obs.detection_confidence for obs in observer.observations]
    
    # 統計計算
    stats = Dict(
        "observation_count" => length(observer.observations),
        "clock_offset_mean_ns" => mean(clock_offsets),
        "clock_offset_std_ns" => std(clock_offsets),
        "clock_offset_range_ns" => maximum(clock_offsets) - minimum(clock_offsets),
        "frequency_offset_mean_ppm" => mean(frequency_offsets),
        "frequency_offset_std_ppm" => std(frequency_offsets),
        "signal_power_mean_dbm" => mean(signal_powers),
        "signal_power_std_dbm" => std(signal_powers),
        "confidence_mean" => mean(confidences),
        "confidence_std" => std(confidences),
        "temperature_mean_c" => mean([obs.temperature_c for obs in observer.observations]),
        "temperature_std_c" => std([obs.temperature_c for obs in observer.observations])
    )
    
    return stats
end

# ===== 観測記録の表示 =====
function display_sync_observations(observer::SyncObserver, max_observations::Int = 10)
    println("=== 端末 $(observer.terminal_id) の同期信号観測記録 ===")
    
    if isempty(observer.observations)
        println("観測記録がありません")
        return
    end
    
    # 統計情報を表示
    stats = analyze_sync_observations(observer)
    println("観測回数: $(stats["observation_count"])")
    println("クロックオフセット: $(round(stats["clock_offset_mean_ns"], digits=2)) ± $(round(stats["clock_offset_std_ns"], digits=2)) ns")
    println("周波数オフセット: $(round(stats["frequency_offset_mean_ppm"], digits=2)) ± $(round(stats["frequency_offset_std_ppm"], digits=2)) ppm")
    println("受信電力: $(round(stats["signal_power_mean_dbm"], digits=1)) ± $(round(stats["signal_power_std_dbm"], digits=1)) dBm")
    println("信頼度: $(round(stats["confidence_mean"], digits=3)) ± $(round(stats["confidence_std"], digits=3))")
    println("温度: $(round(stats["temperature_mean_c"], digits=1)) ± $(round(stats["temperature_std_c"], digits=1)) °C")
    println()
    
    # 個別観測記録を表示（最大10件）
    display_count = min(max_observations, length(observer.observations))
    println("個別観測記録（最新$(display_count)件）:")
    println("ID | グローバル時刻(s) | ローカル時刻(s) | オフセット(ns) | 電力(dBm) | 信頼度")
    println("-" ^ 80)
    
    for i in (length(observer.observations)-display_count+1):length(observer.observations)
        obs = observer.observations[i]
        println("$(obs.observation_id) | $(round(obs.global_detection_time_s, digits=6)) | $(round(obs.local_detection_time_s, digits=6)) | $(round(obs.clock_offset_ns, digits=2)) | $(round(obs.signal_power_dbm, digits=1)) | $(round(obs.detection_confidence, digits=3))")
    end
    println()
end

# ===== クロック同期の実行 =====
function perform_clock_synchronization(observer::SyncObserver, reference_time_s::Float64,
                                      sync_accuracy_ns::Float64 = 100.0)
    
    if isempty(observer.observations)
        println("端末 $(observer.terminal_id): 同期に必要な観測記録がありません")
        return false
    end
    
    # 最新の観測記録を使用
    latest_observation = observer.observations[end]
    local_time_s = latest_observation.local_detection_time_s
    
    # クロック同期を実行
    synchronize_clock(observer.local_clock, reference_time_s, local_time_s, sync_accuracy_ns)
    
    println("端末 $(observer.terminal_id): クロック同期を実行")
    println("  - 基準時刻: $(round(reference_time_s, digits=6)) s")
    println("  - ローカル時刻: $(round(local_time_s, digits=6)) s")
    println("  - 同期精度: $(sync_accuracy_ns) ns")
    
    return true
end

# ===== スロット境界の予測 =====
function predict_slot_boundaries(observer::SyncObserver, reference_sync_time_s::Float64,
                                slot_duration_s::Float64, num_slots::Int)
    
    # ローカルクロックでスロット境界を計算
    local_slot_boundaries = calculate_slot_boundaries(
        observer.local_clock, reference_sync_time_s, slot_duration_s, num_slots
    )
    
    # グローバル時刻に変換
    global_slot_boundaries = Float64[]
    for local_time in local_slot_boundaries
        global_time = convert_to_global_time(observer.local_clock, local_time)
        push!(global_slot_boundaries, global_time)
    end
    
    return local_slot_boundaries, global_slot_boundaries
end

# ===== 観測記録のCSV保存（個別） =====
function save_observations_to_csv(observer::SyncObserver, output_dir::String, filename_prefix::String = "sync_observations")
    if isempty(observer.observations)
        println("端末 $(observer.terminal_id): 保存する観測記録がありません")
        return
    end
    
    # CSVデータを作成
    csv_data = DataFrame(
        observation_id = [obs.observation_id for obs in observer.observations],
        terminal_id = [obs.terminal_id for obs in observer.observations],
        global_detection_time_s = [obs.global_detection_time_s for obs in observer.observations],
        local_detection_time_s = [obs.local_detection_time_s for obs in observer.observations],
        clock_offset_ns = [obs.clock_offset_ns for obs in observer.observations],
        frequency_offset_ppm = [obs.frequency_offset_ppm for obs in observer.observations],
        signal_power_dbm = [obs.signal_power_dbm for obs in observer.observations],
        detection_confidence = [obs.detection_confidence for obs in observer.observations],
        temperature_c = [obs.temperature_c for obs in observer.observations]
    )
    
    # ファイル名を生成
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "$(filename_prefix)_terminal_$(observer.terminal_id)_$(timestamp).csv"
    filepath = joinpath(output_dir, filename)
    
    # CSVファイルに保存
    CSV.write(filepath, csv_data)
    
    println("端末 $(observer.terminal_id): 観測記録を $(filepath) に保存しました")
    
    return filepath
end

# ===== 全端末の観測記録を1つのCSVファイルに保存 =====
function save_all_observations_to_csv(observers::Vector{SyncObserver}, output_dir::String, filename_prefix::String = "sync_observations")
    # 全端末の観測記録を結合
    all_observations = SyncObservation[]
    
    for observer in observers
        append!(all_observations, observer.observations)
    end
    
    if isempty(all_observations)
        println("保存する観測記録がありません")
        return
    end
    
    # CSVデータを作成
    csv_data = DataFrame(
        observation_id = [obs.observation_id for obs in all_observations],
        terminal_id = [obs.terminal_id for obs in all_observations],
        global_detection_time_s = [obs.global_detection_time_s for obs in all_observations],
        local_detection_time_s = [obs.local_detection_time_s for obs in all_observations],
        clock_offset_ns = [obs.clock_offset_ns for obs in all_observations],
        frequency_offset_ppm = [obs.frequency_offset_ppm for obs in all_observations],
        signal_power_dbm = [obs.signal_power_dbm for obs in all_observations],
        detection_confidence = [obs.detection_confidence for obs in all_observations],
        temperature_c = [obs.temperature_c for obs in all_observations]
    )
    
    # 端末IDでソート
    sort!(csv_data, :terminal_id)
    
    # ファイル名を生成
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "$(filename_prefix)_all_terminals_$(timestamp).csv"
    filepath = joinpath(output_dir, filename)
    
    # CSVファイルに保存
    CSV.write(filepath, csv_data)
    
    println("全端末の観測記録を $(filepath) に保存しました")
    println("  - 総観測数: $(length(all_observations))")
    println("  - 端末数: $(length(observers))")
    
    return filepath
end

# ===== 観測時間を端末ごとの列に分けたCSVファイル保存 =====
function save_observation_times_for_plotting(observers::Vector{SyncObserver}, output_dir::String, filename_prefix::String = "sync_times")
    # 総観測数を確認
    total_observations = sum(length(observer.observations) for observer in observers)
    if total_observations == 0
        println("保存する観測記録がありません")
        return
    end

    # 縦持ち（tidy）形式の列を準備
    sample_index_vec = Int[]
    terminal_id_vec = Int[]
    reference_time_s_vec = Float64[]
    local_time_s_vec = Float64[]
    offset_s_vec = Float64[]
    offset_us_vec = Float64[]

    # 各端末の観測を縦持ちで格納（sample_indexは端末内で0始まり）
    for observer in observers
        for (idx, obs) in enumerate(observer.observations)
            push!(sample_index_vec, idx - 1)
            push!(terminal_id_vec, observer.terminal_id)
            push!(reference_time_s_vec, obs.global_detection_time_s)
            push!(local_time_s_vec, obs.local_detection_time_s)
            offset_s = obs.local_detection_time_s - obs.global_detection_time_s
            push!(offset_s_vec, offset_s)
            push!(offset_us_vec, offset_s * 1e6)
        end
    end

    # DataFrameを作成（縦持ち）
    times_data = DataFrame(
        sample_index = sample_index_vec,
        terminal_id = terminal_id_vec,
        reference_time_s = reference_time_s_vec,
        local_time_s = local_time_s_vec,
        offset_s = offset_s_vec,
        offset_us = offset_us_vec,
    )

    # ファイル名を生成
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "$(filename_prefix)_$(timestamp).csv"
    filepath = joinpath(output_dir, filename)

    # CSVファイルに保存
    CSV.write(filepath, times_data)

    println("観測時間データ（縦持ち・グラフ用）を $(filepath) に保存しました")
    println("  - 総観測数: $(nrow(times_data))")
    println("  - 端末数: $(length(observers))")

    return filepath
end

# ===== 端末パラメータを別ファイルに保存 =====
function save_terminal_parameters(observers::Vector{SyncObserver}, terminals::Vector{TerminalInfo}, output_dir::String, filename_prefix::String = "terminal_parameters")
    # 端末パラメータデータを作成
    terminal_params = DataFrame(
        terminal_id = Int[],
        x_m = Float64[],
        y_m = Float64[],
        distance_m = Float64[],
        clock_offset_ns = Float64[],
        frequency_offset_ppm = Float64[],
        temperature_c = Float64[],
        detection_threshold_dbm = Float64[],
        confidence_threshold = Float64[]
    )
    
    for observer in observers
        if !isempty(observer.observations)
            # 最新の観測記録からパラメータを取得
            latest_obs = observer.observations[end]
            term = terminals[observer.terminal_id]
            
            push!(terminal_params, (
                observer.terminal_id,
                term.x_m,
                term.y_m,
                term.distance_m,
                latest_obs.clock_offset_ns,
                latest_obs.frequency_offset_ppm,
                latest_obs.temperature_c,
                10 * log10(observer.detection_threshold_dbm * 1000),  # W → mW → dBm
                observer.confidence_threshold
            ))
        end
    end
    
    # ファイル名を生成
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "$(filename_prefix)_$(timestamp).csv"
    filepath = joinpath(output_dir, filename)
    
    # CSVファイルに保存
    CSV.write(filepath, terminal_params)
    
    println("端末パラメータを $(filepath) に保存しました")
    println("  - 端末数: $(nrow(terminal_params))")
    
    return filepath
end
