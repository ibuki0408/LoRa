# ===== ローカルクロック機能 =====

using Random, Statistics, LinearAlgebra

# ===== クロックパラメータ =====
struct ClockParameters
    initial_frequency_ppm::Float64      # 初期周波数オフセット（ppm）
    frequency_drift_ppm::Float64        # 周波数ドリフト（ppm/時間）
    jitter_std_ns::Float64             # ジッター標準偏差（ns）
    temperature_coeff::Float64          # 温度係数（ppm/°C）
    aging_coeff::Float64               # エイジング係数（ppm/年）
    
    # コンストラクタ
    function ClockParameters(;
        initial_frequency_ppm::Float64,
        frequency_drift_ppm::Float64,
        jitter_std_ns::Float64,
        temperature_coeff::Float64,
        aging_coeff::Float64
    )
        new(initial_frequency_ppm, frequency_drift_ppm, jitter_std_ns, temperature_coeff, aging_coeff)
    end
end

# ===== ローカルクロック =====
mutable struct LocalClock
    terminal_id::Int                    # 端末ID
    reference_frequency_hz::Float64     # 基準周波数（Hz）
    current_frequency_hz::Float64       # 現在の周波数（Hz）
    frequency_offset_ppm::Float64       # 周波数オフセット（ppm）
    phase_offset_ns::Float64           # 位相オフセット（ns）
    accumulated_drift_ns::Float64       # 累積ドリフト（ns）
    last_update_time_s::Float64         # 最後の更新時刻（s）
    clock_parameters::ClockParameters   # クロックパラメータ
    temperature_c::Float64              # 温度（°C）
    age_years::Float64                 # 経過年数
end

# ===== デフォルトクロックパラメータ =====
function create_default_clock_parameters()
    return ClockParameters(
        initial_frequency_ppm = 20.0,      # 20ppm（一般的な水晶発振器）
        frequency_drift_ppm = 0.1,         # 0.1ppm/時間
        jitter_std_ns = 1.0,              # 1ns標準偏差
        temperature_coeff = 0.04,          # 0.04ppm/°C
        aging_coeff = 0.1                  # 0.1ppm/年
    )
end

# ===== ドリフトなしクロックパラメータ =====
function create_no_drift_clock_parameters()
    return ClockParameters(
        initial_frequency_ppm = 0.0,      # 初期オフセットなし
        frequency_drift_ppm = 0.0,        # ドリフトなし
        jitter_std_ns = 0.0,              # ジッターなし
        temperature_coeff = 0.0,          # 温度依存なし
        aging_coeff = 0.0                 # エイジングなし
    )
end

# ===== 高精度クロックパラメータ =====
function create_high_precision_clock_parameters()
    return ClockParameters(
        initial_frequency_ppm = 2.0,       # 2ppm（高精度水晶発振器）
        frequency_drift_ppm = 0.01,        # 0.01ppm/時間
        jitter_std_ns = 0.1,              # 0.1ns標準偏差
        temperature_coeff = 0.01,          # 0.01ppm/°C
        aging_coeff = 0.01                 # 0.01ppm/年
    )
end

# ===== 低精度クロックパラメータ =====
function create_low_precision_clock_parameters()
    return ClockParameters(
        initial_frequency_ppm = 100.0,     # 100ppm（低精度発振器）
        frequency_drift_ppm = 1.0,         # 1ppm/時間
        jitter_std_ns = 10.0,             # 10ns標準偏差
        temperature_coeff = 0.2,           # 0.2ppm/°C
        aging_coeff = 1.0                  # 1ppm/年
    )
end

# ===== ローカルクロック初期化 =====
function initialize_local_clock(terminal_id::Int, reference_frequency_hz::Float64, 
                               clock_params::ClockParameters, temperature_c::Float64 = 25.0)
    
    # 初期周波数オフセットを計算（ランダム + 温度依存）
    temp_offset_ppm = clock_params.temperature_coeff * (temperature_c - 25.0)
    # ドリフトなし設定ではランダム初期オフセットも付与しない
    random_initial_ppm = (clock_params.initial_frequency_ppm == 0.0 && clock_params.temperature_coeff == 0.0) ? 0.0 : randn() * 5.0
    initial_offset_ppm = clock_params.initial_frequency_ppm + temp_offset_ppm + random_initial_ppm
    
    # 初期周波数を計算
    initial_frequency_hz = reference_frequency_hz * (1.0 + initial_offset_ppm / 1e6)
    
    return LocalClock(
        terminal_id,
        reference_frequency_hz,
        initial_frequency_hz,
        initial_offset_ppm,
        0.0,  # 初期位相オフセット
        0.0,  # 初期累積ドリフト
        0.0,  # 初期更新時刻
        clock_params,
        temperature_c,
        0.0   # 初期経過年数
    )
end

# ===== クロック更新 =====
function update_local_clock(clock::LocalClock, current_time_s::Float64, 
                           temperature_c::Union{Float64, Nothing} = nothing)
    
    # 温度更新（指定されていない場合は現在の温度を維持）
    if temperature_c !== nothing
        clock.temperature_c = temperature_c
    end
    
    # 時間経過を計算
    time_elapsed_s = current_time_s - clock.last_update_time_s
    
    if time_elapsed_s <= 0.0
        return clock
    end
    
    # 周波数ドリフトを計算
    temp_drift_ppm = clock.clock_parameters.temperature_coeff * (clock.temperature_c - 25.0)
    aging_drift_ppm = clock.clock_parameters.aging_coeff * clock.age_years
    random_drift_ppm = randn() * clock.clock_parameters.frequency_drift_ppm * sqrt(time_elapsed_s)
    
    total_drift_ppm = temp_drift_ppm + aging_drift_ppm + random_drift_ppm
    
    # 周波数を更新
    clock.frequency_offset_ppm += total_drift_ppm
    clock.current_frequency_hz = clock.reference_frequency_hz * (1.0 + clock.frequency_offset_ppm / 1e6)
    
    # 位相ドリフトを計算
    phase_drift_ns = (clock.frequency_offset_ppm / 1e6) * time_elapsed_s * 1e9
    clock.accumulated_drift_ns += phase_drift_ns
    
    # ジッターを追加
    jitter_ns = randn() * clock.clock_parameters.jitter_std_ns
    clock.phase_offset_ns += jitter_ns
    
    # 時間を更新
    clock.last_update_time_s = current_time_s
    clock.age_years += time_elapsed_s / (365.25 * 24 * 3600)  # 年単位に変換
    
    return clock
end

# ===== 時刻変換 =====
function convert_to_local_time(clock::LocalClock, global_time_s::Float64)
    # クロックを更新
    update_local_clock(clock, global_time_s)
    
    # ローカル時刻を計算
    local_time_s = global_time_s + clock.accumulated_drift_ns / 1e9 + clock.phase_offset_ns / 1e9
    
    return local_time_s
end

function convert_to_global_time(clock::LocalClock, local_time_s::Float64)
    # クロックを更新
    update_local_clock(clock, local_time_s)
    
    # グローバル時刻を計算（逆変換）
    global_time_s = local_time_s - clock.accumulated_drift_ns / 1e9 - clock.phase_offset_ns / 1e9
    
    return global_time_s
end

# ===== 同期信号検出時刻の調整 =====
function adjust_detection_time(clock::LocalClock, detection_time_s::Float64, 
                              global_time_s::Float64)
    # クロックを更新
    update_local_clock(clock, global_time_s)
    
    # ローカル時刻での検出時刻を計算
    local_detection_time_s = convert_to_local_time(clock, detection_time_s)
    
    return local_detection_time_s
end

# ===== スロット境界計算 =====
function calculate_slot_boundaries(clock::LocalClock, sync_time_s::Float64, 
                                  slot_duration_s::Float64, num_slots::Int)
    # クロックを更新
    update_local_clock(clock, sync_time_s)
    
    slot_boundaries = Float64[]
    
    for i in 0:(num_slots-1)
        # グローバル時刻でのスロット境界
        global_slot_time_s = sync_time_s + i * slot_duration_s
        
        # ローカル時刻に変換
        local_slot_time_s = convert_to_local_time(clock, global_slot_time_s)
        push!(slot_boundaries, local_slot_time_s)
    end
    
    return slot_boundaries
end

# ===== クロック同期 =====
function synchronize_clock(clock::LocalClock, reference_time_s::Float64, 
                          local_time_s::Float64, sync_accuracy_ns::Float64 = 100.0)
    # 時刻差を計算
    time_diff_s = reference_time_s - local_time_s
    
    # 同期精度内の場合のみ調整
    if abs(time_diff_s * 1e9) <= sync_accuracy_ns
        # 位相オフセットを調整
        clock.phase_offset_ns += time_diff_s * 1e9
        
        # 周波数オフセットを部分的に調整（過度な調整を避ける）
        frequency_adjustment_ppm = time_diff_s * 1e6 / (reference_time_s - clock.last_update_time_s + 1e-6)
        clock.frequency_offset_ppm += frequency_adjustment_ppm * 0.1  # 10%の調整
        
        # 周波数を更新
        clock.current_frequency_hz = clock.reference_frequency_hz * (1.0 + clock.frequency_offset_ppm / 1e6)
    end
    
    return clock
end

# ===== クロック統計 =====
function analyze_clock_performance(clocks::Vector{LocalClock}, reference_time_s::Float64)
    println("=== ローカルクロック性能分析 ===")
    
    for clock in clocks
        # クロックを更新
        update_local_clock(clock, reference_time_s)
        
        # 統計を計算
        frequency_error_ppm = clock.frequency_offset_ppm
        phase_error_ns = clock.accumulated_drift_ns + clock.phase_offset_ns
        temperature_c = clock.temperature_c
        age_years = clock.age_years
        
        println("端末 $(clock.terminal_id):")
        println("  - 周波数オフセット: $(round(frequency_error_ppm, digits=2)) ppm")
        println("  - 位相エラー: $(round(phase_error_ns, digits=2)) ns")
        println("  - 温度: $(round(temperature_c, digits=1)) °C")
        println("  - 経過年数: $(round(age_years, digits=6)) 年")
        println()
    end
end

# ===== クロック品質評価 =====
function evaluate_clock_quality(clock::LocalClock, reference_time_s::Float64)
    # クロックを更新
    update_local_clock(clock, reference_time_s)
    
    # 品質指標を計算
    frequency_stability = abs(clock.frequency_offset_ppm) < 50.0  # 50ppm以下
    phase_stability = abs(clock.accumulated_drift_ns + clock.phase_offset_ns) < 1000.0  # 1μs以下
    temperature_stability = abs(clock.temperature_c - 25.0) < 10.0  # 25°C±10°C
    
    overall_quality = frequency_stability && phase_stability && temperature_stability
    
    return Dict(
        "frequency_stability" => frequency_stability,
        "phase_stability" => phase_stability,
        "temperature_stability" => temperature_stability,
        "overall_quality" => overall_quality,
        "frequency_error_ppm" => clock.frequency_offset_ppm,
        "phase_error_ns" => clock.accumulated_drift_ns + clock.phase_offset_ns
    )
end
