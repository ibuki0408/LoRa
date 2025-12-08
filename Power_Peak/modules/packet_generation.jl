module PacketGeneration

export generate_next_poisson_time

"""
    generate_next_poisson_time(current_time_ms, mean_interval_ms, drift_factor, next_available_time_ms)

Generate the next packet transmission time based on a Poisson process (exponentially distributed intervals),
while strictly enforcing a Duty Cycle constraint.

# Arguments
- `current_time_ms`: The current simulation time (or the time of the last event).
- `mean_interval_ms`: The mean interval for the Poisson process (1/lambda).
- `drift_factor`: The clock drift factor for the terminal (e.g., 1.000005).
- `next_available_time_ms`: The earliest time the terminal is allowed to transmit due to Duty Cycle.

# Returns
- `actual_next_tx`: The calculated next transmission time (Float64).
"""
function generate_next_poisson_time(current_time_ms::Float64, mean_interval_ms::Float64, drift_factor::Float64, next_available_time_ms::Float64)
    # 1. Generate Poisson interval (Exponential Distribution)
    # interval = -mean * ln(rand)
    # Apply clock drift to the interval measurement
    interval = -mean_interval_ms * log(rand()) * drift_factor
    
    # 2. Calculate target time based on the random interval
    target_time = current_time_ms + interval
    
    # 3. Enforce Duty Cycle Constraint (Clamping)
    # If the target time is earlier than the allowed DC time, we must wait.
    # This effectively creates a "Poisson process with dead time".
    actual_next_tx = max(target_time, next_available_time_ms)
    
    return actual_next_tx
end

end

# ============================================================
# ポアソン過程検証コード（このファイルを直接実行した場合のみ動作）
# ============================================================
if abspath(PROGRAM_FILE) == @__FILE__
    using Statistics
    
    println("="^60)
    println("   Poisson Process Verification")
    println("="^60)
    
    # パラメータ設定
    mean_interval_ms = 60000.0  # 平均60秒間隔
    drift_factor = 1.0          # ドリフトなし（検証のため）
    num_samples = 10000         # サンプル数
    simulation_duration_ms = 600000.0  # 10分間
    
    println("\nParameters:")
    println("  Mean Interval: $(mean_interval_ms) ms ($(mean_interval_ms/1000) seconds)")
    println("  Number of Samples: $(num_samples)")
    println("  Drift Factor: $(drift_factor)")
    println()
    
    # ============================================================
    # Test 1: 間隔の分布が指数分布に従うか
    # ============================================================
    println("Test 1: Interval Distribution (Exponential)")
    println("-"^60)
    
    intervals = Float64[]
    current_time = 0.0
    next_available = 0.0
    
    for i in 1:num_samples
        local next_time = PacketGeneration.generate_next_poisson_time(current_time, mean_interval_ms, drift_factor, next_available)
        local interval = next_time - current_time
        push!(intervals, interval)
        global current_time = next_time
        global next_available = current_time  # DC制約なし（純粋なポアソン過程）
    end
    
    # 統計量計算
    mean_interval = mean(intervals)
    std_interval = std(intervals)
    min_interval = minimum(intervals)
    max_interval = maximum(intervals)
    
    # 理論値（指数分布）
    theoretical_mean = mean_interval_ms
    theoretical_std = mean_interval_ms  # 指数分布では標準偏差 = 平均
    
    println("  Observed Statistics:")
    println("    Mean:   $(round(mean_interval, digits=2)) ms")
    println("    Std:    $(round(std_interval, digits=2)) ms")
    println("    Min:    $(round(min_interval, digits=2)) ms")
    println("    Max:    $(round(max_interval, digits=2)) ms")
    println()
    println("  Theoretical (Exponential Distribution):")
    println("    Mean:   $(round(theoretical_mean, digits=2)) ms")
    println("    Std:    $(round(theoretical_std, digits=2)) ms")
    println()
    
    # 誤差評価
    mean_error = abs(mean_interval - theoretical_mean) / theoretical_mean * 100
    std_error = abs(std_interval - theoretical_std) / theoretical_std * 100
    
    println("  Error:")
    println("    Mean Error: $(round(mean_error, digits=2)) %")
    println("    Std Error:  $(round(std_error, digits=2)) %")
    println()
    
    # 判定
    if mean_error < 5.0 && std_error < 10.0
        println("  ✅ PASS: Intervals follow exponential distribution")
    else
        println("  ❌ FAIL: Intervals do NOT follow exponential distribution")
    end
    println()
    
    # ============================================================
    # Test 2: 到着数の分布がポアソン分布に従うか
    # ============================================================
    println("Test 2: Arrival Count Distribution (Poisson)")
    println("-"^60)
    
    # 一定時間窓での到着数をカウント
    window_duration_ms = 60000.0  # 60秒窓
    num_windows = 100
    
    arrival_counts = Int[]
    
    for w in 1:num_windows
        local current_time = 0.0
        local next_available = 0.0
        local count = 0
        
        while current_time < window_duration_ms
            next_time = PacketGeneration.generate_next_poisson_time(current_time, mean_interval_ms, drift_factor, next_available)
            if next_time < window_duration_ms
                count += 1
                current_time = next_time
                next_available = current_time
            else
                break
            end
        end
        
        push!(arrival_counts, count)
    end
    
    # 統計量
    mean_arrivals = mean(arrival_counts)
    var_arrivals = var(arrival_counts)
    
    # 理論値（ポアソン分布では平均 = 分散）
    theoretical_lambda = window_duration_ms / mean_interval_ms
    
    println("  Window Duration: $(window_duration_ms) ms ($(window_duration_ms/1000) seconds)")
    println("  Number of Windows: $(num_windows)")
    println()
    println("  Observed Statistics:")
    println("    Mean Arrivals: $(round(mean_arrivals, digits=2))")
    println("    Variance:      $(round(var_arrivals, digits=2))")
    println()
    println("  Theoretical (Poisson Distribution):")
    println("    Lambda (Mean): $(round(theoretical_lambda, digits=2))")
    println("    Variance:      $(round(theoretical_lambda, digits=2))")
    println()
    
    # 誤差評価
    mean_arrival_error = abs(mean_arrivals - theoretical_lambda) / theoretical_lambda * 100
    var_ratio = var_arrivals / mean_arrivals  # ポアソン分布では1に近いはず
    
    println("  Error:")
    println("    Mean Error:      $(round(mean_arrival_error, digits=2)) %")
    println("    Variance/Mean:   $(round(var_ratio, digits=2)) (should be ≈ 1.0)")
    println()
    
    # 判定
    if mean_arrival_error < 10.0 && abs(var_ratio - 1.0) < 0.2
        println("  ✅ PASS: Arrival counts follow Poisson distribution")
    else
        println("  ❌ FAIL: Arrival counts do NOT follow Poisson distribution")
    end
    println()
    
    # ============================================================
    # Test 3: ヒストグラム表示（簡易版）
    # ============================================================
    println("Test 3: Interval Distribution Histogram")
    println("-"^60)
    
    # ビン分割（0-5倍平均まで、10ビン）
    num_bins = 10
    max_interval_plot = mean_interval_ms * 5
    bin_width = max_interval_plot / num_bins
    
    histogram = zeros(Int, num_bins)
    
    for interval in intervals
        if interval < max_interval_plot
            bin_idx = min(Int(floor(interval / bin_width)) + 1, num_bins)
            histogram[bin_idx] += 1
        end
    end
    
    # 正規化
    max_count = maximum(histogram)
    bar_width = 40
    
    println("  Interval Range (ms) | Count | Distribution")
    println("  " * "-"^56)
    
    for i in 1:num_bins
        range_start = (i-1) * bin_width
        range_end = i * bin_width
        count = histogram[i]
        bar_length = Int(round(count / max_count * bar_width))
        bar = "█"^bar_length
        
        println("  $(lpad(Int(round(range_start)), 7)) - $(lpad(Int(round(range_end)), 7)) | $(lpad(count, 5)) | $bar")
    end
    
    println()
    println("  Expected: Exponential decay (more samples at lower intervals)")
    println()
    
    # ============================================================
    # まとめ
    # ============================================================
    println("="^60)
    println("Summary:")
    println("  The packet generation follows a Poisson process if:")
    println("  1. Inter-arrival times are exponentially distributed")
    println("  2. Number of arrivals in fixed windows follows Poisson distribution")
    println()
    println("  Both tests should PASS for valid Poisson process.")
    println("="^60)
end
