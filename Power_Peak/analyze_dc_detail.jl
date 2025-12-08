using CSV, DataFrames, Statistics

# より詳細な分析
function detailed_dc_analysis(log_file::String)
    df = CSV.read(log_file, DataFrame)
    
    println("="^60)
    println("Detailed Duty Cycle Analysis")
    println("="^60)
    
    # 端末1の詳細を見る
    term1 = filter(row -> row.terminal_id == 1, df)
    sort!(term1, :start_ms)
    
    println("\nTerminal 1 Transmission Timeline (first 5 packets):")
    println("-"^60)
    
    for i in 1:min(5, nrow(term1))
        row = term1[i, :]
        duration = row.end_ms - row.start_ms
        
        if i > 1
            prev_end = term1[i-1, :end_ms]
            interval = row.start_ms - prev_end
            println("Packet $i:")
            println("  Start: $(round(row.start_ms, digits=2)) ms")
            println("  End: $(round(row.end_ms, digits=2)) ms")
            println("  Duration: $(round(duration, digits=2)) ms")
            println("  Interval from prev: $(round(interval, digits=2)) ms")
            println("  Expected off-period (99x duration): $(round(duration * 99, digits=2)) ms")
            println()
        else
            println("Packet $i:")
            println("  Start: $(round(row.start_ms, digits=2)) ms")
            println("  End: $(round(row.end_ms, digits=2)) ms")
            println("  Duration: $(round(duration, digits=2)) ms")
            println()
        end
    end
    
    # 理論値との比較
    println("-"^60)
    println("Theoretical Analysis:")
    println("-"^60)
    
    airtime_ms = 288.77
    duty_cycle = 0.01
    
    # 理論的な送信周期
    theoretical_period = airtime_ms / duty_cycle
    println("ToA (airtime): $airtime_ms ms")
    println("Duty Cycle: $(duty_cycle * 100)%")
    println("Theoretical Tx Period: $(round(theoretical_period, digits=2)) ms")
    
    # 実際の平均間隔
    if nrow(term1) > 1
        intervals = Float64[]
        for i in 2:nrow(term1)
            push!(intervals, term1[i, :start_ms] - term1[i-1, :start_ms])
        end
        avg_interval = mean(intervals)
        println("Actual Avg Tx Interval: $(round(avg_interval, digits=2)) ms")
        println("Difference: $(round(avg_interval - theoretical_period, digits=2)) ms")
        
        # 実効DC
        actual_dc = airtime_ms / avg_interval * 100
        println("Actual DC (from interval): $(round(actual_dc, digits=3))%")
    end
    
    println("="^60)
end

# 実行
log_file = "result_integrated/integrated_tx_log_2025-12-08_10-23-34.csv"
if isfile(log_file)
    detailed_dc_analysis(log_file)
else
    println("Log file not found: $log_file")
end
