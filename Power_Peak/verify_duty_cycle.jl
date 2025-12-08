using CSV, DataFrames, Statistics

# 最新の送信ログを読み込み
function verify_duty_cycle(log_file::String, sim_duration_ms::Float64, airtime_ms::Float64)
    println("="^60)
    println("Duty Cycle Verification: $log_file")
    println("="^60)
    
    # CSVファイルを探す
    if !isfile(log_file)
        println("Error: File not found: $log_file")
        return
    end
    
    df = CSV.read(log_file, DataFrame)
    
    # 端末ごとに集計
    terminal_ids = sort(unique(df.terminal_id))
    
    println("\nPer-Terminal Duty Cycle Analysis:")
    println("-"^60)
    
    dc_violations = 0
    total_terminals = length(terminal_ids)
    
    for tid in terminal_ids
        term_df = filter(row -> row.terminal_id == tid, df)
        
        # 送信回数
        num_tx = nrow(term_df)
        
        # 総送信時間
        total_tx_time = num_tx * airtime_ms
        
        # 実効Duty Cycle
        actual_dc = total_tx_time / sim_duration_ms * 100
        
        # 1%を超えているかチェック
        is_violation = actual_dc > 1.01  # 1%に若干のマージン
        
        if is_violation
            dc_violations += 1
            println("Terminal $tid: $(num_tx) pkts, DC=$(round(actual_dc, digits=3))% ⚠️ VIOLATION")
        else
            println("Terminal $tid: $(num_tx) pkts, DC=$(round(actual_dc, digits=3))% ✓")
        end
    end
    
    println("-"^60)
    println("Summary:")
    println("  Total Terminals: $total_terminals")
    println("  DC Violations: $dc_violations")
    println("  Compliance Rate: $(round((total_terminals - dc_violations) / total_terminals * 100, digits=2))%")
    println("="^60)
end

# パラメータ
sim_duration_ms = 600000.0  # 10分
airtime_ms = 288.77  # SF10, 10 bytes

# 最新のログファイルを探す
function find_latest_log(dir::String, pattern::String)
    files = readdir(dir, join=true)
    matching = filter(f -> occursin(pattern, f) && endswith(f, ".csv"), files)
    if isempty(matching)
        return nothing
    end
    return sort(matching, by=mtime, rev=true)[1]
end

# CSMA/CA
csma_log = find_latest_log("result_csma_ca", "csma_ca_summary")
if csma_log !== nothing
    println("\nNote: CSMA/CA summary file found, but we need the detailed tx_log.")
    println("Searching for detailed logs...")
end

# Integrated
integrated_log = find_latest_log("result_integrated", "integrated_tx_log")
if integrated_log !== nothing
    verify_duty_cycle(integrated_log, sim_duration_ms, airtime_ms)
else
    println("No integrated tx log found")
end

println("\n")
println("Note: CSMA/CA simulation does not save detailed tx logs by default.")
println("We can only verify from the integrated simulation logs.")
