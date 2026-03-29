using Distributed

# 並列ワーカーの追加
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

# 全ワーカーで環境構築
@everywhere using Statistics, Printf, Dates, Random
@everywhere include("lbt_aloha_sim.jl")

function run_lbt_batch_evaluation()
    # ============================================================
    # 設定: ここで評価パラメータを指定
    # ============================================================
    TERMINAL_COUNTS = [100,200,300,400,500]        # 評価したい端末数リスト
    NUM_TRIALS = 100                # 試行回数 (平均を取るため)
    AREA_SIZE_M = 2000.0           # エリアサイズ (m)
    SF = 8                         # 拡散率
    # ============================================================

    println("="^60)
    println("   LBT-ALOHA Parallel Evaluation")
    println("   SF: $SF, Trials: $NUM_TRIALS")
    println("="^60)

    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    summary_file = "result_parallel_evaluation/LBT/summary_SF$(SF)_$(timestamp).csv"
    mkpath(dirname(summary_file))

    open(summary_file, "w") do io
        println(io, "num_terminals,mean_per,std_per,mean_throughput_bps")
        
        for num_terminals in TERMINAL_COUNTS
            println("\nEvaluating $num_terminals terminals...")
            
            # 並列実行
            results = pmap(1:NUM_TRIALS) do trial_id
                # 乱数シードの設定
                Random.seed!(trial_id + num_terminals * 1000)
                
                # パラメータ生成 (lbt_aloha_sim.jl の関数を使用)
                p = create_csma_params()
                
                # パラメータを上書き
                p = CSMAParameters(
                    p.tx_power_dbm, p.noise_figure_db,
                    num_terminals, AREA_SIZE_M,
                    SF, p.payload_bytes, p.num_channels,
                    p.simulation_duration_ms, p.duty_cycle,
                    p.enable_carrier_sense, p.mean_event_interval_ms,
                    p.max_startup_delay_ms,
                    p.pass_loss_exp, p.shadowing_std_db, p.reference_path_loss_db,
                    p.lbt_duration_ms, p.lbt_sample_interval_ms
                )
                
                # 標準出力を抑制して実行
                res = nothing
                redirect_stdout(devnull) do
                    # lbt_aloha_sim.jl の実行関数
                    res = run_csma_ca_with_params(p)
                end
                return res
            end
            
            # 集計
            pers = [r["per"] for r in results]
            throughputs = [r["throughput_bps"] for r in results]
            
            m_per = mean(pers)
            s_per = std(pers)
            m_thru = mean(throughputs)
            
            @printf("  Mean PER: %.2f %% (std: %.2f)\n", m_per, s_per)
            @printf("  Mean Throughput: %.2f bps\n", m_thru)
            
            println(io, "$num_terminals,$m_per,$s_per,$m_thru")
            flush(io)
        end
    end
    
    println("\nSummary saved to: $summary_file")
end

run_lbt_batch_evaluation()
