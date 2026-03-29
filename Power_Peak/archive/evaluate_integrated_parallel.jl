using Distributed

# 並列ワーカーの追加
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

# 全ワーカーで環境構築
@everywhere using Statistics, Printf, Dates, Random
@everywhere include("integrated_lora_sim.jl")

function run_batch_evaluation()
    # ============================================================
    # 設定: ここで評価パラメータを指定
    # ============================================================
    TERMINAL_COUNTS = [100, 200, 300, 400, 500]
    NUM_TRIALS = 100
    AREA_SIZE_M = 2000.0  # 半径 1000m
    SF = 10
    SLOT_LENGTH_MS = 300.0
    BEACON_INTERVAL_MS = 300.0
    # ============================================================

    println("="^60)
    println("   Integrated Simulation Parallel Evaluation")
    println("   Radius: $(AREA_SIZE_M/2)m, SF: $SF, Trials: $NUM_TRIALS")
    println("="^60)

    summary_file = "result_parallel_evaluation/integrated/summary_radius1km_SF$(SF)_$(Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")).csv"
    mkpath(dirname(summary_file))

    open(summary_file, "w") do io
        println(io, "num_terminals,mean_per,std_per,sync_rate")
        
        for num_terminals in TERMINAL_COUNTS
            println("\nEvaluating $num_terminals terminals...")
            
            # 並列実行
            results = pmap(1:NUM_TRIALS) do trial_id
                # 乱数シードの設定
                Random.seed!(trial_id + num_terminals * 1000)
                
                # パラメータ生成
                p = create_integrated_params()
                
                # 重要パラメータを上書き
                p = IntegratedParameters(
                    p.signal_duration_us, p.signal_bw_mhz, p.terminal_bw_mhz,
                    p.tx_sampling_rate_mhz, p.rx_sampling_rate_mhz, p.tx_power_dbm, p.noise_figure_db,
                    num_terminals, AREA_SIZE_M, SLOT_LENGTH_MS, p.packet_airtime_ms,
                    p.transmission_prob, p.enable_carrier_sense, p.cs_threshold_dbm,
                    SF, p.lora_payload_bytes, p.num_channels,
                    BEACON_INTERVAL_MS, p.simulation_duration_ms,
                    p.max_startup_delay_ms, p.duty_cycle, p.mean_event_interval_ms,
                    p.shadowing_enabled, p.shadowing_std_db, p.pass_loss_exp,
                    p.sync_center_freq_ghz, p.data_center_freq_ghz, p.reference_path_loss_db,
                    p.sync_bs_x_m, p.sync_bs_y_m,
                    p.sync_observation_duration_ms, p.gw_tx_power_dbm,
                    p.noise_floor_window_ms, p.detection_margin_db,
                    p.min_samples, p.debounce_time_ms, p.initial_window_duration_ms,
                    p.tx_jitter_max_ms,
                    p.enable_intermittent_rx, p.intermittent_window_ms, p.initial_search_duration_ms,
                    p.use_deterministic_jitter, p.num_jitter_offsets, p.deterministic_jitter_random_ms,
                    p.collision_model,
                    p.enable_ack, p.max_retries, p.ack_timeout_ms,
                    p.rx1_delay_ms, p.backoff_base_ms,
                    p.lbt_duration_ms, p.lbt_sample_interval_ms,
                    false, false, false # 無効化: file, plot, detailed_logs
                )
                
                # 標準出力を抑制して実行
                res = nothing
                redirect_stdout(devnull) do
                    res = run_integrated_simulation_with_params(p)
                end
                return res
            end
            
            # 集計
            pers = [r["per"] for r in results]
            sync_rates = [r["sync_success_rate"] for r in results]
            
            m_per = mean(pers)
            s_per = std(pers)
            m_sync = mean(sync_rates)
            
            @printf("  Mean PER: %.2f %% (std: %.2f)\n", m_per, s_per)
            @printf("  Sync Rate: %.2f %%\n", m_sync)
            
            println(io, "$num_terminals,$m_per,$s_per,$m_sync")
            flush(io)
        end
    end
    
    println("\nSummary saved to: $summary_file")
end

run_batch_evaluation()
