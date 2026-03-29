using Distributed

# 並列ワーカーの追加（CPUコア数に応じて自動調整）
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)  # 最大8ワーカー
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

# 1. 必要なファイルを読み込み
@everywhere include("Prop_Async.jl")

function evaluate_mean_per_parallel(num_terminals::Int, num_trials::Int)
    # パラメータ生成
    params = create_integrated_params()
    
    # 端末数を上書き
    params = IntegratedParameters(
        params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
        params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
        num_terminals,  # ← 指定された端末数
        params.area_size_m, params.slot_length_ms, params.packet_airtime_ms,
        params.enable_carrier_sense, params.cs_threshold_dbm,
        params.spreading_factor, params.lora_payload_bytes, params.num_channels,
        params.beacon_interval_ms, params.simulation_duration_ms,
        params.max_startup_delay_ms,
        params.duty_cycle, params.mean_event_interval_ms,
        params.shadowing_enabled, params.shadowing_std_db, params.pass_loss_exp,
        params.sync_center_freq_ghz, params.data_center_freq_ghz, params.reference_path_loss_db,
        params.sync_bs_x_m, params.sync_bs_y_m,
        params.sync_observation_duration_ms, params.gw_tx_power_dbm,
        params.noise_floor_window_ms, params.detection_margin_db,
        params.min_samples, params.debounce_time_ms, params.initial_window_duration_ms,
        params.tx_jitter_max_ms,
        params.enable_intermittent_rx, params.intermittent_window_ms, params.initial_search_duration_ms,
        params.use_deterministic_jitter, params.num_jitter_offsets, params.deterministic_jitter_random_ms,
        params.collision_model,
        params.enable_ack, params.max_retries, params.ack_timeout_ms,
        params.rx1_delay_ms, params.backoff_base_ms,
        params.lbt_duration_ms, params.lbt_sample_interval_ms,
        false, false, false, false, 10  # ← ログ・ファイル出力内容・プロット無効, force_async, buffer=10
    )
    
    println("="^60)
    println("  Mean PER Evaluation (Parallel)")
    println("="^60)
    @printf("Terminals:  %d\n", num_terminals)
    @printf("Trials:     %d\n", num_trials)
    @printf("Workers:    %d\n", nworkers())
    println("="^60)
    println()

    # 並列実行
    results_list = pmap(1:num_trials) do i
        return run_integrated_simulation_with_params(params)
    end

    # 結果の集計
    per_list = [r["per"] for r in results_list]
    sync_rate_list = [r["sync_success_rate"] for r in results_list]
    total_success = sum([r["success"] for r in results_list])
    total_packets = sum([r["total_packets"] for r in results_list])

    mean_per = mean(per_list)
    std_per = std(per_list)
    mean_sync = mean(sync_rate_list)
    overall_per = (1.0 - (total_success / total_packets))

    # 結果表示
    println("\n" * "="^60)
    println("        FINAL STATISTICS")
    println("="^60)
    @printf("Terminals:           %d\n", num_terminals)
    @printf("Number of Trials:    %d\n", num_trials)
    @printf("Mean PER:            %.6f\n", mean_per)
    @printf("Std Dev of PER:      %.6f\n", std_per)
    @printf("Overall PER (Total): %.6f\n", overall_per)
    @printf("Mean Sync Rate:      %.2f %%\n", mean_sync)
    println("="^60)

    # 出力ディレクトリ作成
    output_dir = "result_parallel_evaluation/integrated"
    mkpath(output_dir)
    
    # CSV保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    out_file = joinpath(output_dir, "mean_per_$(num_terminals)term_$(num_trials)trials_$(timestamp).csv")
    open(out_file, "w") do io
        println(io, "trial,per,sync_rate")
        for i in 1:num_trials
            println(io, "$i,$(per_list[i]),$(sync_rate_list[i])")
        end
        println(io, "MEAN,$(mean_per),$(mean_sync)")
    end
    println("Results saved to: $out_file")
    println("="^60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    # ============================================================
    # 設定: ここで評価パラメータを指定
    # ============================================================
    TERMINAL_COUNTS = [100,200,300,400,500]  # 評価する端末数のリスト
    NUM_TRIALS = 100         # 各設定での試行回数
    SPREADING_FACTOR = 8    # 拡散率 (7-12)
    SLOT_LENGTH_MS = 100.0   # スロット長 (ms)
    BEACON_INTERVAL_MS = 100.0  # ビーコン間隔 (ms)
    ENABLE_ACK = false      # 再送機能を有効にするかどうか
    MAX_RETRIES = 3         # 最大再送回数
    MAX_BUFFER_SIZE = 10    # バッファサイズ (待ち行列の最大数)
    # ============================================================
    
    # 出力ディレクトリ作成
    out_dir = "result_parallel_evaluation/Prop_Async"
    mkpath(out_dir)
    
    # パラメータ設定をファイルに保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    ack_str = ENABLE_ACK ? "ACK" : "NoACK"
    param_file = joinpath(out_dir, "prop_async_parameters_SF$(SPREADING_FACTOR)_Slot$(Int(SLOT_LENGTH_MS))ms_$(NUM_TRIALS)trials_$(ack_str)_$(timestamp).txt")
    open(param_file, "w") do io
        println(io, "="^60)
        println(io, "Simulation Parameters")
        println(io, "Execution Time: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "="^60)
        println(io, "")
        println(io, "[Evaluation Settings]")
        println(io, "  Terminal Counts:      $(TERMINAL_COUNTS)")
        println(io, "  Number of Trials:     $(NUM_TRIALS)")
        println(io, "")
        println(io, "[LoRa Parameters]")
        println(io, "  Spreading Factor:     $(SPREADING_FACTOR)")
        println(io, "  Payload Size:         10 bytes")
        println(io, "  Bandwidth:            125 kHz")
        println(io, "  Number of Channels:   8")
        println(io, "")
        println(io, "[MAC Layer]")
        println(io, "  Slot Length:          $(SLOT_LENGTH_MS) ms")
        println(io, "  Beacon Interval:      $(BEACON_INTERVAL_MS) ms")
        println(io, "  Initial Window:       $(BEACON_INTERVAL_MS * 2.0) ms (auto-calculated)")
        println(io, "  TX Jitter Max:        50.0 ms")
        println(io, "  Enable ACK:           $(ENABLE_ACK)")
        println(io, "  Max Retries:          $(MAX_RETRIES)")
        println(io, "  Max Buffer Size:      $(MAX_BUFFER_SIZE)")
        println(io, "  LBT Enabled:          true")
        println(io, "  LBT Threshold:        -80.0 dBm")
        println(io, "")
        println(io, "[Environment]")
        println(io, "  Area Size:            2000.0 m")
        println(io, "  Deployment Radius:    1000.0 m")
        println(io, "  Path Loss Exponent:   2.7")
        println(io, "  Shadowing Std Dev:    8.0 dB")
        println(io, "  TX Power:             13.0 dBm")
        println(io, "  Noise Figure:         6.0 dB")
        println(io, "")
        println(io, "[Simulation Control]")
        println(io, "  Duration:             3600000.0 ms (60 min)")
        println(io, "  Max Startup Delay:    60000.0 ms")
        println(io, "  Duty Cycle:           1.0 %")
        println(io, "  Mean Event Interval:  60000.0 ms")
        println(io, "="^60)
    end
    println("Parameter file saved: $param_file")
    
    # サマリー用のデータ構造
    summary_data = []
    
    # 各端末数について評価
    for num_terminals in TERMINAL_COUNTS
        println("\n" * "="^60)
        println("  Evaluating: $num_terminals terminals")
        println("  SF: $SPREADING_FACTOR, Slot: $(SLOT_LENGTH_MS)ms, Beacon: $(BEACON_INTERVAL_MS)ms")
        println("="^60)
        
        # パラメータ生成
        params = create_integrated_params()
        
        # ToA計算（SFに応じて）
        lora_params = create_lora_params(SPREADING_FACTOR, params.lora_payload_bytes)
        airtime_ms = calculate_lora_airtime(lora_params)
        
        # パラメータ上書き
        params = IntegratedParameters(
            params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
            params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
            num_terminals,  # ← 指定された端末数
            params.area_size_m, SLOT_LENGTH_MS, airtime_ms,  # ← スロット長とToAを上書き
            params.enable_carrier_sense, params.cs_threshold_dbm,
            SPREADING_FACTOR, params.lora_payload_bytes, params.num_channels,  # ← SFを上書き
            BEACON_INTERVAL_MS, params.simulation_duration_ms,  # ← ビーコン間隔を上書き
            params.max_startup_delay_ms,
            params.duty_cycle, params.mean_event_interval_ms,
            params.shadowing_enabled, params.shadowing_std_db, params.pass_loss_exp,
            params.sync_center_freq_ghz, params.data_center_freq_ghz, params.reference_path_loss_db,
            params.sync_bs_x_m, params.sync_bs_y_m,
            params.sync_observation_duration_ms, params.gw_tx_power_dbm,
            params.noise_floor_window_ms, params.detection_margin_db,
            params.min_samples, params.debounce_time_ms, params.initial_window_duration_ms,
            params.tx_jitter_max_ms,
            params.enable_intermittent_rx, params.intermittent_window_ms, params.initial_search_duration_ms,
            params.use_deterministic_jitter, params.num_jitter_offsets, params.deterministic_jitter_random_ms,
            params.collision_model,
            ENABLE_ACK, MAX_RETRIES, params.ack_timeout_ms,
            params.rx1_delay_ms, params.backoff_base_ms,
            params.lbt_duration_ms, params.lbt_sample_interval_ms,
            false, false, false, false, MAX_BUFFER_SIZE  # ← ログ無効, force_async=false, buffer上書き
        )
        
        println("Running $NUM_TRIALS trials in parallel...")
        println("  ToA: $(round(airtime_ms, digits=2))ms")
        
        # 並列実行（出力を抑制）
        results_list = pmap(1:NUM_TRIALS) do i
            # 各ワーカーの標準出力を抑制
            original_stdout = stdout
            redirect_stdout(devnull)
            
            result = nothing
            try
                result = run_integrated_simulation_with_params(params)
            finally
                redirect_stdout(original_stdout)
            end
            
            return result
        end
        
        # 結果の集計
        per_list = [r["per"] for r in results_list]
        sync_rate_list = [r["sync_success_rate"] for r in results_list]
        synced_per_list = [r["synced_per"] for r in results_list]
        async_per_list = [r["async_per"] for r in results_list]
        original_per_list = [r["original_per"] for r in results_list]
        residual_per_list = [r["final_reliability_per"] for r in results_list]  # Residual PER（最終不達率）
        total_success = sum([r["success"] for r in results_list])
        total_packets = sum([r["total_packets"] for r in results_list])
        total_collisions = sum([r["collisions"] for r in results_list])
        throughput_bps_list = [r["throughput_bps"] for r in results_list]
        norm_throughput_list = [r["norm_throughput"] for r in results_list]
        
        # Buffer Drop統計
        total_buffer_drops = sum([r["buffer_drops"] for r in results_list])
        total_generated = sum([r["total_generated"] for r in results_list])
        buffer_drop_rate = total_generated > 0 ? (total_buffer_drops / total_generated) : 0.0
        
        # 再送統計
        total_retransmissions = sum([r["total_retransmissions"] for r in results_list])
        retransmission_rate = total_packets > 0 ? (total_retransmissions / total_packets) : 0.0
        
        mean_per = mean(per_list)
        std_per = std(per_list)
        mean_sync = mean(sync_rate_list)
        std_sync = std(sync_rate_list)
        mean_synced_per = mean(synced_per_list)
        std_synced_per = std(synced_per_list)
        mean_async_per = mean(async_per_list)
        std_async_per = std(async_per_list)
        mean_original_per = mean(original_per_list)
        std_original_per = std(original_per_list)
        mean_residual_per = mean(residual_per_list)  # Residual PER平均
        std_residual_per = std(residual_per_list)
        overall_per = (1.0 - (total_success / total_packets))
        
        # 無駄な送信率の計算
        wasted_tx_rate = (total_collisions / total_packets)
        
        # Overall Throughput (全試行の成功パケット数から算出)
        # S = (成功パケット合計 * Airtime) / (シミュレーション時間 * チャネル数 * 試行回数)
        overall_norm_throughput = (total_success * airtime_ms) / (params.simulation_duration_ms * params.num_channels * NUM_TRIALS)
        # bps = (成功ビット合計) / (シミュレーション時間合計)
        overall_throughput_bps = (total_success * params.lora_payload_bytes * 8) / (params.simulation_duration_ms * NUM_TRIALS / 1000.0)
        
        mean_throughput_bps = mean(throughput_bps_list)
        std_throughput_bps = std(throughput_bps_list)
        mean_norm_throughput = mean(norm_throughput_list)
        std_norm_throughput = std(norm_throughput_list)

        # パラメータ設定 (force_async_mode = true に設定)
        params = create_integrated_params()
        params.num_terminals = num_terminals
        params.spreading_factor = SPREADING_FACTOR
        params.slot_length_ms = SLOT_LENGTH_MS
        params.beacon_interval_ms = BEACON_INTERVAL_MS
        params.enable_ack = ENABLE_ACK
        params.max_retries = MAX_RETRIES
        params.max_buffer_size = MAX_BUFFER_SIZE
        params.force_async_mode = true
        params.enable_detailed_logs = false # 高速化のため無効
        # サマリーデータに追加（簡略版）
        push!(summary_data, (
            num_terminals = num_terminals,
            num_trials = NUM_TRIALS,
            mean_per = mean_per,
            std_per = std_per,
            mean_original_per = mean_original_per,
            std_original_per = std_original_per,
            mean_sync = mean_sync,
            mean_throughput_bps = mean_throughput_bps,
            std_throughput_bps = std_throughput_bps,
            mean_norm_throughput = mean_norm_throughput,
            std_norm_throughput = std_norm_throughput,
            total_buffer_drops = total_buffer_drops,
            total_generated = total_generated,
            buffer_drop_rate = buffer_drop_rate,
            total_retransmissions = total_retransmissions,
            retransmission_rate = retransmission_rate
        ))
        
        # 結果表示
        println("\nResults for $num_terminals terminals:")
        @printf("  Mean PER (Total):    %.6f %%\n", mean_per)
        @printf("  Mean PER (Original): %.6f %%\n", mean_original_per)
        @printf("  Overall PER (Total): %.6f %%\n", overall_per)
        @printf("  Throughput (bps):    %.2f (Mean) | %.2f (Overall)\n", mean_throughput_bps, overall_throughput_bps)
        @printf("  Throughput (Norm):   %.5f (Mean) | %.5f (Overall)\n", mean_norm_throughput, overall_norm_throughput)
        @printf("  Mean Sync Rate:      %.2f ± %.2f %%\n", mean_sync, std_sync)
        println("  ├─ Synced PER:       $(round(mean_synced_per, digits=2))% ± $(round(std_synced_per, digits=2))%")
        println("  └─ Async PER:        $(round(mean_async_per, digits=2))% ± $(round(std_async_per, digits=2))%")
    end
    
    # 全体のサマリーを表示
    println("\n" * "="^60)
    println("  SUMMARY: All Terminal Counts")
    println("  SF: $SPREADING_FACTOR, Slot: $(SLOT_LENGTH_MS)ms, Beacon: $(BEACON_INTERVAL_MS)ms")
    println("="^60)
    println("Terminals | Mean PER (%) | Original PER (%) | Thr (bps) | Sync Rate (%)")
    println("-"^75)
    for data in summary_data
        @printf("%9d | %12.4f | %16.4f | %9.1f | %13.2f\n", 
                data.num_terminals, data.mean_per, data.mean_original_per, data.mean_throughput_bps, data.mean_sync)
    end
    println("="^60)
    
    # サマリーCSVを保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    ack_str = ENABLE_ACK ? "ACK" : "NoACK"
    summary_file = joinpath(out_dir, "Async_summary_SF$(SPREADING_FACTOR)_slot$(Int(SLOT_LENGTH_MS))ms_beacon$(Int(BEACON_INTERVAL_MS))ms_$(ack_str)_$(NUM_TRIALS)trials_$(timestamp).csv")
    open(summary_file, "w") do io
        cols = ["num_terminals", "num_trials", 
                "mean_per", "std_per", 
                "mean_original_per", "std_original_per",
                "mean_sync",
                "mean_throughput_bps", "std_throughput_bps",
                "mean_norm_throughput", "std_norm_throughput",
                "total_buffer_drops", "total_generated", "buffer_drop_rate",
                "total_retransmissions", "retransmission_rate"]
        println(io, join(cols, ","))
        for data in summary_data
            @printf(io, "%d,%d,%.6f,%.6f,%.6f,%.6f,%.2f,%.2f,%.2f,%.6f,%.6f,%d,%d,%.6f,%d,%.6f\n",
                    data.num_terminals, data.num_trials,
                    data.mean_per, data.std_per,
                    data.mean_original_per, data.std_original_per,
                    data.mean_sync,
                    data.mean_throughput_bps, data.std_throughput_bps,
                    data.mean_norm_throughput, data.std_norm_throughput,
                    data.total_buffer_drops, data.total_generated, data.buffer_drop_rate,
                    data.total_retransmissions, data.retransmission_rate)
        end
    end
    println("\nSummary saved to: $summary_file")
    println("="^60)
end