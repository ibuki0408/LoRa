using Distributed

# 並列ワーカーの追加（CPUコア数に応じて自動調整）
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)  # 最大8ワーカー
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

# 1. 必要なファイルを読み込み
@everywhere include(joinpath(@__DIR__, "..", "src", "Prop_Async.jl"))

function evaluate_mean_per_parallel(num_terminals::Int, num_trials::Int)
    # パラメータ生成
    params = create_integrated_params()
    
    # 端末数を上書き
    params = IntegratedParameters(
        params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
        params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
        num_terminals,  # ← 指定された端末数
        params.area_size_m, params.slot_length_ms, params.packet_airtime_ms,
        params.cs_threshold_dbm,
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
        false, false, false, true, 10,  # ← ログ無効, force_async=true
        false  # use_waveform_phy
    )
    
    println("="^60)
    println("  Mean PER Evaluation (Parallel - Optimized)")
    println("="^60)
    @printf("Terminals:  %d\n", num_terminals)
    @printf("Trials:     %d\n", num_trials)
    @printf("Workers:    %d\n", nworkers())
    println("="^60)
    println()

    # 並列実行
    results_list = pmap(1:num_trials) do i
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
    total_success = sum([r["success"] for r in results_list])
    total_packets = sum([r["total_packets"] for r in results_list])

    mean_per = mean(per_list)
    std_per = std(per_list)
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
    println("="^60)

    # 出力ディレクトリ作成
    output_dir = joinpath(@__DIR__, "..", "results", "parallel_evaluation", "optimized")
    mkpath(output_dir)
    
    # CSV保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    out_file = joinpath(output_dir, "mean_per_$(num_terminals)term_$(num_trials)trials_$(timestamp).csv")
    open(out_file, "w") do io
        println(io, "trial,per")
        for i in 1:num_trials
            println(io, "$i,$(per_list[i])")
        end
        println(io, "MEAN,$(mean_per)")
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
    ENABLE_ACK = false     # 再送機能を有効にするかどうか (Fair comparison)
    MAX_RETRIES = 3         # 最大再送回数
    MAX_BUFFER_SIZE = 10    # バッファサイズ
    ENABLE_CARRIER_SENSE = true  # キャリアセンスを有効にするか
    ENABLE_CAPTURE_EFFECT = true # キャプチャ効果を有効にするか
    TARGET_SYNC_RATE = 0.95      # 期待される同期成功率 (Fastモード用)
    # ============================================================
    
    # 出力ディレクトリ作成
    out_dir = joinpath(@__DIR__, "..", "results", "parallel_evaluation", "Prop_Async_Optimized")
    mkpath(out_dir)
    
    # パラメータ設定をファイルに保存
    params = create_integrated_params()
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    ack_str = ENABLE_ACK ? "ACK" : "NoACK"
    param_file = joinpath(out_dir, "prop_async_optimized_parameters_SF$(SPREADING_FACTOR)_$(NUM_TRIALS)trials_$(ack_str)_$(timestamp).txt")
    open(param_file, "w") do io
        println(io, "="^60)
        println(io, "Optimized Async Simulation Parameters")
        println(io, "="^60)
        println(io, "Terminal Counts:      $(TERMINAL_COUNTS)")
        println(io, "Number of Trials:     $(NUM_TRIALS)")
        println(io, "SF:                   $(params.spreading_factor)")
        println(io, "Payload:              $(params.lora_payload_bytes) bytes")
        println(io, "Channels:             $(params.num_channels)")
        println(io, "Mean Interval:        $(params.mean_event_interval_ms) ms")
        println(io, "ACK Enabled:          $(params.enable_ack)")
        println(io, "Max Retries:          $(params.max_retries)")
        println(io, "LBT Threshold:        $(params.cs_threshold_dbm) dBm")
        println(io, "Duration:             $(params.simulation_duration_ms) ms")
        println(io, "="^60)
    end
    
    summary_data = []
    
    for num_terminals in TERMINAL_COUNTS
        println("\n" * "="^60)
        println("  Evaluating: $num_terminals terminals")
        println("="^60)
        
        params = create_integrated_params()
        lora_params = create_lora_params(SPREADING_FACTOR, params.lora_payload_bytes)
        airtime_ms = calculate_lora_airtime(lora_params)
        
        params = IntegratedParameters(
            params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
            params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
            num_terminals,
            params.area_size_m, SLOT_LENGTH_MS, airtime_ms,
            params.cs_threshold_dbm,
            SPREADING_FACTOR, params.lora_payload_bytes, params.num_channels,
            BEACON_INTERVAL_MS, params.simulation_duration_ms,
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
            false, false, false, true, MAX_BUFFER_SIZE,
            ENABLE_CARRIER_SENSE, ENABLE_CAPTURE_EFFECT, TARGET_SYNC_RATE,
            false  # use_waveform_phy
        )
        
        results_list = pmap(1:NUM_TRIALS) do i
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
        
        per_list = [r["per"] for r in results_list]
        original_per_list = [r["original_per"] for r in results_list]
        total_success = sum([r["success"] for r in results_list])
        total_packets = sum([r["total_packets"] for r in results_list])
        throughput_bps_list = [r["throughput_bps"] for r in results_list]
        norm_throughput_list = [r["norm_throughput"] for r in results_list]
        total_generated_list = [r["total_generated"] for r in results_list]
        total_packets_list = [r["total_packets"] for r in results_list]
        
        mean_per = mean(per_list)
        mean_original_per = mean(original_per_list)
        mean_throughput_bps = mean(throughput_bps_list)
        mean_norm_throughput = mean(norm_throughput_list)
        mean_total_generated = mean(total_generated_list)
        mean_total_packets = mean(total_packets_list)
        
        push!(summary_data, (
            num_terminals = num_terminals,
            mean_per = mean_per,
            mean_original_per = mean_original_per,
            mean_throughput_bps = mean_throughput_bps,
            mean_norm_throughput = mean_norm_throughput,
            mean_total_generated = mean_total_generated,
            mean_total_packets = mean_total_packets
        ))
        
        @printf("  Mean PER (Total):    %.6f %%\n", mean_per)
        @printf("  Mean PER (Original): %.6f %%\n", mean_original_per)
        @printf("  Throughput (bps):    %.2f\n", mean_throughput_bps)
    end
    
    # サマリーCSVを保存
    summary_file = joinpath(out_dir, "Optimized_Async_summary_SF$(SPREADING_FACTOR)_$(ack_str)_$(timestamp).csv")
    open(summary_file, "w") do io
        println(io, "num_terminals,mean_per,mean_original_per,mean_throughput_bps,mean_norm_throughput,mean_total_generated,mean_total_packets")
        for data in summary_data
            @printf(io, "%d,%.6f,%.6f,%.2f,%.6f,%.2f,%.2f\n",
                    data.num_terminals, data.mean_per, data.mean_original_per, 
                    data.mean_throughput_bps, data.mean_norm_throughput,
                    data.mean_total_generated, data.mean_total_packets)
        end
    end
    println("\nSummary saved to: $summary_file")
end
