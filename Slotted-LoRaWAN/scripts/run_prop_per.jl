using Distributed

# 並列ワーカーの追加（CPUコア数に応じて自動調整）
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)  # 最大8ワーカー
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

# 全ワーカーで環境構築
@everywhere using Statistics, Printf, Dates
@everywhere include(joinpath(@__DIR__, "..", "src", "Prop.jl"))

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
        false, false, false, false, 10,  # ← ログ・ファイル出力内容・プロット無効, force_async, buffer=10
        false  # use_waveform_phy: 電力モデル（高速）
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
    
    # 追加: 同期/非同期それぞれのPER
    synced_per_list = [r["synced_per"] for r in results_list]
    async_per_list = [r["async_per"] for r in results_list]
    
    total_success = sum([r["success"] for r in results_list])
    total_packets = sum([r["total_packets"] for r in results_list])

    mean_per = mean(per_list)
    std_per = std(per_list)
    mean_sync = mean(sync_rate_list)
    mean_synced_per = mean(synced_per_list) # 追加
    mean_async_per = mean(async_per_list)   # 追加
    overall_per = (1.0 - (total_success / total_packets))

    # 結果表示
    println("\n" * "="^60)
    println("        FINAL STATISTICS")
    println("="^60)
    @printf("Terminals:           %d\n", num_terminals)
    @printf("Number of Trials:    %d\n", num_trials)
    @printf("Mean PER:            %.6f\n", mean_per)
    @printf("  - Synced PER:      %.6f\n", mean_synced_per) # 追加
    @printf("  - Async PER:       %.6f\n", mean_async_per) # 追加
    @printf("Std Dev of PER:      %.6f\n", std_per)
    @printf("Overall PER (Total): %.6f\n", overall_per)
    @printf("Mean Sync Rate:      %.2f %%\n", mean_sync)
    println("="^60)

    # 出力ディレクトリ作成
    output_dir = joinpath(@__DIR__, "..", "results", "parallel_evaluation", "Prop")
    mkpath(output_dir)
    
    # CSV保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    out_file = joinpath(output_dir, "mean_per_$(num_terminals)term_$(num_trials)trials_$(timestamp).csv")
    open(out_file, "w") do io
        println(io, "trial,per,sync_rate,synced_per,async_per") # ヘッダー更新
        for i in 1:num_trials
            println(io, "$i,$(per_list[i]),$(sync_rate_list[i]),$(synced_per_list[i]),$(async_per_list[i])")
        end
        println(io, "MEAN,$(mean_per),$(mean_sync),$(mean_synced_per),$(mean_async_per)")
    end
    println("Results saved to: $out_file")
    println("="^60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    # ============================================================
    # 設定: ここで評価パラメータを指定
    # ============================================================
    TERMINAL_COUNTS = [100,200,300,400,500]  # 評価する端末数のリスト
    NUM_TRIALS = 8         # 各設定での試行回数
    SPREADING_FACTOR = 8    # 拡散率 (7-12)
    SLOT_LENGTH_MS = 200.0   # スロット長 (ms)
    BEACON_INTERVAL_MS = 100.0  # ビーコン間隔 (ms)
    ENABLE_ACK = false       # 再送機能を有効にするかどうか
    MAX_RETRIES = 3         # 最大再送回数
    MAX_BUFFER_SIZE = 10    # バッファサイズ (待ち行列の最大数)
    ENABLE_CARRIER_SENSE = true  # キャリアセンスを有効にするか
    ENABLE_CAPTURE_EFFECT = false # キャプチャ効果を有効にするか
    TARGET_SYNC_RATE = 1.0     # 期待される同期成功率 (Fastモード用)
    # ============================================================
    
    # 出力ディレクトリ作成
    output_dir = joinpath(@__DIR__, "..", "results", "parallel_evaluation", "Prop")
    mkpath(output_dir)
    
    # パラメータ設定をファイルに保存（ログ用に先行してparamsを生成）
    params = create_integrated_params()
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    ack_str = ENABLE_ACK ? "ACK" : "NoACK"
    param_file = joinpath(output_dir, "prop_parameters_SF$(SPREADING_FACTOR)_Slot$(Int(SLOT_LENGTH_MS))ms_$(NUM_TRIALS)trials_$(ack_str)_$(timestamp).txt")
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
        println(io, "  Slot Length:          $(params.slot_length_ms) ms")
        println(io, "  Beacon Interval:      $(params.beacon_interval_ms) ms")
        println(io, "  Initial Window:       $(params.initial_window_duration_ms) ms")
        println(io, "  TX Jitter Max:        $(params.tx_jitter_max_ms) ms")
        println(io, "  Enable ACK:           $(params.enable_ack)")
        println(io, "  Max Retries:          $(params.max_retries)")
        println(io, "  Max Buffer Size:      $(params.max_buffer_size)")
        println(io, "  LBT Threshold:        $(params.cs_threshold_dbm) dBm")
        println(io, "  Enable CS:            $(params.enable_carrier_sense)")
        println(io, "  Enable Capture:       $(params.enable_capture_effect)")
        println(io, "")
        println(io, "[Environment]")
        println(io, "  Area Size:            $(params.area_size_m) m")
        println(io, "  Path Loss Exponent:   $(params.pass_loss_exp)")
        println(io, "  Shadowing Std Dev:    $(params.shadowing_std_db) dB")
        println(io, "  TX Power:             $(params.tx_power_dbm) dBm")
        println(io, "  Noise Figure:         $(params.noise_figure_db) dB")
        println(io, "")
        println(io, "[Simulation Control]")
        println(io, "  Duration:             $(params.simulation_duration_ms) ms ($(params.simulation_duration_ms/60000.0) min)")
        println(io, "  Max Startup Delay:    $(params.max_startup_delay_ms) ms")
        println(io, "  Duty Cycle:           $(params.duty_cycle * 100.0) %")
        println(io, "  Mean Event Interval:  $(params.mean_event_interval_ms) ms")
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
            params.cs_threshold_dbm,
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
            false, false, false, false, MAX_BUFFER_SIZE,  # ← ログ無効, force_async=false, buffer上書き
            ENABLE_CARRIER_SENSE, ENABLE_CAPTURE_EFFECT, TARGET_SYNC_RATE,
            false  # use_waveform_phy: 電力モデル（高速）
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
        #total_retransmissions = sum([r["total_retransmissions"] for r in results_list])
        #retransmission_rate = total_packets > 0 ? (total_retransmissions / total_packets) : 0.0
        
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
        
        # 送信パケット数の平均
        total_generated_list = [r["total_generated"] for r in results_list]
        total_packets_list = [r["total_packets"] for r in results_list]
        mean_total_generated = mean(total_generated_list)
        mean_total_packets = mean(total_packets_list)

        # サマリーデータに追加（簡略版）
        push!(summary_data, (
            num_terminals = num_terminals,
            num_trials = NUM_TRIALS,
            mean_per = mean_per,
            std_per = std_per,
            mean_original_per = mean_original_per,
            std_original_per = std_original_per,
            mean_sync = mean_sync,
            mean_synced_per = mean_synced_per,
            std_synced_per = std_synced_per,
            mean_async_per = mean_async_per,
            std_async_per = std_async_per,
            mean_throughput_bps = mean_throughput_bps,
            std_throughput_bps = std_throughput_bps,
            mean_norm_throughput = mean_norm_throughput,
            std_norm_throughput = std_norm_throughput,
            total_buffer_drops = total_buffer_drops,
            total_generated = total_generated,
            buffer_drop_rate = buffer_drop_rate,
            total_retransmissions = 0,
            retransmission_rate = 0.0,
            mean_total_generated = mean_total_generated,
            mean_total_packets = mean_total_packets
        ))
        
        # 結果表示
        println("\nResults for $num_terminals terminals:")
        @printf("  Mean PER (Total):    %.6f\n", mean_per)
        @printf("  Mean PER (Original): %.6f\n", mean_original_per)
        @printf("  Overall PER (Total): %.6f\n", overall_per)
        @printf("  Buffer Drop Rate:    %.6f (%d / %d)\n", buffer_drop_rate, total_buffer_drops, total_generated)
        @printf("  Throughput (bps):    %.2f (Mean) | %.2f (Overall)\n", mean_throughput_bps, overall_throughput_bps)
        @printf("  Throughput (Norm):   %.5f (Mean) | %.5f (Overall)\n", mean_norm_throughput, overall_norm_throughput)
        @printf("  Mean Sync Rate:      %.2f ± %.2f %%\n", mean_sync, std_sync)
        println("  ├─ Synced PER:       $(round(mean_synced_per, digits=6)) ± $(round(std_synced_per, digits=6))")
        println("  └─ Async PER:        $(round(mean_async_per, digits=6)) ± $(round(std_async_per, digits=6))")
    end
    
    # 全体のサマリーを表示
    println("\n" * "="^60)
    println("  SUMMARY: All Terminal Counts")
    println("  SF: $SPREADING_FACTOR, Slot: $(SLOT_LENGTH_MS)ms, Beacon: $(BEACON_INTERVAL_MS)ms")
    println("="^60)
    println("Terminals | Mean PER     | Original PER     | Thr (bps) | Sync Rate (%)")
    println("-"^75)
    for data in summary_data
        @printf("%9d | %12.4f | %16.4f | %9.1f | %13.2f\n", 
                data.num_terminals, data.mean_per, data.mean_original_per, data.mean_throughput_bps, data.mean_sync)
    end
    println("="^60)
    
    # サマリーCSVを保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    ack_str = ENABLE_ACK ? "ACK" : "NoACK"
    summary_file = joinpath(output_dir, "summary_SF$(SPREADING_FACTOR)_slot$(Int(SLOT_LENGTH_MS))ms_beacon$(Int(BEACON_INTERVAL_MS))ms_$(ack_str)_$(NUM_TRIALS)trials_$(timestamp).csv")
    open(summary_file, "w") do io
        cols = ["num_terminals", "num_trials", 
                "mean_per", "std_per", 
                "mean_original_per", "std_original_per",
                "mean_sync",
                "mean_synced_per", "std_synced_per",
                "mean_async_per", "std_async_per",
                "mean_throughput_bps", "std_throughput_bps",
                "mean_norm_throughput", "std_norm_throughput",
                "total_buffer_drops", "total_generated", "buffer_drop_rate",
                "total_retransmissions", "retransmission_rate",
                "mean_total_generated", "mean_total_packets"]
        println(io, join(cols, ","))
        for data in summary_data
            @printf(io, "%d,%d,%.6f,%.6f,%.6f,%.6f,%.2f,%.6f,%.6f,%.6f,%.6f,%.2f,%.2f,%.6f,%.6f,%d,%d,%.6f,%d,%.6f,%.2f,%.2f\n",
                    data.num_terminals, data.num_trials,
                    data.mean_per, data.std_per,
                    data.mean_original_per, data.std_original_per,
                    data.mean_sync,
                    data.mean_synced_per, data.std_synced_per,
                    data.mean_async_per, data.std_async_per,
                    data.mean_throughput_bps, data.std_throughput_bps,
                    data.mean_norm_throughput, data.std_norm_throughput,
                    data.total_buffer_drops, data.total_generated, data.buffer_drop_rate,
                    data.total_retransmissions, data.retransmission_rate,
                    data.mean_total_generated, data.mean_total_packets)
        end
    end
    println("\nSummary saved to: $summary_file")
    println("="^60)
end