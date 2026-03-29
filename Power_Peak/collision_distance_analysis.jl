# 衝突距離分析用の新しい評価スクリプト
# collision_distance_analysis.jl

using Distributed
using Statistics, Printf, Dates, CSV, DataFrames

# 並列ワーカーの追加
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

@everywhere using Statistics
@everywhere include("Prop.jl")
@everywhere include("Prop_Async.jl")

"""
衝突距離分析スクリプト

提案手法（Sync）と比較手法（Async）の両方で、
衝突した端末ペアの距離分布を分析し、
隠れ端末問題を定量化する。
"""

if abspath(PROGRAM_FILE) == @__FILE__
    # ============================================================
    # 設定
    # ============================================================
    NUM_TERMINALS = 100  # テスト用に削減
    NUM_TRIALS = 10      # テスト用に削減
    SPREADING_FACTOR = 8
    SLOT_LENGTH_MS = 200.0
    BEACON_INTERVAL_MS = 100.0
    ENABLE_ACK = false
    MAX_RETRIES = 3
    MAX_BUFFER_SIZE = 10
    
    println("="^60)
    println("  Collision Distance Analysis")
    println("  Terminals: $NUM_TERMINALS, Trials: $NUM_TRIALS")
    println("="^60)
    
    # ============================================================
    # 提案手法（Sync）の実行
    # ============================================================
    println("\n[1/2] Running Synchronized Method...")
    
    params_sync = create_integrated_params()
    lora_params = create_lora_params(SPREADING_FACTOR, params_sync.lora_payload_bytes)
    airtime_ms = calculate_lora_airtime(lora_params)
    
    params_sync = IntegratedParameters(
        params_sync.signal_duration_us, params_sync.signal_bw_mhz, params_sync.terminal_bw_mhz,
        params_sync.tx_sampling_rate_mhz, params_sync.rx_sampling_rate_mhz, params_sync.tx_power_dbm, params_sync.noise_figure_db,
        NUM_TERMINALS,
        params_sync.area_size_m, SLOT_LENGTH_MS, airtime_ms,
        params_sync.enable_carrier_sense, params_sync.cs_threshold_dbm,
        SPREADING_FACTOR, params_sync.lora_payload_bytes, params_sync.num_channels,
        BEACON_INTERVAL_MS, params_sync.simulation_duration_ms,
        params_sync.max_startup_delay_ms,
        params_sync.duty_cycle, params_sync.mean_event_interval_ms,
        params_sync.shadowing_enabled, params_sync.shadowing_std_db, params_sync.pass_loss_exp,
        params_sync.sync_center_freq_ghz, params_sync.data_center_freq_ghz, params_sync.reference_path_loss_db,
        params_sync.sync_bs_x_m, params_sync.sync_bs_y_m,
        params_sync.sync_observation_duration_ms, params_sync.gw_tx_power_dbm,
        params_sync.noise_floor_window_ms, params_sync.detection_margin_db,
        params_sync.min_samples, params_sync.debounce_time_ms, params_sync.initial_window_duration_ms,
        params_sync.tx_jitter_max_ms,
        params_sync.enable_intermittent_rx, params_sync.intermittent_window_ms, params_sync.initial_search_duration_ms,
        params_sync.use_deterministic_jitter, params_sync.num_jitter_offsets, params_sync.deterministic_jitter_random_ms,
        params_sync.collision_model,
        ENABLE_ACK, MAX_RETRIES, params_sync.ack_timeout_ms,
        params_sync.rx1_delay_ms, params_sync.backoff_base_ms,
        params_sync.lbt_duration_ms, params_sync.lbt_sample_interval_ms,
        false, false, false, false, MAX_BUFFER_SIZE
    )
    
    results_sync = pmap(1:NUM_TRIALS) do i
        original_stdout = stdout
        redirect_stdout(devnull)
        result = nothing
        try
            result = run_integrated_simulation_with_params(params_sync)
        finally
            redirect_stdout(original_stdout)
        end
        return result
    end
    
    # ============================================================
    # 比較手法（Async）の実行
    # ============================================================
    println("[2/2] Running Asynchronous Method...")
    
    # Asyncも同じcreate_integrated_params()を使用し、force_async_modeで制御
    params_async = create_integrated_params()
    params_async = IntegratedParameters(
        params_async.signal_duration_us, params_async.signal_bw_mhz, params_async.terminal_bw_mhz,
        params_async.tx_sampling_rate_mhz, params_async.rx_sampling_rate_mhz, params_async.tx_power_dbm, params_async.noise_figure_db,
        NUM_TERMINALS,
        params_async.area_size_m, SLOT_LENGTH_MS, airtime_ms,
        params_async.enable_carrier_sense, params_async.cs_threshold_dbm,
        SPREADING_FACTOR, params_async.lora_payload_bytes, params_async.num_channels,
        BEACON_INTERVAL_MS, params_async.simulation_duration_ms,
        params_async.max_startup_delay_ms,
        params_async.duty_cycle, params_async.mean_event_interval_ms,
        params_async.shadowing_enabled, params_async.shadowing_std_db, params_async.pass_loss_exp,
        params_async.sync_center_freq_ghz, params_async.data_center_freq_ghz, params_async.reference_path_loss_db,
        params_async.sync_bs_x_m, params_async.sync_bs_y_m,
        params_async.sync_observation_duration_ms, params_async.gw_tx_power_dbm,
        params_async.noise_floor_window_ms, params_async.detection_margin_db,
        params_async.min_samples, params_async.debounce_time_ms, params_async.initial_window_duration_ms,
        params_async.tx_jitter_max_ms,
        params_async.enable_intermittent_rx, params_async.intermittent_window_ms, params_async.initial_search_duration_ms,
        params_async.use_deterministic_jitter, params_async.num_jitter_offsets, params_async.deterministic_jitter_random_ms,
        params_async.collision_model,
        ENABLE_ACK, MAX_RETRIES, params_async.ack_timeout_ms,
        params_async.rx1_delay_ms, params_async.backoff_base_ms,
        params_async.lbt_duration_ms, params_async.lbt_sample_interval_ms,
        true, false, false, true, MAX_BUFFER_SIZE  # force_async_mode = true
    )
    
    results_async = pmap(1:NUM_TRIALS) do i
        original_stdout = stdout
        redirect_stdout(devnull)
        result = nothing
        try
            result = run_integrated_simulation_with_params(params_async)
        finally
            redirect_stdout(original_stdout)
        end
        return result
    end
    
    # ============================================================
    # 衝突距離データの集計
    # ============================================================
    println("\nAggregating collision distance data...")
    
    # 距離ビンの定義
    distance_bins = collect(0:100:2000)
    num_bins = length(distance_bins)
    
    # Sync の集計
    sync_collision_matrix = zeros(Int, NUM_TRIALS, num_bins)
    for (trial_idx, result) in enumerate(results_sync)
        counts = result["collision_counts_per_bin"]
        sync_collision_matrix[trial_idx, :] = counts
    end
    
    # 正規化: 各試行ごとに正規化してから平均を取る
    sync_rate_matrix = zeros(Float64, NUM_TRIALS, num_bins)
    for trial in 1:NUM_TRIALS
        trial_total = sum(sync_collision_matrix[trial, :])
        if trial_total > 0
            sync_rate_matrix[trial, :] = (sync_collision_matrix[trial, :] ./ trial_total) .* 100.0
        end
    end
    
    sync_collision_rate = vec(mean(sync_rate_matrix, dims=1))
    sync_std_rate = vec(std(sync_rate_matrix, dims=1))
    
    # Async の集計
    async_collision_matrix = zeros(Int, NUM_TRIALS, num_bins)
    for (trial_idx, result) in enumerate(results_async)
        counts = result["collision_counts_per_bin"]
        async_collision_matrix[trial_idx, :] = counts
    end
    
    # 正規化: 各試行ごとに正規化してから平均を取る
    async_rate_matrix = zeros(Float64, NUM_TRIALS, num_bins)
    for trial in 1:NUM_TRIALS
        trial_total = sum(async_collision_matrix[trial, :])
        if trial_total > 0
            async_rate_matrix[trial, :] = (async_collision_matrix[trial, :] ./ trial_total) .* 100.0
        end
    end
    
    async_collision_rate = vec(mean(async_rate_matrix, dims=1))
    async_std_rate = vec(std(async_rate_matrix, dims=1))
    
    # ============================================================
    # CSV出力
    # ============================================================
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    output_dir = "collision_distance_analysis"
    mkpath(output_dir)
    
    csv_file = joinpath(output_dir, "collision_distance_N$(NUM_TERMINALS)_$(NUM_TRIALS)trials_$(timestamp).csv")
    
    open(csv_file, "w") do io
        println(io, "distance_bin_center,sync_collision_rate,sync_std_rate,async_collision_rate,async_std_rate")
        for i in 1:num_bins
            bin_center = distance_bins[i] + 100  # ビンの中心値
            @printf(io, "%d,%.2f,%.2f,%.2f,%.2f\n",
                    bin_center,
                    sync_collision_rate[i], sync_std_rate[i],
                    async_collision_rate[i], async_std_rate[i])
        end
    end
    
    println("\n" * "="^60)
    println("  Analysis Complete!")
    println("  CSV saved to: $csv_file")
    println("="^60)
    
    # 結果のサマリー表示
    println("\nCollision Distance Summary (Normalized %):")
    println("Distance (m) | Sync (rate±std) | Async (rate±std)")
    println("-"^60)
    for i in 1:num_bins
        bin_center = distance_bins[i] + 100
        @printf("%d-%d m    | %.2f%% ± %.2f%%   | %.2f%% ± %.2f%%\n",
                distance_bins[i], distance_bins[i]+200,
                sync_collision_rate[i], sync_std_rate[i],
                async_collision_rate[i], async_std_rate[i])
    end
    println("="^60)
end
