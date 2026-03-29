
using Distributed
using Dates
using Printf
using Statistics
using CSV
using DataFrames

# 並列ワーカーの追加
if nprocs() == 1
    addprocs(4) # 4並列
end

@everywhere using Statistics
@everywhere using Random
@everywhere using Dates
@everywhere using Printf

# プロジェクトのルートディレクトリ
@everywhere begin
    # カレントディレクトリをプロジェクトルートとして扱う
    # 必要であれば cd("/path/to/project") を実行するが、
    # 基本的には実行時のディレクトリを使用する
    include("Prop.jl")
end

@everywhere module AsyncModule
    using Statistics
    using Random
    using Dates
    using Printf
    using DataFrames
    using CSV
    
    include("Prop_Async.jl")
end

# シミュレーション設定
@everywhere const NUM_TERMINALS = 100
@everywhere const NUM_TRIALS = 10
@everywhere const AREA_SIZE = 2000.0

@everywhere function calculate_pair_distance_distribution(terminals, distance_bins)
    num_bins = length(distance_bins) - 1  # 0:500:2000 は5要素(0,500,1000,1500,2000) -> 4ビン
    step_size = distance_bins[2] - distance_bins[1]
    pair_counts = zeros(Int, num_bins)
    
    n = length(terminals)
    for i in 1:n
        for j in (i+1):n
            t1 = terminals[i]
            t2 = terminals[j]
            dist = sqrt((t1.x_m - t2.x_m)^2 + (t1.y_m - t2.y_m)^2)
            
            # ビン計算
            if dist >= distance_bins[1] && dist < distance_bins[end]
                bin_idx = Int(floor((dist - distance_bins[1]) / step_size)) + 1
                if bin_idx <= num_bins
                    pair_counts[bin_idx] += 1
                end
            end
        end
    end
    
    return pair_counts
end

# ... (omitted code) ...

function aggregate_risk(results, distance_bins)
    num_bins = length(distance_bins) - 1
    step_size = distance_bins[2] - distance_bins[1]
    
    total_collisions = zeros(Int, num_bins)
    total_pairs = zeros(Int, num_bins)
    
    for res in results
        # ペア数加算
        total_pairs .+= res["pair_counts"]
        
        # 衝突数加算
        dists = res["collision_distances"]
        for d in dists
            if d >= distance_bins[1] && d < distance_bins[end]
                bin_idx = Int(floor((d - distance_bins[1]) / step_size)) + 1
                if bin_idx <= num_bins
                    total_collisions[bin_idx] += 1
                end
            end
        end
    end
    
    # リスク算出 (Collisions / Pairs)
    risk = zeros(Float64, num_bins)
    for i in 1:num_bins
        if total_pairs[i] > 0
            risk[i] = total_collisions[i] / total_pairs[i]
        else
            risk[i] = 0.0
        end
    end
    
    return risk, total_collisions, total_pairs
end





@everywhere function generate_terminals_sync(params)
    ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
    dep_p = TerminalDeploymentParameters(
        "random_fixed", 
        0.0, 
        params.num_terminals, 
        params.area_size_m, 
        10.0, 
        params.area_size_m/2, 
        params.sync_center_freq_ghz*1e9, 
        params.pass_loss_exp, 
        1.0, 
        ref_pl_5g,
        params.sync_bs_x_m,
        params.sync_bs_y_m
    )
    return deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
end

@everywhere function generate_terminals_async(params)
    ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
    dep_p = AsyncModule.TerminalDeploymentParameters(
        "random_fixed", 
        0.0, 
        params.num_terminals, 
        params.area_size_m, 
        10.0, 
        params.area_size_m/2, 
        params.sync_center_freq_ghz*1e9, 
        params.pass_loss_exp, 
        1.0, 
        ref_pl_5g,
        params.sync_bs_x_m,
        params.sync_bs_y_m
    )
    return AsyncModule.deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
end

@everywhere function run_sync_simulation_risk(trial_idx, distance_bins)
    Random.seed!(1234 + trial_idx)
    params = create_integrated_params()
    params = IntegratedParameters(
        params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
        params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
        NUM_TERMINALS, AREA_SIZE, params.slot_length_ms,
        params.packet_airtime_ms, params.enable_carrier_sense, params.cs_threshold_dbm,
        params.spreading_factor, params.lora_payload_bytes, params.num_channels,
        params.beacon_interval_ms, params.simulation_duration_ms,
        params.max_startup_delay_ms, params.duty_cycle, params.mean_event_interval_ms,
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
        false, false, false, # Log off
        params.force_async_mode,
        params.max_buffer_size
    )
    
    # 1. ペア分布計算（同じシード）
    Random.seed!(1234 + trial_idx) # シードリセット
    terminals = generate_terminals_sync(params)
    pair_counts = calculate_pair_distance_distribution(terminals, distance_bins)
    
    # 2. シミュレーション実行（同じシード）
    Random.seed!(1234 + trial_idx) # シードリセット
    result = run_integrated_simulation_with_params(params)
    
    return Dict(
        "collision_distances" => result["collision_distances"],
        "pair_counts" => pair_counts
    )
end

@everywhere function run_async_simulation_risk(trial_idx, distance_bins)
    Random.seed!(5678 + trial_idx)
    params = AsyncModule.create_integrated_params_async()
    params = AsyncModule.IntegratedParameters(
        params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
        params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
        NUM_TERMINALS, AREA_SIZE, params.slot_length_ms,
        params.packet_airtime_ms, params.enable_carrier_sense, params.cs_threshold_dbm,
        params.spreading_factor, params.lora_payload_bytes, params.num_channels,
        params.beacon_interval_ms, params.simulation_duration_ms,
        params.max_startup_delay_ms, params.duty_cycle, params.mean_event_interval_ms,
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
        false, false, false,
        true,
        params.max_buffer_size
    )
    
    # 1. ペア分布計算
    Random.seed!(5678 + trial_idx)
    terminals = generate_terminals_async(params)
    pair_counts = calculate_pair_distance_distribution(terminals, distance_bins)
    
    # 2. 実行
    Random.seed!(5678 + trial_idx)
    result = AsyncModule.run_integrated_simulation_with_params(params)
    
    return Dict(
        "collision_distances" => result["collision_distances"],
        "pair_counts" => pair_counts
    )
end




function main()
    println("============================================================")
    println("  Normalized Collision Risk Analysis (Parallel)")
    println("  Terminals: $NUM_TERMINALS, Trials: $NUM_TRIALS")
    println("  Workers: $(nworkers())")
    println("============================================================")
    
    # 500m刻みに変更
    distance_bins = collect(0:500:2000)
    
    println("\n[1/2] Running Synchronized Method (Parallel)...")
    sync_results = pmap(i -> run_sync_simulation_risk(i, distance_bins), 1:NUM_TRIALS)
    sync_risk, sync_col, sync_pairs = aggregate_risk(sync_results, distance_bins)
    
    println("[2/2] Running Asynchronous Method (Parallel)...")
    async_results = pmap(i -> run_async_simulation_risk(i, distance_bins), 1:NUM_TRIALS)
    async_risk, async_col, async_pairs = aggregate_risk(async_results, distance_bins)
    
    # 出力
    output_dir = "collision_distance_analysis"
    mkpath(output_dir)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    csv_file = joinpath(output_dir, "collision_risk_N$(NUM_TERMINALS)_$(NUM_TRIALS)trials_$(timestamp).csv")
    
    open(csv_file, "w") do io
        println(io, "distance_bin_start,distance_bin_end,sync_risk_prob,async_risk_prob,sync_collisions,async_collisions,total_pairs")
        for i in 1:(length(distance_bins)-1)
            @printf(io, "%d,%d,%.6f,%.6f,%d,%d,%d\n",
                    distance_bins[i], distance_bins[i+1],
                    sync_risk[i], async_risk[i],
                    sync_col[i], async_col[i], sync_pairs[i])
        end
    end
    
    println("\n" * "="^60)
    println("  Analysis Complete!")
    println("  CSV saved to: $csv_file")
    println("="^60)
    
    println("\nNormalized Collision Risk Summary (Probability per Pair):")
    println("Distance (m)      | Sync Risk      | Async Risk")
    println("-"^60)
    for i in 1:(length(distance_bins)-1)
        mark = sync_pairs[i] < 50 ? "*" : "" # サンプル少
        
        @printf("%4d - %4d m | %.6f %s | %.6f %s\n",
                distance_bins[i], distance_bins[i+1],
                sync_risk[i], mark,
                async_risk[i], mark)
    end
    println("(*: Low pair count < 50)")
    println("="^60)
end

main()
