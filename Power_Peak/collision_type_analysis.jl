# using Distributed
using Dates
using Printf
using Statistics
using CSV
using DataFrames

# 並列処理のセットアップ
# 既存のワーカーがあれば削除
# if nprocs() > 1
#     rmprocs(workers())
# end
# 4プロセス追加
# addprocs(4)

# 必要なモジュールを全プロセスにロード
# @everywhere begin
    using Statistics
    using Random
    using Dates
    using Printf
    
    # プロジェクトのルートディレクトリ
    PROJECT_ROOT = "/Users/kimparaibuki/Cursor/Power_Peak"
    cd(PROJECT_ROOT)
    
    include("Prop.jl")
    
    module AsyncModule
        using Statistics
        using Random
        using Dates
        using Printf
        using DataFrames
        using CSV
        
        # 必要な標準ライブラリを再エクスポートまたは使用
        
        include("Prop_Async.jl")
    end
# end

# シミュレーション設定
const NUM_TERMINALS = 100
const NUM_TRIALS = 10  # 試行回数
const AREA_SIZE = 2000.0

# @everywhere begin
    function run_sync_simulation_for_collision_type(trial_idx)
        # 固定シード (再現性のため)
        Random.seed!(1234 + trial_idx)
        
        # パラメータ設定
        params = create_integrated_params()
        
        # 評価用にパラメータ上書き
        # 端末数100, エリア2000m
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
            false, false, false, # ファイル出力などは無効化
            params.force_async_mode,
            params.max_buffer_size
        )
        
        # 実行
        result = run_integrated_simulation_with_params(params)
        
        # 必要なデータのみ返す
        return Dict(
            "collision_distances" => result["collision_distances"],
            "collision_types" => result["collision_types"]
        )
    end

    function run_async_simulation_for_collision_type(trial_idx)
        # 固定シード
        Random.seed!(5678 + trial_idx)
        
        # パラメータ設定 (Async)
        params = AsyncModule.create_integrated_params_async() # Prop_Async.jl の関数
        
        # 評価用にパラメータ上書き (Syncと同じ条件)
        # AsyncModule.IntegratedParameters を使う必要がある
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
            false, false, false, # ファイル出力などは無効化
            true, # force_async_mode = true (念のため)
            params.max_buffer_size
        )
        
        # 実行 (AsyncModule内の関数)
        # 注意: run_integrated_simulation_async も AsyncModule 内にある
        # ただしその関数名は Prop_Async.jl で定義されているもの ('run_integrated_simulation' か 'run_integrated_simulation_async' か確認必要)
        # Prop_Async.jl では 'run_integrated_simulation_with_integrated_params' 的な名前ではなく、'run_integrated_simulation_async' ?
        # 前のステップで create_integrated_params_async に変更したが、実行関数名は変更していないはず。
        # Prop.jl と同じ 'run_integrated_simulation_with_params' の可能性がある。
        
        result = AsyncModule.run_integrated_simulation_with_params(params)
        
        return Dict(
            "collision_distances" => result["collision_distances"],
            "collision_types" => result["collision_types"]
        )
    end
# end

function aggregate_collision_types(results, distance_bins)
    num_bins = length(distance_bins)
    
    # [Trial][Bin] -> Count
    hidden_counts = zeros(Int, NUM_TRIALS, num_bins)
    exposed_counts = zeros(Int, NUM_TRIALS, num_bins)
    total_counts = zeros(Int, NUM_TRIALS, num_bins)
    
    for (t_idx, res) in enumerate(results)
        dists = res["collision_distances"]
        types = res["collision_types"]
        
        for i in 1:length(dists)
            d = dists[i]
            t = types[i]
            
            # ビン計算 (100m刻み)
            bin_idx = min(Int(floor(d / 100)) + 1, num_bins)
            
            total_counts[t_idx, bin_idx] += 1
            if t == "Hidden"
                hidden_counts[t_idx, bin_idx] += 1
            else
                exposed_counts[t_idx, bin_idx] += 1
            end
        end
    end
    
    # 平均と標準偏差を計算
    avg_hidden_rate = zeros(Float64, num_bins)
    std_hidden_rate = zeros(Float64, num_bins)
    
    for bin in 1:num_bins
        rates = Float64[]
        for t_idx in 1:NUM_TRIALS
            total = total_counts[t_idx, bin]
            if total > 0
                push!(rates, (hidden_counts[t_idx, bin] / total) * 100.0)
            else
                # push!(rates, 0.0) # データなしは除外するか0にするか... データなしは除外が安全
            end
        end
        
        if !isempty(rates)
            avg_hidden_rate[bin] = mean(rates)
            std_hidden_rate[bin] = std(rates)
        else
            avg_hidden_rate[bin] = 0.0
            std_hidden_rate[bin] = 0.0
        end
    end
    
    return avg_hidden_rate, std_hidden_rate, sum(total_counts, dims=1)[1, :]
end

function main()
    println("============================================================")
    println("  Collision Type Analysis (Hidden vs Exposed)")
    println("  Terminals: $NUM_TERMINALS, Trials: $NUM_TRIALS")
    println("============================================================")
    
    # 距離ビン (0-2000m, 100m刻み)
    distance_bins = collect(0:100:2000)
    
    # --- Sync ---
    println("\n[1/2] Running Synchronized Method...")
    sync_results = []
    for i in 1:NUM_TRIALS
        push!(sync_results, run_sync_simulation_for_collision_type(i))
    end
    sync_hidden_avg, sync_hidden_std, sync_total_counts = aggregate_collision_types(sync_results, distance_bins)
    
    # --- Async ---
    println("[2/2] Running Asynchronous Method...")
    async_results = []
    for i in 1:NUM_TRIALS
        push!(async_results, run_async_simulation_for_collision_type(i))
    end
    async_hidden_avg, async_hidden_std, async_total_counts = aggregate_collision_types(async_results, distance_bins)
    
    # --- 結果出力 ---
    output_dir = "collision_distance_analysis"
    mkpath(output_dir)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    csv_file = joinpath(output_dir, "collision_type_N$(NUM_TERMINALS)_$(NUM_TRIALS)trials_$(timestamp).csv")
    
    open(csv_file, "w") do io
        println(io, "distance_bin_start,sync_hidden_rate,sync_hidden_std,sync_total_collisions,async_hidden_rate,async_hidden_std,async_total_collisions")
        for i in 1:length(distance_bins)
            @printf(io, "%d,%.2f,%.2f,%d,%.2f,%.2f,%d\n",
                    distance_bins[i],
                    sync_hidden_avg[i], sync_hidden_std[i], sync_total_counts[i],
                    async_hidden_avg[i], async_hidden_std[i], async_total_counts[i])
        end
    end
    
    println("\n" * "="^60)
    println("  Analysis Complete!")
    println("  CSV saved to: $csv_file")
    println("="^60)
    
    println("\nHidden Terminal Rate Summary:")
    println("Distance (m) | Sync Hidden Rate | Async Hidden Rate")
    println("-"^60)
    for i in 1:length(distance_bins)
        # データが少なすぎるビンは表示スキップしてもよいが全て出す
        count_sync = sync_total_counts[i]
        count_async = async_total_counts[i]
        
        # 信頼性目安: サンプル数10以上
        mark_sync = count_sync < 10 ? "*" : ""
        mark_async = count_async < 10 ? "*" : ""
        
        @printf("%4d - %4d m | %6.2f%% ± %5.2f%% %s | %6.2f%% ± %5.2f%% %s\n",
                distance_bins[i], distance_bins[i]+100,
                sync_hidden_avg[i], sync_hidden_std[i], mark_sync,
                async_hidden_avg[i], async_hidden_std[i], mark_async)
    end
    println("(*: Low sample size < 10)")
    println("="^60)
end

main()
