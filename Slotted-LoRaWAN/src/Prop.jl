using Random, Statistics, Printf, DataFrames, CSV, Plots, Dates, LinearAlgebra, DSP

# ==========================================
include(joinpath(@__DIR__, "SlottedLoRaWAN.jl"))
using .SlottedLoRaWAN

using .PacketGeneration


# ==========================================
# Protocol Scenario Settings
# ==========================================
# (Removed local const, now using params.enable_carrier_sense)

"""
    run_integrated_simulation_with_params(params::IntegratedParameters)

パラメータを受け取ってシミュレーションを実行する（評価スクリプト用）
"""
function run_integrated_simulation_with_params(params::IntegratedParameters)
    println("="^60)
    println("   Integrated LoRa Simulation")
    println("   (Hi-Fi PHY Sync -> Time-Sorted MAC)")
    println("="^60)
    
    # ★ ToAが未計算の場合は計算 ★
    if params.packet_airtime_ms == 0.0
        lora_params = create_lora_params(params.spreading_factor, params.lora_payload_bytes)
        params.packet_airtime_ms = calculate_lora_airtime(lora_params)
    end
    
    # ★ スロット長が未設定の場合 (0.0) は ToA + 100ms に設定 ★
    if params.slot_length_ms == 0.0
        params.slot_length_ms = params.packet_airtime_ms + 100.0
    end
    
    println("\nLoRa 設定:")
    println("  SF: $(params.spreading_factor)")
    println("  ペイロード: $(params.lora_payload_bytes) bytes")
    println("  計算された ToA: $(round(params.packet_airtime_ms, digits=2)) ms")
    println("  スロット長: $(params.slot_length_ms) ms")
    println()
    
    # 1. 端末配置 (5G帯でのパスロス計算)
    # 5G帯の基準パスロスを計算
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
        params.sync_bs_x_m,  # 同期BS X座標
        params.sync_bs_y_m   # 同期BS Y座標
    )
    terminals = deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
    
    # ★ integrated_lora_sim.jl と合わせるため、クロックドリフトを0に強制上書き ★
    #terminals = [TerminalInfo(t.terminal_id, t.x_m, t.y_m, t.distance_m, t.path_loss_db, t.shadowing_db, t.total_loss_db, t.rx_power_dbm, 0.0, 1.0) for t in terminals]
    

    
    # 端末情報表示（5G帯のみ） -> 非表示
    # println("\n端末配置（5G同期信号）:")
    # for t in terminals
    #     println("• 端末$(t.terminal_id): 位置($(round(t.x_m, digits=1)), $(round(t.y_m, digits=1))) m, 距離$(round(t.distance_m, digits=1)) m")
    #     println("  - パスロス: $(round(t.path_loss_db, digits=2)) dB")
    #     println("  - シャドウイング: $(round(t.shadowing_db, digits=1)) dB")
    #     println("  - 総損失: $(round(t.total_loss_db, digits=2)) dB")
    #     println("  - 受信電力: $(round(t.rx_power_dbm, digits=1)) dBm")
    # end
    
    # 2. Phase 1: 高精度同期を実行
    # ここで main_simulation.jl のロジックが走り、各端末の開始時刻が決まる
    # 実行場所に関わらず、スクリプトのあるフォルダ(Power_Peak)の下に保存する
    output_dir = joinpath(@__DIR__, "result_Prop")
    mkpath(output_dir)
    # 3. Phase 1: High-Fidelity Synchronization (On-Demand)
    println("Phase 1: High-Fidelity Synchronization per Terminal (On-Demand)...")
    
    # 信号生成用パラメータ
    sig_gen_params = SignalParameters(
        params.signal_duration_us*1e-6,
        params.sync_center_freq_ghz*1e9,
        params.signal_bw_mhz*1e6,
        params.tx_sampling_rate_mhz*1e6,
        params.gw_tx_power_dbm
    )
    
    # 送信スケジュール計算 (軽量)
    sync_start_delay_ms = 10.0 + rand() * 20.0
    beacon_times = calculate_beacon_times(params.beacon_interval_ms, params.sync_observation_duration_ms, sync_start_delay_ms)
    println("Sync Signal Schedule: $(length(beacon_times)) beacons in $(params.sync_observation_duration_ms) ms")
    
    # 同期実行
    terminal_sync_infos, sync_log_df = perform_hifi_synchronization(terminals, beacon_times, sig_gen_params, params, output_dir)
    
    # 3. Phase 2: 初期パケット発生イベントのスケジュール
    println("Phase 2: Initializing Queues & Scheduling Packet Generation...")
    
    event_queue = MACEvent[]
    next_packet_id = 1
    
    # 端末ごとの状態管理
    terminal_queues = Dict{Int, Vector{Int}}(t.terminal_id => Int[] for t in terminals)
    terminal_is_transmitting = Dict{Int, Bool}(t.terminal_id => false for t in terminals)
    next_available_time = Dict{Int, Float64}(t.terminal_id => 0.0 for t in terminals)
    
    for t in terminals
        # 最初のパケット発生時刻を決定
        # 起動時刻直後またはランダムな遅延後
        startup_ms = rand() * params.max_startup_delay_ms
        initial_interval = -params.mean_event_interval_ms * log(rand()) * t.clock_drift_factor
        first_gen_time = startup_ms + initial_interval
        
        if first_gen_time < params.simulation_duration_ms
            new_event = MACEvent(first_gen_time, PACKET_GEN, t.terminal_id, 0, nothing)
            push!(event_queue, new_event)
        end
    end
    
    # 時刻順にソート (降順: 末尾からpopするため)
    sort!(event_queue, by = x -> x.time_ms, rev=true)
    
    # 同期成功端末のIDセットを作成
    synced_terminal_ids = Set([tid for (tid, start_ms) in terminal_sync_infos if start_ms !== nothing])
    
    # イベント挿入用ヘルパー
    function push_event!(q, ev)
        idx = searchsortedfirst(q, ev, by=x->x.time_ms, rev=true)
        insert!(q, idx, ev)
    end

    # 4. Phase 3: MAC層実行 (Discrete Event Simulation)
    println("Phase 3: Running MAC Layer (Queue-based with Independent Generation)...")
    
    # ノイズパラメータ
    noise_bw_hz = params.terminal_bw_mhz * 1e6
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    
    # CS/LBT メトリクス
    cs_attempts = 0
    cs_blocked = 0
    cs_blocked_distances = Float64[]
    
    # SINR/SNR 閾値 (拡散率SFに依存)
    required_sir_db = 6.0
    required_snr_db = -15.0 # SF10 default
    if params.spreading_factor == 7 required_snr_db = -7.5
    elseif params.spreading_factor == 8 required_snr_db = -10.0
    elseif params.spreading_factor == 9 required_snr_db = -12.5
    elseif params.spreading_factor == 11 required_snr_db = -17.5
    elseif params.spreading_factor == 12 required_snr_db = -20.0
    end
    
    active_tx = TransmissionRecord[]      # 現在空中に存在しているパケット
    all_tx_records = TransmissionRecord[]  # 履歴（すべての試行）
    packet_success_map = Dict{Int, Bool}()
    
    # 統計用
    buffer_drops = 0
    total_generated_packets = 0
    
    # 衝突距離分析用
    collision_distances = Float64[]
    collision_types = String[]  # "Hidden" or "Exposed"
    
    # メインループ
    while !isempty(event_queue)
        event = pop!(event_queue)
        curr_t = event.time_ms
        me_id = event.terminal_id
        me = terminals[me_id]
        
        if event.type == PACKET_GEN
            # --- パケット生起イベント ---
            total_generated_packets += 1
            pid = next_packet_id
            next_packet_id += 1
            
            # バッファに空きがあるか確認
            if length(terminal_queues[me_id]) < params.max_buffer_size
                push!(terminal_queues[me_id], pid)
                
                # 端末が送信中でなければ送信開始
                if !terminal_is_transmitting[me_id]
                    terminal_is_transmitting[me_id] = true
                    # Duty Cycleを守りつつ最短時刻で送信
                    base_tx_time = max(curr_t, next_available_time[me_id])
                    schedule_tx_start!(event_queue, me, pid, base_tx_time, params, terminal_sync_infos, synced_terminal_ids)
                end
            else
                # バッファフルによる破棄
                buffer_drops += 1
                # 破棄記録 (解析用ダミーレコード)
                dummy_rec = TransmissionRecord(me_id, curr_t, curr_t, 0.0, me.x_m, me.y_m, "BufferDrop", -Inf, 0, 0, false, pid, "BufferFull", -Inf, -Inf, 0, false, "", "")
                push!(all_tx_records, dummy_rec)
            end
            
            # 次のパケット発生をスケジュール
            next_gen_interval = -params.mean_event_interval_ms * log(rand()) * me.clock_drift_factor
            next_gen_time = curr_t + next_gen_interval
            if next_gen_time < params.simulation_duration_ms
                new_gen_event = MACEvent(next_gen_time, PACKET_GEN, me_id, 0, nothing)
                idx = searchsortedfirst(event_queue, new_gen_event, by=x->x.time_ms, rev=true)
                insert!(event_queue, idx, new_gen_event)
            end
            
        elseif event.type == TX_START
            # --- 送信開始イベント ---
            cand = event.data
            pid = event.packet_id
            
            # Duty Cycle チェック (予備)
            if curr_t < next_available_time[me_id]
                # まだ待つ必要がある場合はスケジュールし直し
                schedule_tx_start!(event_queue, me, pid, next_available_time[me_id], params, terminal_sync_infos, synced_terminal_ids)
                continue
            end
            
            # キャリアセンス (LBT)
            is_busy = false
            if params.enable_carrier_sense
                cs_attempts += 1
                # 端末間RSSベースの正確なLBT
                for other in active_tx
                    if other.channel == cand.channel
                        other_term = terminals[other.terminal_id]
                        # 端末間の距離を計算
                        dx = me.x_m - other_term.x_m
                        dy = me.y_m - other_term.y_m
                        dist_peer = sqrt(dx^2 + dy^2)
                        
                        # 物理チャネルモデルを適用（距離減衰 + シャドウイング）
                        pl_peer = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist_peer, 1.0))
                        sh_peer = randn() * params.shadowing_std_db # LBT判定用の独立的シャドウイング
                        p_peer_dbm = params.tx_power_dbm - (pl_peer + sh_peer)
                        
                        if p_peer_dbm >= params.cs_threshold_dbm
                            is_busy = true
                            cs_blocked += 1
                            push!(cs_blocked_distances, dist_peer)
                            break
                        end
                    end
                end
            end
            
            if !is_busy
                # 送信実行
                # 送信実行
                
                dur = params.packet_airtime_ms
                dist_gw = sqrt(me.x_m^2 + me.y_m^2)
                pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist_gw, 1.0))
                rx_gw = params.tx_power_dbm - (pl_gw + me.shadowing_db)

                rec = TransmissionRecord(
                    me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, 
                    me.x_m, me.y_m, "Processing", rx_gw, cand.channel,
                    0, false, pid,
                    "None", 0.0, 0.0, 0, false, "", ""
                )
                
                push!(active_tx, rec)
                push!(all_tx_records, rec)
                
                # 送信終了イベント
                push_event!(event_queue, MACEvent(curr_t + dur, TX_END, me_id, pid, rec))
                
                # 次回送信可能時刻 (Duty Cycle)
                off_period = dur * (1.0 / params.duty_cycle - 1.0)
                next_available_time[me_id] = curr_t + dur + off_period
            else
                # Busyならバックオフしてリトライ
                backoff = me_id in synced_terminal_ids ? (rand(1:3) * params.slot_length_ms) : (params.lbt_duration_ms + rand()*100.0)
                schedule_tx_start!(event_queue, me, pid, curr_t + backoff, params, terminal_sync_infos, synced_terminal_ids)
            end
            
        elseif event.type == TX_END
            # --- 送信終了イベント ---
            rec = event.data
            pid = event.packet_id
            filter!(r -> r !== rec, active_tx)

            # --- 成功判定 (SINR/SNR) ---
            interferers = filter(other -> (other !== rec && other.channel == rec.channel && 
                                          max(rec.start_ms, other.start_ms) < min(rec.end_ms, other.end_ms)), all_tx_records)
            interference_w = sum([10^(inf.rx_power_at_gw/10)*1e-3 for inf in interferers], init=0.0)
            signal_w = 10^(rec.rx_power_at_gw/10)*1e-3
            noise_w = 10^(noise_power_dbm/10)*1e-3
            
            sinr_db = interference_w > 0 ? 10*log10(signal_w / interference_w) : Inf
            snr_db = 10*log10(signal_w / noise_w)
            
            rec.sinr_db = sinr_db
            rec.snr_db = snr_db
            rec.num_interferers = length(interferers)

            if params.enable_capture_effect
                # キャプチャ効果あり (SINR/SNRベース)
                is_successful = (sinr_db >= required_sir_db) && (snr_db >= required_snr_db)
            else
                # キャプチャ効果なし (衝突があれば失敗)
                is_successful = isempty(interferers)
            end
            
            if is_successful
                rec.status = "Success"
                packet_success_map[pid] = true
            else
                rec.status = "Collision"
                # 衝突分析
                for inf in interferers
                    inf_term_id = inf.terminal_id
                    if inf_term_id != me_id
                        dx = terminals[me_id].x_m - terminals[inf_term_id].x_m
                        dy = terminals[me_id].y_m - terminals[inf_term_id].y_m
                        distance = sqrt(dx^2 + dy^2)
                        push!(collision_distances, distance)
                        
                        pl_p = PathLossParameters(distance, params.data_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, params.reference_path_loss_db)
                        pl_db = calculate_path_loss(pl_p)
                        rx_power_from_inf = params.tx_power_dbm - (pl_db + terminals[inf_term_id].shadowing_db)
                        
                        if rx_power_from_inf < params.cs_threshold_dbm
                            push!(collision_types, "Hidden")
                        else
                            push!(collision_types, "Exposed")
                        end
                    end
                end
            end

            # パケットをキューから削除 (再送なし)
            if !isempty(terminal_queues[me_id]) && terminal_queues[me_id][1] == pid
                popfirst!(terminal_queues[me_id])
            end
            
            # 次のパケットがあれば送信開始
            if !isempty(terminal_queues[me_id])
                next_pid = terminal_queues[me_id][1]
                schedule_tx_start!(event_queue, me, next_pid, max(curr_t + 1.0, next_available_time[me_id]), params, terminal_sync_infos, synced_terminal_ids)
            else
                terminal_is_transmitting[me_id] = false
            end
        end
    end
    
    # 5. Summary & Results
    finished_tx = all_tx_records
    total_terminals = params.num_terminals
    
    # 成功パケットを取得 (重複IDを排除)
    success_pid_set = Set{Int}()
    for r in finished_tx
        if r.status == "Success"
            push!(success_pid_set, r.packet_id)
        end
    end
    
    # 全生成パケット数 (PACKET_GENでカウント済み)
    original_total = total_generated_packets
    original_success = length(success_pid_set)
    original_per = original_total > 0 ? (1.0 - original_success / max(original_total, 1)) : 0.0
    
    # 全試行ベースの統計 (BufferDropを除く)
    attempt_tx = filter(r -> r.status != "BufferDrop", finished_tx)
    total_attempts = length(attempt_tx)
    total_success_attempts = count(r -> r.status == "Success", attempt_tx)
    total_failed_attempts = total_attempts - total_success_attempts
    attempt_per = total_attempts > 0 ? (total_failed_attempts / max(total_attempts, 1)) : 0.0
    
    # スループット計算
    sim_duration_s = params.simulation_duration_ms / 1000.0
    total_data_bits = original_success * params.lora_payload_bytes * 8
    throughput_bps = total_data_bits / sim_duration_s
    
    # 正規化スループット
    norm_throughput = (original_success * params.packet_airtime_ms) / (params.simulation_duration_ms * params.num_channels)
    
    # 同期/非同期別のPER計算
    synced_attempts = 0
    synced_failures = 0
    async_attempts = 0
    async_failures = 0
    
    for r in attempt_tx
        if r.terminal_id in synced_terminal_ids
            synced_attempts += 1
            if r.status != "Success"
                synced_failures += 1
            end
        else
            async_attempts += 1
            if r.status != "Success"
                async_failures += 1
            end
        end
    end
    
    synced_per = synced_attempts > 0 ? (synced_failures / synced_attempts) : 0.0
    async_per = async_attempts > 0 ? (async_failures / async_attempts) : 0.0
    
    # 衝突距離のビニング
    distance_bins = collect(0:100:2000)
    collision_counts_per_bin = zeros(Int, length(distance_bins))
    
    for dist in collision_distances
        bin_idx = min(Int(floor(dist / 100)) + 1, length(distance_bins))
        collision_counts_per_bin[bin_idx] += 1
    end
    
    println("\nSimulation Summary:")
    println("  Total Generated:     $original_total")
    println("  Dropped (Buffer):    $buffer_drops")
    println("  Original Success:    $original_success")
    println("  Original PER:        $(round(original_per, digits=6))")
    println("  Total Attempts:      $total_attempts")
    println("  Total Collisions:    $total_failed_attempts")
    println("  Attempt PER:         $(round(attempt_per, digits=6))")
    println("  Throughput (bps):    $(round(throughput_bps, digits=2))")
    println("  Throughput (Norm):   $(round(norm_throughput, digits=5))")
    
    println("\nCarrier Sense (LBT) Summary:")
    println("  CS Enabled:         $(params.enable_carrier_sense)")
    println("  Total CS Attempts:  $cs_attempts")
    println("  Busy (Blocked):     $cs_blocked")
    println("  LBT Success Rate:   $(cs_attempts > 0 ? round(100.0 * (1 - cs_blocked/cs_attempts), digits=1) : 100.0) %")
    println("  Total Transmitted:  $(cs_attempts - cs_blocked)")
    
    if cs_blocked > 0
        println("\nCarrier Sense Sensing Distance Stats (m):")
        println("  Mean Sensing Dist:  $(round(mean(cs_blocked_distances), digits=1)) m")
        println("  Max Sensing Dist:   $(round(maximum(cs_blocked_distances), digits=1)) m")
        println("  Min Sensing Dist:   $(round(minimum(cs_blocked_distances), digits=1)) m")
        
        # Binned Histogram for CS Blocking
        dist_bins = collect(0:500:2000)
        bin_counts = zeros(Int, length(dist_bins))
        for d in cs_blocked_distances
            b_idx = min(Int(floor(d / 500)) + 1, length(dist_bins))
            bin_counts[b_idx] += 1
        end
        print("  Distribution:       ")
        for (i, b) in enumerate(dist_bins)
            print("$(Int(b))m: $(bin_counts[i]) | ")
        end
        println()
    end
    
    # CSV保存 (簡易版)
    if params.enable_file_output && !isempty(finished_tx)
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        params_str = "N$(params.num_terminals)_SF$(params.spreading_factor)_Slot$(round(Int, params.slot_length_ms))"
        tx_df = DataFrame(
            terminal_id = [r.terminal_id for r in finished_tx],
            start_ms = [r.start_ms for r in finished_tx],
            status = [r.status for r in finished_tx],
            packet_id = [r.packet_id for r in finished_tx]
        )
        CSV.write(joinpath(output_dir, "retrans_tx_log_$(params_str)_$(timestamp).csv"), tx_df)

        # サマリーCSVの保存
        summary_df = DataFrame(
            Metric = [
                "Total Generated", "Dropped (Buffer)", "Original Success", "Original PER",
                "Total Attempts", "Total Collisions", "Attempt PER",
                "Throughput (bps)", "Throughput (Norm)", "Sync Success Rate (%)",
                "Synced PER", "Async PER", "CS Attempts", "CS Blocked", "LBT Success Rate (%)"
            ],
            Value = [
                original_total, buffer_drops, original_success, original_per,
                total_attempts, total_failed_attempts, attempt_per,
                throughput_bps, norm_throughput, (length(synced_terminal_ids) / total_terminals) * 100.0,
                synced_per, async_per, cs_attempts, cs_blocked, (cs_attempts > 0 ? 100.0 * (1 - cs_blocked/cs_attempts) : 100.0)
            ]
        )
        CSV.write(joinpath(output_dir, "summary_$(params_str)_$(timestamp).csv"), summary_df)
    end
    
    return Dict(
        "total_packets" => total_attempts,           # 全送信試行数 (Prop.jlと互換)
        "success" => total_success_attempts,         # 全成功試行数 (Prop.jlと互換)
        "collisions" => total_failed_attempts,       # 全衝突数
        "per" => attempt_per,                       # 試行ベースのPER
        "original_per" => original_per,             # メッセージベースの不達率 (生成数に対するロス)
        "final_reliability_per" => original_per,   # 最終的な不達率 (互換性用)
        "message_loss_rate" => original_per,        # 明示的な名前
        "buffer_drops" => buffer_drops,
        "total_generated" => total_generated_packets,
        "throughput_bps" => throughput_bps,
        "norm_throughput" => norm_throughput,
        "sync_success_rate" => (length(synced_terminal_ids) / total_terminals) * 100.0,
        "synced_per" => synced_per,
        "async_per" => async_per,
        "cs_attempts" => cs_attempts,
        "cs_blocked" => cs_blocked,
        "cs_blocked_distances" => cs_blocked_distances,
        "lbt_success_rate" => (cs_attempts > 0 ? 100.0 * (1 - cs_blocked/cs_attempts) : 100.0),
        "collision_distance_bins" => distance_bins,
        "collision_counts_per_bin" => collision_counts_per_bin,
        "collision_distances" => collision_distances,
        "collision_types" => collision_types
    )
end

"""
    run_integrated_simulation()
    
デフォルトパラメータでシミュレーションを実行する（既存の動作を維持）
"""
function run_integrated_simulation()
    params = create_integrated_params()
    return run_integrated_simulation_with_params(params)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_integrated_simulation()
end
