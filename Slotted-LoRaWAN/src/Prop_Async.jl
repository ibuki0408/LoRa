using Random, Statistics, Printf, DataFrames, CSV, Plots, Dates, LinearAlgebra, DSP

# ==========================================
include(joinpath(@__DIR__, "SlottedLoRaWAN.jl"))
using .SlottedLoRaWAN

# Redundant includes removed (now part of SlottedLoRaWAN)

using .PacketGeneration

# Redundant types removed (now defined in SimulationTypes via SlottedLoRaWAN)

# ==========================================
# 4. ヘルパー関数: パケット送信のスケジュール
# ==========================================
"""
    schedule_tx_start!(event_queue, me, packet_id, base_time, params)

非同期モード用: 指定された base_time で送信イベントをスケジュールする。
（同期スロットへのスナップ処理を削除）
"""
function schedule_tx_start!(event_queue, me, packet_id, base_time, params)
    # 非同期端末: base_time をそのまま使用
    tx_time_final = base_time
    
    if tx_time_final < params.simulation_duration_ms
        channel = rand(1:params.num_channels)
        new_cand = CandidateSlot(tx_time_final, me, 1, channel)
        new_event = MACEvent(tx_time_final, TX_START, me.terminal_id, packet_id, new_cand)
        
        idx = searchsortedfirst(event_queue, new_event, by=x->x.time_ms, rev=true)
        insert!(event_queue, idx, new_event)
    end
end

# ==========================================
# Protocol Scenario Settings
# ==========================================
# (Removed local const, now using params.enable_carrier_sense)

function run_integrated_simulation_with_params(params::IntegratedParameters)
    println("="^60)
    println("   Optimized Integrated LoRa Simulation (Async-Only)")
    println("="^60)
    
    if params.packet_airtime_ms == 0.0
        lora_params = create_lora_params(params.spreading_factor, params.lora_payload_bytes)
        params.packet_airtime_ms = calculate_lora_airtime(lora_params)
    end
    
    if params.slot_length_ms == 0.0
        params.slot_length_ms = params.packet_airtime_ms + 100.0
    end
    
    println("\nLoRa 設定:")
    println("  SF: $(params.spreading_factor)")
    println("  ペイロード: $(params.lora_payload_bytes) bytes")
    println("  計算された ToA: $(round(params.packet_airtime_ms, digits=2)) ms")
    println()
    
    ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
    dep_p = TerminalDeploymentParameters(
        "random_fixed", 0.0, params.num_terminals, params.area_size_m, 10.0, params.area_size_m/2, 
        params.sync_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, ref_pl_5g,
        params.sync_bs_x_m, params.sync_bs_y_m
    )
    terminals = deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
    
    output_dir = joinpath(@__DIR__, "result_Prop_Async")
    mkpath(output_dir)
    
    println("Phase 1: Initializing Asynchronous Mode...")
    # 同期処理をスキップ
    
    println("Phase 2: Initializing Queues & Scheduling Packet Generation...")
    event_queue = MACEvent[]
    next_packet_id = 1
    
    terminal_queues = Dict{Int, Vector{Int}}(t.terminal_id => Int[] for t in terminals)
    terminal_is_transmitting = Dict{Int, Bool}(t.terminal_id => false for t in terminals)
    next_available_time = Dict{Int, Float64}(t.terminal_id => 0.0 for t in terminals)
    
    for t in terminals
        startup_ms = rand() * params.max_startup_delay_ms
        initial_interval = -params.mean_event_interval_ms * log(rand()) * t.clock_drift_factor
        first_gen_time = startup_ms + initial_interval
        
        if first_gen_time < params.simulation_duration_ms
            new_event = MACEvent(first_gen_time, PACKET_GEN, t.terminal_id, 0, nothing)
            push!(event_queue, new_event)
        end
    end
    sort!(event_queue, by = x -> x.time_ms, rev=true)
    
    function push_event!(q, ev)
        idx = searchsortedfirst(q, ev, by=x->x.time_ms, rev=true)
        insert!(q, idx, ev)
    end

    println("Phase 3: Running MAC Layer...")
    noise_bw_hz = params.terminal_bw_mhz * 1e6
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    
    required_sir_db = 6.0
    required_snr_db = -15.0 # SF10 default
    if params.spreading_factor == 7 required_snr_db = -7.5
    elseif params.spreading_factor == 8 required_snr_db = -10.0
    elseif params.spreading_factor == 9 required_snr_db = -12.5
    elseif params.spreading_factor == 11 required_snr_db = -17.5
    elseif params.spreading_factor == 12 required_snr_db = -20.0
    end
    
    active_tx = TransmissionRecord[]
    all_tx_records = TransmissionRecord[]
    packet_success_map = Dict{Int, Bool}()
    
    buffer_drops = 0
    total_generated_packets = 0
    collision_distances = Float64[]
    
    while !isempty(event_queue)
        event = pop!(event_queue)
        curr_t = event.time_ms
        me_id = event.terminal_id
        me = terminals[me_id]
        
        if event.type == PACKET_GEN
            total_generated_packets += 1
            pid = next_packet_id
            next_packet_id += 1
            
            if length(terminal_queues[me_id]) < params.max_buffer_size
                push!(terminal_queues[me_id], pid)
                if !terminal_is_transmitting[me_id]
                    terminal_is_transmitting[me_id] = true
                    base_tx_time = max(curr_t, next_available_time[me_id])
                    schedule_tx_start!(event_queue, me, pid, base_tx_time, params)
                end
            else
                buffer_drops += 1
                dummy_rec = TransmissionRecord(me_id, curr_t, curr_t, 0.0, me.x_m, me.y_m, "BufferDrop", -Inf, 0, pid, "BufferFull", -Inf, -Inf, 0, false, "", "")
                push!(all_tx_records, dummy_rec)
            end
            
            next_gen_interval = -params.mean_event_interval_ms * log(rand()) * me.clock_drift_factor
            next_gen_time = curr_t + next_gen_interval
            if next_gen_time < params.simulation_duration_ms
                push_event!(event_queue, MACEvent(next_gen_time, PACKET_GEN, me_id, 0, nothing))
            end
            
        elseif event.type == TX_START
            cand = event.data
            pid = event.packet_id
            
            if curr_t < next_available_time[me_id]
                schedule_tx_start!(event_queue, me, pid, next_available_time[me_id], params)
                continue
            end
            
            # キャリアセンス (LBT) - CS-ALOHA設定時
            is_busy = false
            if params.enable_carrier_sense
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
                        sh_peer = randn() * params.shadowing_std_db
                        p_peer_dbm = params.tx_power_dbm - (pl_peer + sh_peer)
                        
                        if p_peer_dbm >= params.cs_threshold_dbm
                            is_busy = true
                            break
                        end
                    end
                end
            end
            
            if !is_busy
                dur = params.packet_airtime_ms
                dist_gw = sqrt(me.x_m^2 + me.y_m^2)
                #pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist_gw, 1.0))
                pl_gw=0.0
                rx_gw = params.tx_power_dbm - (pl_gw + me.shadowing_db)
                
                rec = TransmissionRecord(
                    me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, 
                    me.x_m, me.y_m, "Processing", rx_gw, cand.channel,
                    0, false, pid,
                    "None", 0.0, 0.0, 0, false, "", ""
                )
                
                push!(active_tx, rec)
                push!(all_tx_records, rec)
                
                push_event!(event_queue, MACEvent(curr_t + dur, TX_END, me_id, pid, rec))
                
                off_period = dur * (1.0 / params.duty_cycle - 1.0)
                next_available_time[me_id] = curr_t + dur + off_period
            else
                backoff = 100.0 + rand() * 300.0
                schedule_tx_start!(event_queue, me, pid, curr_t + backoff, params)
            end
            
        elseif event.type == TX_END
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
            
            if params.enable_capture_effect
                is_successful = (sinr_db >= required_sir_db) && (snr_db >= required_snr_db)
            else
                is_successful = isempty(interferers)
            end
            
            rec.sinr_db = sinr_db
            rec.snr_db = snr_db
            rec.num_interferers = length(interferers)
            
            if is_successful
                rec.status = "Success"
                packet_success_map[pid] = true
            else
                rec.status = "Collision"
                for inf in interferers
                    inf_term_id = inf.terminal_id
                    if inf_term_id != me_id
                        dx = terminals[me_id].x_m - terminals[inf_term_id].x_m
                        dy = terminals[me_id].y_m - terminals[inf_term_id].y_m
                        push!(collision_distances, sqrt(dx^2 + dy^2))
                    end
                end
            end

            # パケットをキューから削除（再送なしなので成功・失敗に関わらず削除）
            if !isempty(terminal_queues[me_id]) && terminal_queues[me_id][1] == pid
                popfirst!(terminal_queues[me_id])
            end
            
            # 次のパケットがあればスケジュール
            if !isempty(terminal_queues[me_id])
                next_pid = terminal_queues[me_id][1]
                schedule_tx_start!(event_queue, me, next_pid, max(curr_t + 1.0, next_available_time[me_id]), params)
            else
                terminal_is_transmitting[me_id] = false
            end
        end
    end
    
    finished_tx = all_tx_records
    total_terminals = params.num_terminals
    success_pid_set = Set{Int}([r.packet_id for r in finished_tx if r.status == "Success"])
    
    original_total = total_generated_packets
    original_success = length(success_pid_set)
    original_per = original_total > 0 ? (1.0 - original_success / max(original_total, 1)) : 0.0
    
    attempt_tx = filter(r -> r.status != "BufferDrop", finished_tx)
    total_attempts = length(attempt_tx)
    total_success_attempts = count(r -> r.status == "Success", attempt_tx)
    total_failed_attempts = total_attempts - total_success_attempts
    attempt_per = total_attempts > 0 ? (total_failed_attempts / max(total_attempts, 1)) : 0.0
    
    sim_duration_s = params.simulation_duration_ms / 1000.0
    total_data_bits = original_success * params.lora_payload_bytes * 8
    throughput_bps = total_data_bits / sim_duration_s
    norm_throughput = (original_success * params.packet_airtime_ms) / (params.simulation_duration_ms * params.num_channels)
    
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
    
    if params.enable_detailed_logs && params.enable_file_output && !isempty(finished_tx)
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        params_str = "N$(params.num_terminals)_SF$(params.spreading_factor)"
        tx_df = DataFrame(
            terminal_id = [r.terminal_id for r in finished_tx],
            start_ms = [r.start_ms for r in finished_tx],
            status = [r.status for r in finished_tx],
            packet_id = [r.packet_id for r in finished_tx]
        )
        CSV.write(joinpath(output_dir, "optimized_tx_log_$(params_str)_$(timestamp).csv"), tx_df)

        # サマリーCSVの保存
        summary_df = DataFrame(
            Metric = [
                "Total Generated", "Dropped (Buffer)", "Original Success", "Original PER",
                "Total Attempts", "Total Collisions", "Attempt PER",
                "Throughput (bps)", "Throughput (Norm)"
            ],
            Value = [
                original_total, buffer_drops, original_success, original_per,
                total_attempts, total_failed_attempts, attempt_per,
                throughput_bps, norm_throughput
            ]
        )
        CSV.write(joinpath(output_dir, "summary_$(params_str)_$(timestamp).csv"), summary_df)
    end
    
    return Dict(
        "total_packets" => total_attempts,
        "success" => total_success_attempts,
        "collisions" => total_failed_attempts,
        "per" => attempt_per,
        "original_per" => original_per,
        "final_reliability_per" => original_per,
        "message_loss_rate" => original_per,
        "buffer_drops" => buffer_drops,
        "total_generated" => total_generated_packets,
        "throughput_bps" => throughput_bps,
        "norm_throughput" => norm_throughput,
        "sync_success_rate" => 0.0,
        "synced_per" => 0.0,
        "async_per" => 0.0,
        "total_retransmissions" => 0,
        "retransmission_rate" => 0.0,
        "collision_distance_bins" => [],
        "collision_counts_per_bin" => []
    )
end

function run_integrated_simulation()
    params = create_integrated_params()
    return run_integrated_simulation_with_params(params)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_integrated_simulation()
end
