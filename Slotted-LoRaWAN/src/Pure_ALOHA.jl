using Random, Statistics, Printf, DataFrames, CSV, Dates, LinearAlgebra, DSP

# ==========================================
# 1. ライブラリ読み込み
# ==========================================
include(joinpath(@__DIR__, "SlottedLoRaWAN.jl"))
using .SlottedLoRaWAN

using .PacketGeneration

# Redundant types removed (now defined in SimulationTypes via SlottedLoRaWAN)

# ==========================================
# 4. ヘルパー関数: パケット送信のスケジュール
# ==========================================
function schedule_tx_start!(event_queue, me, packet_id, base_time, params)
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
    println("   Pure ALOHA Simulation (Async-Only)")
    println("="^60)
    
    if params.packet_airtime_ms == 0.0
        lora_params = create_lora_params(params.spreading_factor, params.lora_payload_bytes)
        params.packet_airtime_ms = calculate_lora_airtime(lora_params)
    end
    
    if params.slot_length_ms == 0.0
        params.slot_length_ms = params.packet_airtime_ms + 100.0
    end
    
    ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
    dep_p = TerminalDeploymentParameters(
        "random_fixed", 0.0, params.num_terminals, params.area_size_m, 10.0, params.area_size_m/2, 
        params.sync_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, ref_pl_5g,
        params.sync_bs_x_m, params.sync_bs_y_m
    )
    terminals = deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
    
    output_dir = joinpath(@__DIR__, "result_Pure_ALOHA")
    mkpath(output_dir)
    
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

    noise_bw_hz = params.terminal_bw_mhz * 1e6
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    
    required_sir_db = 6.0
    required_snr_db = -10.0 # SF8
    
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
                    schedule_tx_start!(event_queue, me, pid, max(curr_t, next_available_time[me_id]), params)
                end
            else
                buffer_drops += 1
            end
            
            next_gen_interval = -params.mean_event_interval_ms * log(rand()) * me.clock_drift_factor
            next_gen_time = curr_t + next_gen_interval
            if next_gen_time < params.simulation_duration_ms
                push_event!(event_queue, MACEvent(next_gen_time, PACKET_GEN, me_id, 0, nothing))
            end
            
        elseif event.type == TX_START
            cand = event.data
            pid = event.packet_id
            
            # キャリアセンス (LBT) - CS-ALOHA設定時
            is_busy = false
            if params.enable_carrier_sense
                # 端末間RSSベースの正確なLBT
                for other in active_tx
                    if other.channel == cand.channel
                        other_term = terminals[other.terminal_id]
                        dx = me.x_m - other_term.x_m
                        dy = me.y_m - other_term.y_m
                        dist_peer = sqrt(dx^2 + dy^2)
                        
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
            
            if is_busy
                # Busyならバックオフしてリトライ
                backoff = 100.0 + rand() * 300.0
                schedule_tx_start!(event_queue, me, pid, curr_t + backoff, params)
                continue
            end
            
            # 送信開始
            dur = params.packet_airtime_ms
            dist_gw = sqrt(me.x_m^2 + me.y_m^2)
            pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist_gw, 1.0))
            rx_gw = params.tx_power_dbm - (pl_gw + me.shadowing_db)
            
            prev_attempts = filter(r -> r.packet_id == pid, all_tx_records)
            current_retry = isempty(prev_attempts) ? 0 : maximum([r.retry_count for r in prev_attempts]) + 1
            
            rec = TransmissionRecord(
                me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, 
                me.x_m, me.y_m, "Processing", rx_gw, cand.channel,
                current_retry, current_retry > 0, pid,
                "None", 0.0, 0.0, 0, false, "", ""
            )
            
            push!(active_tx, rec)
            push!(all_tx_records, rec)
            
            push_event!(event_queue, MACEvent(curr_t + dur, TX_END, me_id, pid, rec))
            push_event!(event_queue, MACEvent(curr_t + dur + 1.0, ACK_CHECK, me_id, pid, rec))
            
            off_period = dur * (1.0 / params.duty_cycle - 1.0)
            next_available_time[me_id] = curr_t + dur + off_period
            
        elseif event.type == TX_END
            rec = event.data
            filter!(r -> r !== rec, active_tx)
            
        elseif event.type == ACK_CHECK
            rec = event.data
            pid = event.packet_id
            
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
            
            if is_successful
                rec.status = "Success"
                packet_success_map[pid] = true
                popfirst!(terminal_queues[me_id])
                if !isempty(terminal_queues[me_id])
                    schedule_tx_start!(event_queue, me, terminal_queues[me_id][1], max(curr_t + 1.0, next_available_time[me_id]), params)
                else
                    terminal_is_transmitting[me_id] = false
                end
            else
                rec.status = "Collision"
                if params.enable_ack && rec.retry_count < params.max_retries
                    backoff = rand() * 4000.0
                    schedule_tx_start!(event_queue, me, pid, max(curr_t + backoff, next_available_time[me_id]), params)
                else
                    rec.status = "Dropped"
                    popfirst!(terminal_queues[me_id])
                    if !isempty(terminal_queues[me_id])
                        schedule_tx_start!(event_queue, me, terminal_queues[me_id][1], max(curr_t + 1.0, next_available_time[me_id]), params)
                    else
                        terminal_is_transmitting[me_id] = false
                    end
                end
            end
        end
    end
    
    finished_tx = all_tx_records
    original_success = length(Set([r.packet_id for r in finished_tx if r.status == "Success"]))
    total_generated = total_generated_packets
    
    println("\nSimulation Summary (Pure ALOHA):")
    println("  Total Generated:     $total_generated")
    println("  Original Success:    $original_success")
    val_per = (1 - original_success/total_generated) * 100
    println("  Original PER:        $(round(val_per, digits=4)) %")
    
    # CSV保存
    if params.enable_file_output && !isempty(finished_tx)
        output_dir = joinpath(@__DIR__, "result_PureALOHA")
        mkpath(output_dir)
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        params_str = "N$(params.num_terminals)_SF$(params.spreading_factor)"
        
        tx_df = DataFrame(
            terminal_id = [r.terminal_id for r in finished_tx],
            start_ms = [r.start_ms for r in finished_tx],
            status = [r.status for r in finished_tx],
            packet_id = [r.packet_id for r in finished_tx]
        )
        CSV.write(joinpath(output_dir, "purealoha_tx_log_$(params_str)_$(timestamp).csv"), tx_df)

        # サマリーCSVの保存
        summary_df = DataFrame(
            Metric = [
                "Total Generated", "Dropped (Buffer)", "Original Success", "Original PER (%)",
                "Throughput (bps)", "Throughput (Norm)"
            ],
            Value = [
                total_generated, buffer_drops, original_success, val_per,
                (original_success * params.lora_payload_bytes * 8) / (params.simulation_duration_ms / 1000.0),
                (original_success * params.packet_airtime_ms) / (params.simulation_duration_ms * params.num_channels)
            ]
        )
        CSV.write(joinpath(output_dir, "summary_$(params_str)_$(timestamp).csv"), summary_df)
    end
    
    return Dict(
        "total_packets" => length(finished_tx),
        "success" => count(r -> r.status == "Success", finished_tx),
        "collisions" => count(r -> r.status == "Collision" || r.status == "Dropped", finished_tx),
        "per" => (1 - original_success/total_generated)*100,
        "original_per" => (1 - original_success/total_generated)*100,
        "throughput_bps" => (original_success * params.lora_payload_bytes * 8) / (params.simulation_duration_ms / 1000.0),
        "norm_throughput" => (original_success * params.packet_airtime_ms) / (params.simulation_duration_ms * params.num_channels),
        "sync_success_rate" => 0.0,
        "total_retransmissions" => count(r -> r.is_retransmission, finished_tx),
        "buffer_drops" => buffer_drops,
        "total_generated" => total_generated
    )
end

function run_integrated_simulation()
    params = create_integrated_params()
    return run_integrated_simulation_with_params(params)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_integrated_simulation()
end
