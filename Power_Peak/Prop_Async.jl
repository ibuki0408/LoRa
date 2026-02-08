using Random, Statistics, Printf, DataFrames, CSV, Plots, Dates, LinearAlgebra, DSP

# ==========================================
# 1. 統合パラメータ構造体
# ==========================================
mutable struct IntegratedParameters
    # === 物理層 (PHY) ===
    signal_duration_us::Float64
    signal_bw_mhz::Float64
    terminal_bw_mhz::Float64
    tx_sampling_rate_mhz::Float64
    rx_sampling_rate_mhz::Float64
    tx_power_dbm::Float64
    noise_figure_db::Float64
    
    # === MAC層 ===
    num_terminals::Int
    area_size_m::Float64
    slot_length_ms::Float64
    packet_airtime_ms::Float64
    enable_carrier_sense::Bool
    cs_threshold_dbm::Float64
    
    # === LoRa固有 ===
    spreading_factor::Int
    lora_payload_bytes::Int
    num_channels::Int  # マルチチャネル数（1-16、AS923準拠）
    
    # === シミュレーション制御 ===
    beacon_interval_ms::Float64
    simulation_duration_ms::Float64
    max_startup_delay_ms::Float64  # 起動時刻の上限（一様分布の最大値）
    duty_cycle::Float64
    mean_event_interval_ms::Float64 # ポアソン分布の平均送信間隔
    
    # === 環境モデル ===
    shadowing_enabled::Bool
    shadowing_std_db::Float64
    pass_loss_exp::Float64
    
    # === Out-of-band同期 ===
    sync_center_freq_ghz::Float64
    data_center_freq_ghz::Float64
    reference_path_loss_db::Float64
    
    # === 同期基地局位置 ===
    sync_bs_x_m::Float64  # 同期信号送信基地局のX座標 (m)
    sync_bs_y_m::Float64  # 同期信号送信基地局のY座標 (m)
    
    # === 同期検出 ===
    sync_observation_duration_ms::Float64
    gw_tx_power_dbm::Float64
    noise_floor_window_ms::Float64
    detection_margin_db::Float64
    min_samples::Int
    debounce_time_ms::Float64
    initial_window_duration_ms::Float64  # スロット開始前の窓開放時間(ms)
    tx_jitter_max_ms::Float64            # 送信開始タイミングのランダムジッタ(ms)
    
    # === 間欠受信（省電力化） ===
    enable_intermittent_rx::Bool         # 間欠受信の有効/無効
    intermittent_window_ms::Float64      # 各ビーコン周辺のサンプリング窓（±1ms → 2ms）
    initial_search_duration_ms::Float64  # 最初のビーコン探索時間（連続受信）
    
    # === 决定的ジッタ（Deterministic Jitter） ===
    use_deterministic_jitter::Bool       # 決定的ジッタの有効化（端末ID依存のオフセット）
    num_jitter_offsets::Int              # オフセット位置の数（例: 20）
    deterministic_jitter_random_ms::Float64  # 微調整用のランダム成分(ms)
    
    collision_model::Symbol              # :sinr (Captureあり) or :overlap (単純衝突)
    
    # === ACK/再送 ===
    enable_ack::Bool                     # ACK機能の有効化
    max_retries::Int                     # 最大再送回数
    ack_timeout_ms::Float64              # ACKタイムアウト
    rx1_delay_ms::Float64                # RX1窓の遅延
    backoff_base_ms::Float64             # 基本バックオフ時間
    
    # === LBT (Listen Before Talk) ===
    lbt_duration_ms::Float64             # LBT期間（ARIB: 5ms）
    lbt_sample_interval_ms::Float64      # サンプリング間隔
    
    # === 出力制御 ===
    enable_file_output::Bool             # ファイル出力の有効/無効
    enable_plot_output::Bool             # プロット生成の有効/無効
    enable_detailed_logs::Bool           # 詳細ログ（同期ログ、送信ログ）の有効/無効
    force_async_mode::Bool               # 強制非同期モード（全端末を同期失敗扱い）
    max_buffer_size::Int                 # パケットバッファ（待ち行列）の最大サイズ
end

function create_integrated_params()
    # ========================================
    # パラメータ設定
    # ========================================
    
    # --- キャリアセンス設定 ---
    enable_cs = true
    
    # --- LoRa設定 (ToA計算に使用) ---
    sf = 8
    payload_bytes = 10
    
    # --- 周波数設定 ---
    sync_freq_ghz = 3.7   # 5G同期信号
    data_freq_ghz = 0.92  # LoRaデータ
    
    # --- ビーコン間隔設定 ---
    beacon_interval_ms = 200.0  # ビーコン送信間隔
    
    # --- ACK/再送設定 ---
    enable_ack = true     # true: 再送あり, false: 再送なし
    
    # --- 比較実験用: 非同期モード ---
    force_async_mode = true  # 非同期版なので常にtrue扱い
    
    return IntegratedParameters(
        # === 物理層 ===
        66.67,    # signal_duration_us
        1.905,      # signal_bw_mhz
        0.125,    # terminal_bw_mhz
        3.84,     # tx_sampling_rate_mhz
        0.256,      # rx_sampling_rate_mhz
        13.0,     # tx_power_dbm
        6.0,      # noise_figure_db
        
        # === MAC層 ===
        100,        # num_terminals
        2000.0,    # area_size_m
        100.0,      # slot_length_ms
        0.0,      # packet_airtime_ms (自動計算)
        true,  # enable_carrier_sense
        -80.0,   # cs_threshold_dbm
        
        # === LoRa固有パラメータ ===
        sf,
        payload_bytes,
        8,        # num_channels
        
        # === シミュレーション制御 ===
        beacon_interval_ms,
        3600000.0, # simulation_duration_ms (60分)
        60000.0,  # max_startup_delay_ms
        0.01,     # duty_cycle (1%)
        30000.0,  # mean_event_interval_ms
        
        # === 環境モデル ===
        true,     # shadowing_enabled
        8.0,      # shadowing_std_db
        2.7,      # pass_loss_exp
        
        # === Out-of-band同期 ===
        sync_freq_ghz,
        data_freq_ghz,
        20*log10(data_freq_ghz*1e9) - 147.55,
        
        # === 同期基地局位置 ===
        0.0,      # sync_bs_x_m
        0.0,      # sync_bs_y_m
        
        # === 同期検出 ===
        65000.0,  # sync_observation_duration_ms
        40.0,     # gw_tx_power_dbm
        10.0,      # noise_floor_window_ms
        15.0,     # detection_margin_db
        3,        # min_samples
        1.0,      # debounce_time_ms
        200.0,    # initial_window_duration_ms
        10.0,    # tx_jitter_max_ms
        
        # === 間欠受信 ===
        true,      # enable_intermittent_rx
        6.0,      # intermittent_window_ms
        60.0,      # initial_search_duration_ms
        
        # === 決定的ジッタ ===
        false,    # use_deterministic_jitter
        20,       # num_jitter_offsets
        5.0,      # deterministic_jitter_random_ms
        
        :sinr,    # collision_model
        
        # === ACK/再送 ===
        enable_ack,
        3,        # max_retries
        2000.0,   # ack_timeout_ms
        1000.0,   # rx1_delay_ms
        1000.0,   # backoff_base_ms
        
        # === LBT (Listen Before Talk) ===
        5.0,      # lbt_duration_ms
        1.0,      # lbt_sample_interval_ms
        
        # === 出力制御 ===
        true,     # enable_file_output
        true,     # enable_plot_output
        true,     # enable_detailed_logs
        force_async_mode,
        10        # max_buffer_size
    )
end

# ==========================================
# 2. 既存モジュールの読み込み
# ==========================================
include("modules/path_loss.jl")
include("modules/shadowing.jl")
include("modules/noise_generation.jl")
include("modules/terminal_deployment.jl")
include("modules/local_clock.jl")
include("modules/lora_airtime.jl")
include("modules/collision_detection.jl")
include("modules/packet_generation.jl")

using .PacketGeneration

# ==========================================
# 3. イベント管理構造体
# ==========================================
struct CandidateSlot
    time_global_ms::Float64
    terminal_node
    slot_index::Int
    channel::Int
end

mutable struct TransmissionRecord
    terminal_id::Int
    start_ms::Float64
    end_ms::Float64
    tx_power_dbm::Float64
    x::Float64
    y::Float64
    status::String
    rx_power_at_gw::Float64
    channel::Int
    retry_count::Int
    is_retransmission::Bool
    original_packet_id::Int
    failure_reason::String
    sinr_db::Float64
    snr_db::Float64
    num_interferers::Int
    cs_detected::Bool
    interferer_ids::String
    interferer_distances::String
end

@enum EventType begin
    TX_START
    TX_END
    ACK_CHECK
    PACKET_GEN
end

struct MACEvent
    time_ms::Float64
    type::EventType
    terminal_id::Int
    packet_id::Int
    data::Any
end

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
# 5. メイン処理：Phase 2 & 3 (MAC実行)
# ==========================================
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
    
    output_dir = joinpath(@__DIR__, "result_Prop_Async_Optimized")
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
                dummy_rec = TransmissionRecord(me_id, curr_t, curr_t, 0.0, me.x_m, me.y_m, "BufferDrop", -Inf, 0, 0, false, pid, "BufferFull", -Inf, -Inf, 0, false, "", "")
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
            
            is_busy = false
            if params.enable_carrier_sense
                for other in active_tx
                    if other.channel == cand.channel && other.start_ms <= curr_t && other.end_ms > curr_t
                        if other.rx_power_at_gw > params.cs_threshold_dbm
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
                
                prev_attempts = filter(r -> r.original_packet_id == pid && r.status != "BufferDrop", all_tx_records)
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
            else
                backoff = 100.0 + rand() * 300.0
                schedule_tx_start!(event_queue, me, pid, curr_t + backoff, params)
            end
            
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
            
            is_successful = (sinr_db >= required_sir_db) && (snr_db >= required_snr_db)
            
            if is_successful
                rec.status = "Success"
                packet_success_map[pid] = true
                if !isempty(terminal_queues[me_id]) && terminal_queues[me_id][1] == pid
                    popfirst!(terminal_queues[me_id])
                end
                if !isempty(terminal_queues[me_id])
                    next_pid = terminal_queues[me_id][1]
                    schedule_tx_start!(event_queue, me, next_pid, max(curr_t + 1.0, next_available_time[me_id]), params)
                else
                    terminal_is_transmitting[me_id] = false
                end
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
                
                if params.enable_ack && rec.retry_count < params.max_retries
                    backoff = rand() * 4000.0
                    schedule_tx_start!(event_queue, me, pid, max(curr_t + backoff, next_available_time[me_id]), params)
                else
                    rec.status = "Dropped"
                    if !isempty(terminal_queues[me_id]) && terminal_queues[me_id][1] == pid
                        popfirst!(terminal_queues[me_id])
                    end
                    if !isempty(terminal_queues[me_id])
                        next_pid = terminal_queues[me_id][1]
                        schedule_tx_start!(event_queue, me, next_pid, max(curr_t + 1.0, next_available_time[me_id]), params)
                    else
                        terminal_is_transmitting[me_id] = false
                    end
                end
            end
        end
    end
    
    finished_tx = all_tx_records
    total_terminals = params.num_terminals
    success_pid_set = Set{Int}([r.original_packet_id for r in finished_tx if r.status == "Success"])
    
    original_total = total_generated_packets
    original_success = length(success_pid_set)
    original_per = original_total > 0 ? (1.0 - original_success / max(original_total, 1)) : 0.0
    
    attempt_tx = filter(r -> r.status != "BufferDrop", finished_tx)
    total_attempts = length(attempt_tx)
    total_success_attempts = count(r -> r.status == "Success", attempt_tx)
    total_failed_attempts = total_attempts - total_success_attempts
    attempt_per = total_attempts > 0 ? (total_failed_attempts / max(total_attempts, 1)) : 0.0
    
    distance_bins = collect(0:100:1000)
    collision_counts_per_bin = zeros(Int, length(distance_bins))
    for dist in collision_distances
        bin_idx = min(Int(floor(dist / 100)) + 1, length(distance_bins))
        collision_counts_per_bin[bin_idx] += 1
    end
    
    total_retransmissions = count(r -> r.is_retransmission && r.status != "BufferDrop", finished_tx)
    retransmission_rate = total_attempts > 0 ? (total_retransmissions / total_attempts) * 100.0 : 0.0
    
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
            retry_count = [r.retry_count for r in finished_tx],
            packet_id = [r.original_packet_id for r in finished_tx]
        )
        CSV.write(joinpath(output_dir, "optimized_tx_log_$(params_str)_$(timestamp).csv"), tx_df)
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
        "total_retransmissions" => total_retransmissions,
        "retransmission_rate" => retransmission_rate,
        "collision_distance_bins" => distance_bins,
        "collision_counts_per_bin" => collision_counts_per_bin
    )
end

function run_integrated_simulation()
    params = create_integrated_params()
    return run_integrated_simulation_with_params(params)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_integrated_simulation()
end
