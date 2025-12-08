# ============================================================
# Optimized Parameterized Simulation Wrapper
# 評価スクリプト用の最適化ラッパー関数
# ============================================================

# 元のシミュレーションを読み込み
include("integrated_lora_sim.jl")

"""
    run_simulation_with_params_fast(num_terminals, num_channels, seed)

パラメータを指定してシミュレーションを高速実行する（出力・ファイルI/O最小化）

# Arguments
- `num_terminals::Int`: 端末数
- `num_channels::Int`: チャネル数
- `seed::Int`: 乱数シード

# Returns
- `Dict`: 結果辞書 (total_packets, success, collisions, per)
"""
function run_simulation_with_params_fast(num_terminals::Int, num_channels::Int, seed::Int)
    # 全出力を抑制
    original_stdout = stdout
    original_stderr = stderr
    redirect_stdout(devnull)
    redirect_stderr(devnull)
    
    try
        # 乱数シードを設定
        Random.seed!(seed)
    
    # パラメータを作成
    sf = 10
    payload_bytes = 10
    sync_freq_ghz = 3.7
    data_freq_ghz = 0.92
    
    params = IntegratedParameters(
        # === 物理層 ===
        66.67, 3.6, 0.125, 7.68, 2.0, 13.0, 6.0,
        
        # === MAC層 ===
        num_terminals,  # ← 可変
        1000.0,    # area_size_m (ユーザーの変更を反映)
        400.0,     # slot_length_ms
        0.0,       # packet_airtime_ms (自動計算)
        0.1,       # transmission_prob
        true,      # enable_carrier_sense
        -80.0,     # cs_threshold_dbm
        
        # === LoRa固有 ===
        sf, payload_bytes,
        num_channels,  # ← 可変
        
        # === シミュレーション制御 ===
        20.0,      # beacon_interval_ms
        600000.0,  # simulation_duration_ms
        30000.0,   # max_startup_delay_ms
        15000.0,   # mean_startup_delay_ms
        0.01,      # duty_cycle
        60000.0,   # mean_event_interval_ms
        
        # === 環境モデル ===
        true,      # shadowing_enabled
        8.0,       # shadowing_std_db
        2.7,       # pass_loss_exp
        
        # === Out-of-band同期 ===
        sync_freq_ghz, data_freq_ghz,
        20*log10(sync_freq_ghz*1e9) - 147.55,
        
        # === 同期検出 ===
        35000.0,   # sync_observation_duration_ms (ポアソン起動対応)
        43.0,      # gw_tx_power_dbm
        9.0,       # noise_floor_window_ms
        10.0,      # detection_margin_db
        2,         # min_samples
        1.0,       # debounce_time_ms
        110.0      # initial_wait_ms
    )
    
    # ToA計算
    lora_params = create_lora_params(params.spreading_factor, params.lora_payload_bytes)
    airtime_ms = calculate_lora_airtime(lora_params)
    
    # パラメータ更新（mutableではないので新しいインスタンス作成）
    params = IntegratedParameters(
        params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
        params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
        params.num_terminals, params.area_size_m, params.slot_length_ms,
        airtime_ms,  # ← 更新
        params.transmission_prob, params.enable_carrier_sense, params.cs_threshold_dbm,
        params.spreading_factor, params.lora_payload_bytes, params.num_channels,
        params.beacon_interval_ms, params.simulation_duration_ms,
        params.max_startup_delay_ms, params.mean_startup_delay_ms,
        params.duty_cycle, params.mean_event_interval_ms,
        params.shadowing_enabled, params.shadowing_std_db, params.pass_loss_exp,
        params.sync_center_freq_ghz, params.data_center_freq_ghz, params.reference_path_loss_db,
        params.sync_observation_duration_ms, params.gw_tx_power_dbm,
        params.noise_floor_window_ms, params.detection_margin_db,
        params.min_samples, params.debounce_time_ms, params.initial_wait_ms
    )
    
    # 端末配置
    ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
    dep_p = TerminalDeploymentParameters(
        "random_fixed", 0.0, params.num_terminals, params.area_size_m, 10.0,
        params.area_size_m/2, params.sync_center_freq_ghz*1e9,
        params.pass_loss_exp, 1.0, ref_pl_5g
    )
    terminals = deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
    
    # 同期実行（ファイル出力なし）
    sig_params = SignalParameters(
        params.signal_duration_us*1e-6, params.sync_center_freq_ghz*1e9,
        params.signal_bw_mhz*1e6, params.tx_sampling_rate_mhz*1e6, params.gw_tx_power_dbm
    )
    
    sync_duration_ms = params.sync_observation_duration_ms
    time_tx, sig_tx_high, _, ideal_beacon_times = generate_periodic_sync_signals(sig_params, params.beacon_interval_ms, sync_duration_ms)
    
    sync_results = Dict{Int, Union{Float64, Nothing}}()
    
    for t in terminals
        # ポアソン分布起動
        startup_ms = -params.mean_startup_delay_ms * log(rand())
        startup_ms = min(startup_ms, params.max_startup_delay_ms)
        
        # 信号受信シミュレーション（簡略版）
        dist = sqrt(t.x_m^2 + t.y_m^2)
        pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist, 1.0))
        
        if params.shadowing_enabled
            pl += randn() * params.shadowing_std_db
        end
        
        rx_power_dbm = params.gw_tx_power_dbm - pl
        noise_floor_dbm = -174 + 10*log10(params.signal_bw_mhz*1e6) + params.noise_figure_db
        detection_threshold_dbm = noise_floor_dbm + params.detection_margin_db
        
        # 同期成功判定（簡略版）
        if rx_power_dbm > detection_threshold_dbm
            # 最初のビーコンを検出
            first_beacon_idx = findfirst(t -> t > startup_ms, ideal_beacon_times)
            if first_beacon_idx !== nothing
                detected_time = ideal_beacon_times[first_beacon_idx]
                first_slot_start_ms = detected_time + params.initial_wait_ms
                sync_results[t.terminal_id] = first_slot_start_ms
            else
                sync_results[t.terminal_id] = nothing
            end
        else
            sync_results[t.terminal_id] = nothing
        end
    end
    
    # スロット生成（ポアソン分布）
    candidates = CandidateSlot[]
    next_available_time = Dict{Int, Float64}()
    
    for t in terminals
        start_ms = sync_results[t.terminal_id]
        if start_ms === nothing
            continue
        end
        
        next_available_time[t.terminal_id] = start_ms
        current_time_ms = start_ms
        
        while current_time_ms < params.simulation_duration_ms
            # ポアソン分布でパケット生成
            next_tx_time = PacketGeneration.generate_next_poisson_time(
                current_time_ms,
                params.mean_event_interval_ms,
                t.clock_drift_factor,
                next_available_time[t.terminal_id]
            )
            
            if next_tx_time >= params.simulation_duration_ms
                break
            end
            
            # スロット境界にスナップ
            slot_len = params.slot_length_ms * t.clock_drift_factor
            slot_idx = Int(floor((next_tx_time - start_ms) / slot_len))
            snapped_time = start_ms + slot_idx * slot_len
            
            if snapped_time < params.simulation_duration_ms
                channel = rand(1:params.num_channels)
                push!(candidates, CandidateSlot(snapped_time, t, slot_idx, channel))
            end
            
            current_time_ms = next_tx_time
        end
    end
    
    sort!(candidates, by=x->x.time_global_ms, rev=true)
    
    # MAC層シミュレーション
    active_tx = TransmissionRecord[]
    finished_tx = TransmissionRecord[]
    
    for t in terminals
        if haskey(next_available_time, t.terminal_id)
            next_available_time[t.terminal_id] = sync_results[t.terminal_id]
        end
    end
    
    while !isempty(candidates)
        cand = pop!(candidates)
        curr_t = cand.time_global_ms
        me = cand.terminal_node
        
        filter!(x -> x.end_ms > curr_t, active_tx)
        
        if curr_t < next_available_time[me.terminal_id]
            continue
        end
        
        # Carrier Sense
        is_busy = false
        if params.enable_carrier_sense
            for other in active_tx
                if other.channel != cand.channel
                    continue
                end
                
                dist = sqrt((me.x_m - other.x)^2 + (me.y_m - other.y)^2)
                pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist, 1.0))
                rssi = other.tx_power_dbm - pl
                
                if rssi > params.cs_threshold_dbm
                    is_busy = true
                    break
                end
            end
        end
        
        if !is_busy
            dur = params.packet_airtime_ms
            off_period = dur * (1.0 / params.duty_cycle - 1.0)
            next_available_time[me.terminal_id] = curr_t + dur + off_period
            
            dist_gw = sqrt(me.x_m^2 + me.y_m^2)
            pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist_gw, 1.0))
            rx_gw = params.tx_power_dbm - pl_gw
            
            rec = TransmissionRecord(me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, me.x_m, me.y_m, "Success", rx_gw, cand.channel)
            push!(active_tx, rec)
            push!(finished_tx, rec)
        else
            # スロットベースバックオフ + チャネル再選択
            backoff_slots = rand(1:5)
            slot_len = params.slot_length_ms * me.clock_drift_factor
            new_time = curr_t + backoff_slots * slot_len
            
            if new_time < params.simulation_duration_ms
                new_channel = rand(1:params.num_channels)
                new_cand = CandidateSlot(new_time, me, cand.slot_index, new_channel)
                idx = searchsortedfirst(candidates, new_cand, by=x->x.time_global_ms, rev=true)
                insert!(candidates, idx, new_cand)
            end
        end
    end
    
    # 衝突判定
    noise_bw_hz = 125000.0
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    (success, collisions) = detect_collisions_sinr(finished_tx, params.spreading_factor, noise_power_dbm)
    
        return Dict(
            "total_packets" => length(finished_tx),
            "success" => success,
            "collisions" => collisions,
            "per" => isempty(finished_tx) ? 0.0 : collisions / length(finished_tx) * 100
        )
    catch e
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        rethrow(e)
    finally
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
    end
end

# 互換性のため元の関数名でもエクスポート
const run_simulation_with_params = run_simulation_with_params_fast
