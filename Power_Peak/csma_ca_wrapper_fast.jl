# ============================================================
# CSMA/CA Simulation Wrapper (Fast)
# 評価スクリプト用の最適化ラッパー関数
# ============================================================

include("csma_ca_sim.jl")

"""
    run_csma_ca_with_params_fast(num_terminals, num_channels, seed)

CSMA/CAシミュレーションを高速実行する

# Arguments
- `num_terminals::Int`: 端末数
- `num_channels::Int`: チャネル数
- `seed::Int`: 乱数シード

# Returns
- `Dict`: 結果辞書 (total_packets, success, collisions, per)
"""
function run_csma_ca_with_params_fast(num_terminals::Int, num_channels::Int, seed::Int)
    # 全出力を抑制
    original_stdout = stdout
    original_stderr = stderr
    redirect_stdout(devnull)
    redirect_stderr(devnull)
    
    try
        Random.seed!(seed)
        
        # パラメータ作成
        sf = 10
        payload_bytes = 10
        data_freq_ghz = 0.92
        ref_pl_db = 20*log10(data_freq_ghz*1e9) - 147.55
        
        params = CSMAParameters(
            13.0,     # tx_power_dbm
            6.0,      # noise_figure_db
            num_terminals,  # num_terminals
            1000.0,   # area_size_m
            sf, payload_bytes,
            num_channels,  # num_channels
            600000.0, # simulation_duration_ms
            0.01,     # duty_cycle
            true,     # enable_carrier_sense
            60000.0,  # mean_event_interval_ms
            2.7,      # pass_loss_exp
            8.0,      # shadowing_std_db
            ref_pl_db # reference_path_loss_db
        )
        
        # ToA計算
        lora_params = create_lora_params(params.sf, params.payload_bytes)
        airtime_ms = calculate_lora_airtime(lora_params)
        
        # 端末配置
        terminals = []
        for i in 1:params.num_terminals
            x = (rand() - 0.5) * params.area_size_m
            y = (rand() - 0.5) * params.area_size_m
            drift_ppm = rand(-20.0:0.1:20.0)
            drift_factor = 1.0 + drift_ppm / 1e6
            push!(terminals, (id=i, x=x, y=y, drift_ppm=drift_ppm, drift_factor=drift_factor))
        end
        
        # CSMA/CA シミュレーション
        cs_threshold_dbm = -80.0
        events = TransmissionEvent[]
        
        for t in terminals
            initial_time = rand() * params.mean_event_interval_ms * 2
            push!(events, TransmissionEvent(initial_time, t, 0))
        end
        
        sort!(events, by=x->x.time_ms, rev=true)
        
        active_tx = TransmissionRecord[]
        finished_tx = TransmissionRecord[]
        next_available_time = Dict{Int, Float64}()
        
        for t in terminals
            next_available_time[t.id] = 0.0
        end
        
        while !isempty(events)
            event = pop!(events)
            curr_t = event.time_ms
            t = event.terminal
            
            filter!(x -> x.end_ms > curr_t, active_tx)
            
            if curr_t < next_available_time[t.id]
                delayed_time = next_available_time[t.id]
                new_event = TransmissionEvent(delayed_time, t, event.attempt)
                idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
                insert!(events, idx, new_event)
                continue
            end
            
            # Carrier Sense
            is_busy = false
            if params.enable_carrier_sense
                for other in active_tx
                    if other.terminal_id == t.id
                        continue
                    end
                    
                    if other.channel != 1
                        continue
                    end
                    
                    dist = sqrt((t.x - other.x)^2 + (t.y - other.y)^2)
                    pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist, 1.0))
                    rssi = other.tx_power_dbm - pl
                    
                    if rssi > cs_threshold_dbm
                        is_busy = true
                        break
                    end
                end
            end
            
            if !is_busy
                end_time = curr_t + airtime_ms
                off_period = airtime_ms * (1.0 / params.duty_cycle - 1.0)
                next_available_time[t.id] = end_time + off_period
                
                dist_gw = sqrt(t.x^2 + t.y^2)
                pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist_gw, 1.0))
                rx_gw = params.tx_power_dbm - pl_gw
                
                channel = rand(1:params.num_channels)
                rec = TransmissionRecord(t.id, curr_t, end_time, params.tx_power_dbm, t.x, t.y, "Success", rx_gw, channel)
                push!(active_tx, rec)
                push!(finished_tx, rec)
                
                interval = -params.mean_event_interval_ms * log(rand()) * t.drift_factor
                next_time = curr_t + interval
                
                if next_time < params.simulation_duration_ms
                    new_event = TransmissionEvent(next_time, t, 0)
                    idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
                    insert!(events, idx, new_event)
                end
            else
                backoff_ms = 10.0 + rand() * 90.0
                new_time = curr_t + backoff_ms
                
                if new_time < params.simulation_duration_ms
                    new_event = TransmissionEvent(new_time, t, event.attempt + 1)
                    idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
                    insert!(events, idx, new_event)
                end
            end
        end
        
        # 衝突判定
        noise_bw_hz = 125000.0
        noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
        (success, collisions) = detect_collisions_sinr(finished_tx, params.sf, noise_power_dbm)
        
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
