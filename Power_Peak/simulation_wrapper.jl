# ============================================================
# Parameterized Simulation Wrapper
# 評価スクリプト用のラッパー関数
# ============================================================

# 元のシミュレーションを読み込み
include("integrated_lora_sim.jl")

"""
    run_simulation_with_params(num_terminals, num_channels, seed)

パラメータを指定してシミュレーションを実行する

# Arguments
- `num_terminals::Int`: 端末数
- `num_channels::Int`: チャネル数
- `seed::Int`: 乱数シード

# Returns
- `Dict`: 結果辞書 (total_packets, success, collisions, per)
"""
function run_simulation_with_params(num_terminals::Int, num_channels::Int, seed::Int)
    # グローバル変数を一時的に変更（ハック）
    # より良い方法: create_integrated_params() を引数付きに変更
    
    # 乱数シードを設定
    Random.seed!(seed)
    
    # パラメータを作成（元の関数をコピーして修正）
    sf = 10
    payload_bytes = 10
    sync_freq_ghz = 3.7
    data_freq_ghz = 0.92
    
    params = IntegratedParameters(
        # === 物理層 ===
        66.67, 3.6, 0.125, 7.68, 2.0, 13.0, 6.0,
        
        # === MAC層 ===
        num_terminals,  # ← 可変
        500.0, 400.0, 0.0, 0.1, true, -120.0,
        
        # === LoRa固有 ===
        sf, payload_bytes,
        num_channels,  # ← 可変
        
        # === シミュレーション制御 ===
        20.0, 600000.0, 200.0, 0.01,
        
        # === 環境モデル ===
        true, 0.0, 2.5,
        
        # === Out-of-band同期 ===
        sync_freq_ghz, data_freq_ghz,
        20*log10(sync_freq_ghz*1e9) - 147.55,
        
        # === 同期検出 ===
        500.0, 43.0, 9.0, 10.0, 3, 1.0, 110.0
    )
    
    # シミュレーション実行（出力抑制）
    original_stdout = stdout
    redirect_stdout(devnull)
    
    try
        # 端末配置
        ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
        dep_p = TerminalDeploymentParameters("random_fixed", 0.0, params.num_terminals, params.area_size_m, 10.0, params.area_size_m/2, params.sync_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, ref_pl_5g)
        terminals = deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
        
        # 同期実行
        output_dir = joinpath(@__DIR__, "result_integrated")
        mkpath(output_dir)
        sync_results = perform_hifi_synchronization(params, terminals; output_dir=output_dir)
        
        # スロット生成
        candidates = CandidateSlot[]
        for t in terminals
            start_ms = sync_results[t.terminal_id]
            if start_ms === nothing continue end
            
            slot_offset = rand(0:9)
            curr_ms = start_ms + slot_offset * params.slot_length_ms
            actual_slot_length = params.slot_length_ms * t.clock_drift_factor
            
            idx = 1
            while curr_ms < params.simulation_duration_ms
                channel = rand(1:params.num_channels)
                push!(candidates, CandidateSlot(curr_ms, t, idx, channel))
                curr_ms += actual_slot_length
                idx += 1
            end
        end
        
        sort!(candidates, by=x->x.time_global_ms, rev=true)
        
        # MAC層シミュレーション（簡略版）
        active_tx = TransmissionRecord[]
        finished_tx = TransmissionRecord[]
        next_available_time = Dict{Int, Float64}()
        for t in terminals
            next_available_time[t.terminal_id] = 0.0
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
                    pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(dist)
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
                pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(dist_gw)
                rx_gw = params.tx_power_dbm - pl_gw
                
                rec = TransmissionRecord(me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, me.x_m, me.y_m, "Success", rx_gw, cand.channel)
                push!(active_tx, rec)
                push!(finished_tx, rec)
            end
        end
        
        # 衝突判定
        noise_bw_hz = 125000.0
        noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
        (success, collisions) = detect_collisions_sinr(finished_tx, params.spreading_factor, noise_power_dbm)
        
        redirect_stdout(original_stdout)
        
        return Dict(
            "total_packets" => length(finished_tx),
            "success" => success,
            "collisions" => collisions,
            "per" => isempty(finished_tx) ? 0.0 : collisions / length(finished_tx) * 100
        )
    catch e
        redirect_stdout(original_stdout)
        rethrow(e)
    end
end
