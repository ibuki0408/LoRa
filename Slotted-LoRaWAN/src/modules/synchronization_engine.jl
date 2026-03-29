module SynchronizationEngine

using DataFrames, Statistics, Dates, CSV, DSP
using ..Parameters
using ..Utilities
using ..SignalGeneration
using ..PathLoss
using ..NoiseGeneration
using ..SimulationTypes

export perform_hifi_synchronization

"""
    perform_hifi_synchronization(terminals, beacon_times, sig_gen_params, params, output_dir)

メインシミュレーション関数 (オンデマンド同期版)
必要な区間だけ信号を生成して同期判定を行う（メモリ・速度最適化）
"""
function perform_hifi_synchronization(terminals, beacon_times, sig_gen_params, params, output_dir)
    # 結果格納用
    terminal_sync_infos = Dict()
    sync_results = SyncResult[]

    # 同期ログ保存用データフレーム
    sync_log_df = DataFrame(
        terminal_id = Int[],
        status = String[],
        distance_m = Float64[],
        snr_db = Float64[],
        sync_error_ms = Union{Float64, Missing}[],
        start_ms = Union{Float64, Missing}[],
        raw_peak_time = Union{Float64, Missing}[],
        raw_peak_power = Union{Float64, Missing}[],
        dummy = Union{Float64, Missing}[]
    )
    
    # 各端末で処理
    # Ratio for Downsampling (High -> Low)
    downsample_ratio = Int(round(params.tx_sampling_rate_mhz / params.rx_sampling_rate_mhz))

    # ノイズパラメータ
    bw_hz = params.terminal_bw_mhz * 1e6
    noise_dbm = -174 + 10*log10(bw_hz) + params.noise_figure_db
    
    for t in terminals
        # 起動時刻決定: 0 ~ max_startup_delay_ms の一様分布
        startup_ms = rand() * params.max_startup_delay_ms
        
        # コンソール出力 (詳細ログが有効な場合のみ)
        if params.enable_detailed_logs
            println("  [Term $(t.terminal_id)] Startup Time: $(round(startup_ms, digits=4)) ms")
        end
        
        # 受信窓の定義
        window_duration = params.initial_window_duration_ms
        window_end_ms = startup_ms + window_duration
        
        # 信号生成区間 (マージン込み)
        gen_start_ms = max(0.0, startup_ms - 5.0)
        gen_end_ms = window_end_ms + 10.0 # 少し余裕を持つ
        
        # ★ オンデマンド信号生成 (High Rate) ★
        (sig_tx_segment, num_samples_high) = generate_sync_signal_segment(sig_gen_params, beacon_times, gen_start_ms, gen_end_ms)
        
        # ★ ダウンサンプリング (High -> Low) ★
        # use_waveform_phy = true のとき DSP.resample() でアンチエイリアシング付きリサンプリング
        # use_waveform_phy = false (デフォルト) のとき単純間引き（高速）
        sig_rx_segment = if params.use_waveform_phy
            rate_ratio = params.rx_sampling_rate_mhz / params.tx_sampling_rate_mhz
            resample(sig_tx_segment, rate_ratio)
        else
            sig_tx_segment[1:downsample_ratio:end]
        end
        
        # 時間軸 (RXレート)、gen_start_ms 起点
        rx_len = length(sig_rx_segment)
        time_rx = collect(0:rx_len-1) * (1.0 / (params.rx_sampling_rate_mhz * 1e6)) * 1000.0 .+ gen_start_ms
        
        # --- 受信信号処理 ---
        
        # パスロス適用 (5G基準)
        ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
        pl_p = PathLossParameters(t.distance_m, params.sync_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, ref_pl_5g)
        pl_db = calculate_path_loss(pl_p)
        
        # 帯域幅ミスマッチ損失
        # use_waveform_phy = true (DSP.resample使用) の場合は、リサンプリング時のLPFによって
        # 自然に帯域外電力が落ちるため、ここで手動で引くと二重減衰になってしまう
        bw_mismatch_loss_db = if params.use_waveform_phy
            0.0
        else
            10 * log10(params.signal_bw_mhz / params.terminal_bw_mhz)
        end
        
        total_loss_lin = 10^(-(pl_db + t.shadowing_db + bw_mismatch_loss_db)/10)
        
        # 信号減衰
        sig_rx = sig_rx_segment * sqrt(total_loss_lin)
        
        # 信号電力を保存（ノイズ加算前）
        # sig_power = abs2.(sig_rx) # 未使用
        
        # ノイズ加算
        noise = generate_awgn_noise(noise_dbm, length(sig_rx))
        rx_signal_noisy = sig_rx + noise
        rx_power = abs2.(rx_signal_noisy)
        
        # ★ 間欠受信（省電力化）★
        if params.enable_intermittent_rx
            # 全区間でビーコン検出を試みる
            nf_dbm_full = estimate_noise_floor_integrated(rx_power, time_rx; pre_signal_end_ms=params.noise_floor_window_ms)
            # 間欠受信の初回フルサーチでも動的マージンを適用
            actual_margin_full_db = params.detection_margin_db
            thresh_dbm_full = nf_dbm_full + actual_margin_full_db
            thresh_w_full = 10^(thresh_dbm_full/10) * 1e-3
            
            full_cross = detect_crossings_integrated(rx_power, time_rx, thresh_w_full)
            
            if !isempty(full_cross[:end_times])
                # 最初のビーコンが見つかった！
                first_beacon_time = full_cross[:end_times][1]
                
                # 間欠受信マスク作成
                rx_mask = zeros(Bool, length(time_rx))
                
                # 最初のビーコンまでは全て有効（連続受信）
                rx_mask[time_rx .<= first_beacon_time] .= true
                
                # 予測ビーコン時刻周辺のみ有効
                k = 1
                while true
                    pred_time = first_beacon_time + k * params.beacon_interval_ms
                    if pred_time > window_end_ms
                        break
                    end
                    
                    # ±(intermittent_window_ms/2)の窓
                    half_window = params.intermittent_window_ms / 2
                    in_window = (time_rx .>= pred_time - half_window) .& (time_rx .<= pred_time + half_window)
                    rx_mask .|= in_window
                    k += 1
                end
                
                # マスク適用（窓外をノイズレベルに設定）
                noise_level_w = 10^(noise_dbm/10) / 1000
                rx_power[.!rx_mask] .= noise_level_w
            end
        end
        
        # --- 検出ロジック ---
        # 以前は波形PHY時にマージンを下げていたが、元の15dBに戻す
        actual_margin_db = params.detection_margin_db

        if params.enable_intermittent_rx && @isdefined(full_cross) && !isempty(full_cross[:end_times])
            raw_cross = full_cross
        else
            nf_dbm = estimate_noise_floor_integrated(rx_power, time_rx; pre_signal_end_ms=params.noise_floor_window_ms)
            thresh_dbm = nf_dbm + actual_margin_db
            thresh_w = 10^(thresh_dbm/10) * 1e-3
            raw_cross = detect_crossings_integrated(rx_power, time_rx, thresh_w)
        end
        
        # サンプル数フィルタ
        sampling_interval_ms = 1.0 / (params.rx_sampling_rate_mhz * 1e6) * 1000
        min_samples = params.min_samples
        
        filtered_indices = Int[]
        if !isempty(raw_cross[:start_times])
            for i in 1:length(raw_cross[:start_times])
                dur = raw_cross[:end_times][i] - raw_cross[:start_times][i]
                if ceil(dur / sampling_interval_ms) >= min_samples
                    push!(filtered_indices, i)
                end
            end
        end
        
        crossings_filtered = Dict{Symbol, Vector{Float64}}(
            k => (isempty(raw_cross[k]) ? Float64[] : raw_cross[k][filtered_indices])
            for k in keys(raw_cross)
        )
        
        # デバウンス
        final_cross = debounce_integrated(crossings_filtered, params.debounce_time_ms)
        
        # 窓内ビーコン判定
        valid_times = final_cross[:end_times]
        beacons_in_window = filter(t_val -> t_val >= startup_ms && t_val <= window_end_ms, valid_times)
        
        if !isempty(beacons_in_window)
            # 成功
            slot_start_beacon_ms = beacons_in_window[end]
            first_slot_start_ms = slot_start_beacon_ms
            terminal_sync_infos[t.terminal_id] = first_slot_start_ms
            
            signal_dur_ms = params.signal_duration_us / 1000.0
            min_diff = 1e9
            nearest_ideal_time = 0.0
            
            for bt in beacon_times
                ideal_end = bt + signal_dur_ms + (t.distance_m / 3e8 * 1000)
                diff = abs(slot_start_beacon_ms - ideal_end)
                if diff < min_diff
                    min_diff = diff
                    nearest_ideal_time = ideal_end
                end
            end
            
            sync_error = min_diff
            push!(sync_results, SyncResult(t.terminal_id, true, slot_start_beacon_ms, sync_error))
            push!(sync_log_df, (t.terminal_id, "Success", t.distance_m, 10*log10(mean(rx_power)*1000), sync_error, slot_start_beacon_ms, missing, missing, missing))

            if params.enable_detailed_logs
                println("  [Term $(t.terminal_id)] Sync Success: Collected $(length(beacons_in_window)) beacons in $(window_duration)ms window")
            end
            
            # 端末1データ保存
            if t.terminal_id == 1 && params.enable_detailed_logs && params.enable_file_output
                save_end_time = window_end_ms + 200.0
                indices_to_save = findall(x -> x >= startup_ms && x <= save_end_time, time_rx)
                if !isempty(indices_to_save)
                    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
                    out_path = joinpath(output_dir, "integrated_term1_success_$(timestamp).csv")
                    power_vals = rx_power[indices_to_save] * 1000
                    time_vals = time_rx[indices_to_save] .- startup_ms
                    open(out_path, "w") do io
                        println(io, "time_ms,power_mw,power_dbm")
                        for k in 1:length(indices_to_save)
                            p_mw = power_vals[k]
                            p_dbm = 10*log10(p_mw + 1e-20)
                            println(io, "$(time_vals[k]),$p_mw,$p_dbm")
                        end
                    end
                end
            end
        else
            # 失敗
            terminal_sync_infos[t.terminal_id] = nothing
            push!(sync_results, SyncResult(t.terminal_id, false, 0.0, 0.0))
            push!(sync_log_df, (t.terminal_id, "Failed", t.distance_m, 10*log10(mean(rx_power)*1000), missing, missing, missing, missing, missing))
            
            if params.enable_detailed_logs
                println("  [Term $(t.terminal_id)] Sync Failed: No beacons detected")
            end
            
            if t.terminal_id == 1 && params.enable_file_output
                save_end_time = window_end_ms
                indices_to_save = findall(x -> x >= startup_ms && x <= save_end_time, time_rx)
                if !isempty(indices_to_save)
                   timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
                   out_path = joinpath(output_dir, "integrated_term1_fail_$(timestamp).csv")
                   power_vals = rx_power[indices_to_save] * 1000
                   time_vals = time_rx[indices_to_save] .- startup_ms
                   open(out_path, "w") do io
                       println(io, "time_ms,power_mw,power_dbm")
                       for k in 1:length(indices_to_save)
                           p_mw = power_vals[k]
                           p_dbm = 10*log10(p_mw + 1e-20)
                           println(io, "$(time_vals[k]),$p_mw,$p_dbm")
                       end
                   end
                end
            end
        end
    end
    
    # 同期サマリー表示
    synced_terminals = count(x -> x !== nothing, values(terminal_sync_infos))
    sync_rate = (synced_terminals / length(terminals)) * 100
    
    println("\n" * "="^60)
    println("Synchronization Summary:")
    println("  Total Terminals:    $(length(terminals))")
    println("  Synced Terminals:   $synced_terminals")
    println("  Sync Success Rate:  $(round(sync_rate, digits=2)) %")
    println("="^60)
    
    return terminal_sync_infos, sync_log_df
end

end # module
