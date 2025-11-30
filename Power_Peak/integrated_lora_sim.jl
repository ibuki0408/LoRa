using Random, Statistics, Printf, DataFrames, CSV, Plots, Dates, LinearAlgebra, DSP

# ==========================================
# 1. 統合パラメータ構造体
# ==========================================
struct IntegratedParameters
    # --- 信号・PHYパラメータ (main_simulation.jl 由来) ---
    signal_duration_us::Float64
    center_freq_ghz::Float64
    signal_bw_mhz::Float64
    terminal_bw_mhz::Float64
    tx_sampling_rate_mhz::Float64
    rx_sampling_rate_mhz::Float64
    tx_power_dbm::Float64
    noise_figure_db::Float64
    
    # --- ネットワーク・MACパラメータ (advance_lora_sim.jl 由来) ---
    num_terminals::Int
    area_size_m::Float64
    slot_length_ms::Float64
    packet_airtime_ms::Float64
    transmission_prob::Float64
    enable_carrier_sense::Bool
    cs_threshold_dbm::Float64
    
    # --- 制御パラメータ ---
    beacon_interval_ms::Float64
    total_duration_ms::Float64
    max_startup_delay_ms::Float64
    
    # --- 環境 ---
    shadowing_enabled::Bool
    shadowing_std_db::Float64
end

function create_integrated_params()
    return IntegratedParameters(
        # --- PHY (LoRa風物理層パラメータ) ---
        66.67,   # signal_duration_us: 信号長 (μs)
        0.92,    # center_freq_ghz: 中心周波数 (GHz) - 920MHz帯
        3.6,     # signal_bw_mhz: 送信信号帯域幅 (MHz) - 広帯域
        0.125,   # terminal_bw_mhz: 端末受信帯域幅 (MHz) - 125kHz
        7.68,    # tx_sampling_rate_mhz: 送信サンプリングレート (MHz)
        0.25,    # rx_sampling_rate_mhz: 受信サンプリングレート (MHz)
        13.0,    # tx_power_dbm: 送信電力 (dBm)
        5.0,     # noise_figure_db: 受信機雑音指数 (dB)

        # --- MAC (ネットワーク・制御パラメータ) ---
        2,      # num_terminals: 端末数
        1000.0,  # area_size_m: エリアサイズ (m) - 1000m x 1000m
        100.0,   # slot_length_ms: 1スロットの長さ (ms)
        50.0,    # packet_airtime_ms: パケット送信時間 (ms)
        0.1,     # transmission_prob: 送信確率 (各スロットで送信する確率)
        true,    # enable_carrier_sense: キャリアセンス有効化 (true/false)
        -120.0,  # cs_threshold_dbm: キャリアセンス閾値 (dBm)

        # --- Control (シミュレーション制御) ---
        20.0,    # beacon_interval_ms: 同期ビーコン間隔 (ms)
        5000.0,  # total_duration_ms: シミュレーション総時間 (ms)
        200.0,   # max_startup_delay_ms: 最大起動遅延 (ms) - ランダムな起動ズレ

        # --- Env (環境パラメータ) ---
        true,    # shadowing_enabled: シャドウイング有効化
        0.0      # shadowing_std_db: シャドウイング標準偏差 (dB)
    )
end

# ==========================================
# 2. 既存モジュールの読み込み
# ==========================================
include("modules/signal_generation.jl")
include("modules/path_loss.jl")
include("modules/shadowing.jl")
include("modules/noise_generation.jl")
include("modules/terminal_deployment.jl")
include("modules/local_clock.jl")

# main_simulation.jl から重要な関数を再定義・統合
# (依存関係を断ち切るため、必要なロジックをここに移植します)

# --- A. ノイズフロア推定 ---
# --- A. ノイズフロア推定 (main_simulation.jl 準拠) ---
function estimate_noise_floor_integrated(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}; pre_signal_end_ms::Float64=9.0)
    idxs = findall(t -> t < pre_signal_end_ms, time_axis_ms)
    if isempty(idxs)
        # フォールバック: 全区間の中央値
        noise_power_w = median(rx_power)
        return 10 * log10(noise_power_w * 1000)
    end
    segment = rx_power[idxs]
    # ロバストに中央値を使用
    noise_power_w = median(segment)
    return 10 * log10(noise_power_w * 1000)
end

# --- B. 閾値超過検出 (線形補間付き) ---
# --- B. 閾値超過検出 (main_simulation.jl 準拠: 線形補間・詳細記録) ---
function detect_crossings_integrated(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}, threshold_w::Float64)
    crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[],      # 閾値上抜け時刻
        :end_times => Float64[],        # 閾値下抜け時刻
        :peak_times => Float64[],       # ピーク時刻
        :peak_powers_w => Float64[],    # ピーク電力（W）
        :peak_powers_dbm => Float64[]   # ピーク電力（dBm）
    )
    
    is_above_threshold = false
    segment_start_idx = 0
    segment_peak_idx = 0
    segment_peak_power = 0.0
    
    for i in 2:length(rx_power)
        p_prev = rx_power[i-1]
        p_curr = rx_power[i]
        t_prev = time_axis_ms[i-1]
        t_curr = time_axis_ms[i]
        
        # 上抜け検出
        if !is_above_threshold && p_prev <= threshold_w && p_curr > threshold_w
            # 線形補間で正確な上抜け時刻を計算
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            
            push!(crossings[:start_times], crossing_time)
            is_above_threshold = true
            segment_start_idx = i
            segment_peak_idx = i
            segment_peak_power = p_curr
        end
        
        # 閾値超過区間内でピークを追跡
        if is_above_threshold
            if p_curr > segment_peak_power
                segment_peak_idx = i
                segment_peak_power = p_curr
            end
        end
        
        # 下抜け検出
        if is_above_threshold && p_prev > threshold_w && p_curr <= threshold_w
            # 線形補間で正確な下抜け時刻を計算
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            
            push!(crossings[:end_times], crossing_time)
            
            # この区間のピークを記録
            peak_time = time_axis_ms[segment_peak_idx]
            push!(crossings[:peak_times], peak_time)
            push!(crossings[:peak_powers_w], segment_peak_power)
            push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
            
            is_above_threshold = false
        end
    end
    
    # 最後まで閾値を超えていた場合の処理
    if is_above_threshold
        if !isempty(time_axis_ms)
            push!(crossings[:end_times], time_axis_ms[end])
            peak_time = time_axis_ms[segment_peak_idx]
            push!(crossings[:peak_times], peak_time)
            push!(crossings[:peak_powers_w], segment_peak_power)
            push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
        end
    end
    
    return crossings
end

# --- C. デバウンス処理 ---
# --- C. デバウンス処理 (main_simulation.jl 準拠) ---
function debounce_integrated(crossings::Dict{Symbol, Vector{Float64}}, debounce_time_ms::Float64)
    if isempty(crossings[:peak_times])
        return crossings
    end

    filtered_crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[],
        :end_times => Float64[],
        :peak_times => Float64[],
        :peak_powers_w => Float64[],
        :peak_powers_dbm => Float64[]
    )

    last_accepted_peak_time_ms = -Inf

    for i in 1:length(crossings[:peak_times])
        current_peak_time = crossings[:peak_times][i]

        if current_peak_time >= last_accepted_peak_time_ms + debounce_time_ms
            push!(filtered_crossings[:start_times], crossings[:start_times][i])
            push!(filtered_crossings[:end_times], crossings[:end_times][i])
            push!(filtered_crossings[:peak_times], current_peak_time)
            push!(filtered_crossings[:peak_powers_w], crossings[:peak_powers_w][i])
            push!(filtered_crossings[:peak_powers_dbm], crossings[:peak_powers_dbm][i])

            last_accepted_peak_time_ms = current_peak_time
        end
    end

    return filtered_crossings
end

# ==========================================
# 3. イベント管理構造体
# ==========================================
struct CandidateSlot
    time_global_ms::Float64  # グローバル時刻での開始
    terminal_node            # 端末オブジェクト
    slot_index::Int
end

mutable struct TransmissionRecord
    terminal_id::Int
    start_ms::Float64
    end_ms::Float64
    tx_power_dbm::Float64
    x::Float64; y::Float64
    status::String
    rx_power_at_gw::Float64
end

# ==========================================
# 4. メイン処理：Phase 1 (物理層同期)
# ==========================================
"""
main_simulation.jl のロジックを使って、各端末の「最初のスロット開始時刻」を決定する
"""
function perform_hifi_synchronization(params::IntegratedParameters, terminals; output_dir::String="result_integrated")
    println("Phase 1: High-Fidelity Synchronization per Terminal...")
    
    # 1. 理想的な送信信号生成 (高レート)
    sig_params = SignalParameters(params.signal_duration_us*1e-6, params.center_freq_ghz*1e9, 
                                  params.signal_bw_mhz*1e6, params.tx_sampling_rate_mhz*1e6, 43.0) # GW Power 24dBm
    
    # 同期用に最初の少しの期間だけ信号を作る (計算軽量化のため250ms分だけ)
    sync_duration_ms = 250.0 
    time_tx, sig_tx_high, _ = generate_periodic_sync_signals(sig_params, params.beacon_interval_ms, sync_duration_ms)
    
    # 2. ダウンサンプリング
    ratio = params.rx_sampling_rate_mhz / params.tx_sampling_rate_mhz
    sig_tx_low = resample(sig_tx_high, ratio)
    time_rx = collect((0:length(sig_tx_low)-1) / (params.rx_sampling_rate_mhz*1e6) * 1000)
    
    terminal_sync_infos = Dict() # {id => first_slot_start_ms}
    
    # 同期ログ保存用データフレーム
    sync_log_df = DataFrame(
        terminal_id = Int[],
        status = String[],
        distance_m = Float64[],
        rx_power_dbm = Float64[],
        detected_time_ms = Union{Float64, Missing}[],
        slot_start_ms = Union{Float64, Missing}[]
    )

    # ノイズパラメータ
    bw_hz = params.terminal_bw_mhz * 1e6
    noise_dbm = -174 + 10*log10(bw_hz) + params.noise_figure_db
    noise_p = NoiseParameters(noise_dbm, 100.0)
    shad_p = ShadowingParameters(params.shadowing_enabled, params.shadowing_std_db, 50.0, 0.5)

    # 各端末で受信シミュレーション
    for t in terminals
        # A. パスロス・ノイズ付加
        pl_p = PathLossParameters(t.distance_m, params.center_freq_ghz*1e9, 3.0, 1.0, 0.0)
        pl_db = calculate_path_loss(pl_p)
        total_loss_lin = 10^(-(pl_db + t.shadowing_db)/10)
        sig_rx = sig_tx_low * sqrt(total_loss_lin)
        noise = generate_awgn_noise(noise_dbm, length(sig_rx))
        rx_power = abs2.(sig_rx + noise)
        
        # ★ 端末1の受信電力データを保存 (main_simulation.jl 準拠) ★
        if t.terminal_id == 1
            println("Saving Terminal 1 Received Power Data...")
            timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
            min_len = min(length(time_rx), length(rx_power))
            power_df = DataFrame(
                time_ms = time_rx[1:min_len],
                power_mw = rx_power[1:min_len] * 1000,
                power_dbm = 10*log10.(rx_power[1:min_len] * 1000)
            )
            mkpath(output_dir)
            out_path = joinpath(output_dir, "integrated_term1_power_$(timestamp).csv")
            CSV.write(out_path, power_df)
            println(" -> Saved to $out_path")
        end

        # B. 信号検出 (main_simulation.jl の高精度ロジック)
        nf_dbm = estimate_noise_floor_integrated(rx_power, time_rx; pre_signal_end_ms=9.0)
        thresh_w = 10^((nf_dbm + 15.0)/10) * 1e-3 # Margin 15dB
        
        # 1. 生の閾値超過検出 (線形補間あり)
        raw_cross = detect_crossings_integrated(rx_power, time_rx, thresh_w)
        
        # 2. サンプル数フィルタ (2サンプル以上)
        sampling_interval_ms = 1.0 / (params.rx_sampling_rate_mhz * 1e6) * 1000
        min_samples = 2
        
        filtered_indices = Int[]
        for i in 1:length(raw_cross[:start_times])
            dur = raw_cross[:end_times][i] - raw_cross[:start_times][i]
            if ceil(dur / sampling_interval_ms) >= min_samples
                push!(filtered_indices, i)
            end
        end
        
        # フィルタ後の辞書再構築
        crossings_filtered = Dict{Symbol, Vector{Float64}}(
            k => (isempty(raw_cross[k]) ? Float64[] : raw_cross[k][filtered_indices])
            for k in keys(raw_cross)
        )

        # 3. デバウンス処理
        final_cross = debounce_integrated(crossings_filtered, 1.0)
        
        # C. 起動遅延を考慮して、最初に掴むべきビーコンを決定
        startup_ms = rand() * params.max_startup_delay_ms
        
        # 起動時刻以降で最初に見つかったビーコン (同期ポイントとして start_times を使用)
        # main_simulation.jl では start_times が補間された正確な上抜け時刻
        valid_times = final_cross[:start_times]
        
        first_beacon_idx = findfirst(t -> t >= startup_ms, valid_times)
        
        if first_beacon_idx !== nothing
            beacon_arrival_ms = valid_times[first_beacon_idx]
            
            # 同期完了！スロット開始時刻を計算
            initial_wait_ms = 110.0
            first_slot_start_ms = beacon_arrival_ms + initial_wait_ms
            
            terminal_sync_infos[t.terminal_id] = first_slot_start_ms
            
            # ログ記録
            push!(sync_log_df, (t.terminal_id, "Success", t.distance_m, 10*log10(mean(rx_power)*1000), beacon_arrival_ms, first_slot_start_ms))
            println("  [Term $(t.terminal_id)] Sync Success: Detected at $(round(beacon_arrival_ms, digits=4)) ms (Dist: $(round(t.distance_m, digits=1))m)")
        else
            # 同期失敗 (カバレッジ外など)
            terminal_sync_infos[t.terminal_id] = nothing
            
            # ログ記録
            push!(sync_log_df, (t.terminal_id, "Failed", t.distance_m, 10*log10(mean(rx_power)*1000), missing, missing))
            println("  [Term $(t.terminal_id)] Sync Failed (Dist: $(round(t.distance_m, digits=1))m)")
        end
    end
    
    # 同期ログをCSV保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    mkpath(output_dir)
    out_path = joinpath(output_dir, "integrated_sync_log_$(timestamp).csv")
    CSV.write(out_path, sync_log_df)
    println("Sync Log saved to $out_path")
    
    return terminal_sync_infos
end

# ==========================================
# 5. メイン処理：Phase 2 & 3 (MAC実行)
# ==========================================
function run_integrated_simulation()
    println("="^60)
    println("   Integrated LoRa Simulation")
    println("   (Hi-Fi PHY Sync -> Time-Sorted MAC)")
    println("="^60)
    
    params = create_integrated_params()
    
    # 1. 端末配置
    dep_p = TerminalDeploymentParameters("random_fixed", 0.0, params.num_terminals, params.area_size_m, 10.0, params.area_size_m/2, params.center_freq_ghz*1e9, 3.0, 1.0, 0.0)
    terminals = deploy_terminals(dep_p, 0.0, false, params.tx_power_dbm)
    
    # 2. Phase 1: 高精度同期を実行
    # ここで main_simulation.jl のロジックが走り、各端末の開始時刻が決まる
    # 実行場所に関わらず、スクリプトのあるフォルダ(Power_Peak)の下に保存する
    output_dir = joinpath(@__DIR__, "result_integrated")
    mkpath(output_dir)
    sync_results = perform_hifi_synchronization(params, terminals; output_dir=output_dir)
    
    # 3. Phase 2: スロット候補の生成とソート
    println("Phase 2: Scheduling & Sorting...")
    candidates = CandidateSlot[]
    
    for t in terminals
        start_ms = sync_results[t.terminal_id]
        if start_ms === nothing continue end # 同期失敗端末はスキップ
        
        # シミュレーション終了までスロットを生成
        curr_ms = start_ms
        idx = 1
        while curr_ms < params.total_duration_ms
            push!(candidates, CandidateSlot(curr_ms, t, idx))
            curr_ms += params.slot_length_ms
            idx += 1
        end
    end
    
    # ★ここが重要：時刻順にソート★
    sort!(candidates, by = x -> x.time_global_ms)
    println("   -> Total $(length(candidates)) slots scheduled.")
    
    # 4. Phase 3: MAC層シミュレーション (CSMA/CA)
    println("Phase 3: Running MAC Layer...")
    
    active_tx = TransmissionRecord[]
    finished_tx = TransmissionRecord[]
    
    for cand in candidates
        curr_t = cand.time_global_ms
        me = cand.terminal_node
        
        # A. 終わった通信を掃除
        filter!(x -> x.end_ms > curr_t, active_tx)
        
        # B. キャリアセンス (Advance_LoRa ロジック)
        is_busy = false
        if params.enable_carrier_sense
            for other in active_tx
                # 簡易RSSI計算
                dist = sqrt((me.x_m - other.x)^2 + (me.y_m - other.y)^2)
                pl = 20*log10(dist) + 20*log10(params.center_freq_ghz*1e9) - 147.55
                rssi = other.tx_power_dbm - pl # 簡易PathLoss
                
                if rssi > params.cs_threshold_dbm
                    is_busy = true
                    break
                end
            end
        end
        
        # C. 送信判定
        if !is_busy
            if rand() < params.transmission_prob
                # 送信実行
                dur = params.packet_airtime_ms
                
                # GWでの受信電力 (判定用)
                dist_gw = sqrt(me.x_m^2 + me.y_m^2)
                pl_gw = 20*log10(dist_gw) + 20*log10(params.center_freq_ghz*1e9) - 147.55
                rx_gw = params.tx_power_dbm - pl_gw
                
                rec = TransmissionRecord(me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, me.x_m, me.y_m, "Success", rx_gw)
                push!(active_tx, rec)
                push!(finished_tx, rec)
            end
        end
    end
    
    # 5. Phase 4: 結果解析 (衝突判定)
    println("Phase 4: Analyzing Results...")
    
    collisions = 0
    success = 0
    
    # 簡易衝突判定 (時間重複 & 電力差なし)
    for i in 1:length(finished_tx)
        p1 = finished_tx[i]
        is_collided = false
        
        for j in 1:length(finished_tx)
            if i == j continue end
            p2 = finished_tx[j]
            
            # 時間重複
            if max(p1.start_ms, p2.start_ms) < min(p1.end_ms, p2.end_ms)
                # Capture Effect (6dB)
                if p1.rx_power_at_gw < p2.rx_power_at_gw + 6.0
                    is_collided = true # 負けた、または引き分け
                end
            end
        end
        
        if is_collided
            p1.status = "Collision"
            collisions += 1
        else
            success += 1
        end
    end
    
    println("-"^30)
    println("Result Summary:")
    println("  Total Packets: $(length(finished_tx))")
    println("  Success:       $success")
    println("  Collisions:    $collisions")
    println("  PER:           $(round(collisions/length(finished_tx)*100, digits=2)) %")
    println("-"^30)

    # ★ 送信結果のCSV保存 ★
    if !isempty(finished_tx)
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        
        # 1. パケット詳細リスト
        tx_df = DataFrame(
            terminal_id = [r.terminal_id for r in finished_tx],
            start_ms = [r.start_ms for r in finished_tx],
            end_ms = [r.end_ms for r in finished_tx],
            status = [r.status for r in finished_tx],
            tx_power_dbm = [r.tx_power_dbm for r in finished_tx],
            rx_power_gw_dbm = [r.rx_power_at_gw for r in finished_tx],
            x_m = [r.x for r in finished_tx],
            y_m = [r.y for r in finished_tx]
        )
        out_tx = joinpath(output_dir, "integrated_tx_log_$(timestamp).csv")
        CSV.write(out_tx, tx_df)
        println("Tx Log saved to $out_tx")
        
        # 2. 端末ごとの集計
        term_ids = sort(unique(tx_df.terminal_id))
        summary_df = DataFrame(
            terminal_id = Int[],
            total_tx = Int[],
            success = Int[],
            collision = Int[],
            per_percent = Float64[]
        )
        
        for tid in term_ids
            sub = filter(row -> row.terminal_id == tid, tx_df)
            tot = nrow(sub)
            suc = count(x -> x == "Success", sub.status)
            col = count(x -> x == "Collision", sub.status)
            per = tot > 0 ? (col / tot) * 100 : 0.0
            push!(summary_df, (tid, tot, suc, col, per))
        end
        out_sum = joinpath(output_dir, "integrated_summary_$(timestamp).csv")
        CSV.write(out_sum, summary_df)
        println("Summary saved to $out_sum")
    end
    
    # 6. 可視化
    if length(finished_tx) > 0
        p = plot(title="Integrated LoRa Sim", xlabel="Time (ms)", ylabel="Node ID")
        limit_ms = 2000.0
        subset = filter(x -> x.end_ms < limit_ms, finished_tx)
        
        for r in subset
            c = r.status == "Success" ? :blue : :red
            plot!(p, [r.start_ms, r.end_ms, r.end_ms, r.start_ms], [r.terminal_id-0.4, r.terminal_id-0.4, r.terminal_id+0.4, r.terminal_id+0.4], seriestype=:shape, color=c, opacity=0.5, label="", linecolor=:white)
        end
        out_png = joinpath(output_dir, "integrated_result.png")
        savefig(p, out_png)
        println("Plot saved to $out_png")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_integrated_simulation()
end