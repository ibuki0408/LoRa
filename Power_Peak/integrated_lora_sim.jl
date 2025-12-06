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
    transmission_prob::Float64
    enable_carrier_sense::Bool
    cs_threshold_dbm::Float64
    
    # === LoRa固有 ===
    spreading_factor::Int
    lora_payload_bytes::Int
    
    # === シミュレーション制御 ===
    beacon_interval_ms::Float64
    simulation_duration_ms::Float64
    max_startup_delay_ms::Float64
    duty_cycle::Float64
    
    # === 環境モデル ===
    shadowing_enabled::Bool
    shadowing_std_db::Float64
    pass_loss_exp::Float64
    
    # === Out-of-band同期 ===
    sync_center_freq_ghz::Float64
    data_center_freq_ghz::Float64
    reference_path_loss_db::Float64
    
    # === 同期検出 ===
    sync_observation_duration_ms::Float64
    gw_tx_power_dbm::Float64
    noise_floor_window_ms::Float64
    detection_margin_db::Float64
    min_samples::Int
    debounce_time_ms::Float64
    initial_wait_ms::Float64
end

function create_integrated_params()
    # ========================================
    # パラメータ設定
    # ========================================
    
    # --- キャリアセンス設定 ---
    # true:  LBT有効（指数バックオフ、最大5回再試行）→ 衝突削減
    # false: 純粋ALOHA（CS無効）→ 衝突発生
    enable_cs = true
    
    # --- LoRa設定 (ToA計算に使用) ---
    sf = 10
    payload_bytes = 10
    
    # --- 周波数設定 ---
    sync_freq_ghz = 3.7   # 5G同期信号
    data_freq_ghz = 0.92  # LoRaデータ
    
    return IntegratedParameters(
        # === 物理層 ===
        66.67,    # signal_duration_us
        3.6,      # signal_bw_mhz
        0.125,    # terminal_bw_mhz
        7.68,     # tx_sampling_rate_mhz
        2.0,      # rx_sampling_rate_mhz
        13.0,     # tx_power_dbm
        6.0,      # noise_figure_db
        
        # === MAC層 ===
        50,        # num_terminals
        500.0,    # area_size_m
        400.0,    # slot_length_ms
        0.0,      # packet_airtime_ms (自動計算)
        0.1,      # transmission_prob
        true,  # enable_carrier_sense (true: LBT有効, false: 純粋ALOHA)
        -120.0,   # cs_threshold_dbm
        
        # === LoRa固有 ===
        sf,
        payload_bytes,
        
        # === シミュレーション制御 ===
        20.0,  # beacon_interval_ms
        600000.0, # simulation_duration_ms (10分)
        200.0,    # max_startup_delay_ms
        0.01,     # duty_cycle (1%)
        
        # === 環境モデル ===
        true,     # shadowing_enabled
        0.0,      # shadowing_std_db
        2.5,      # pass_loss_exp
        
        # === Out-of-band同期 ===
        sync_freq_ghz,
        data_freq_ghz,
        20*log10(sync_freq_ghz*1e9) - 147.55,  # reference_path_loss_db
        
        # === 同期検出 ===
        500.0,    # sync_observation_duration_ms
        43.0,     # gw_tx_power_dbm
        9.0,      # noise_floor_window_ms
        10.0,     # detection_margin_db
        3,        # min_samples
        1.0,      # debounce_time_ms
        110.0     # initial_wait_ms
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
include("modules/lora_airtime.jl")
include("modules/collision_detection.jl")

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
    
    # 1. 理想的な送信信号生成 (高レート) - 5G帯の同期信号
    sig_params = SignalParameters(params.signal_duration_us*1e-6, params.sync_center_freq_ghz*1e9, 
                                  params.signal_bw_mhz*1e6, params.tx_sampling_rate_mhz*1e6, params.gw_tx_power_dbm)
    
    # 同期用の信号を生成（観察時間はパラメータで指定）
    sync_duration_ms = params.sync_observation_duration_ms
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

    # 各端末で受信シミュレーション (5G帯の同期信号)
    for t in terminals
        # A. パスロス・ノイズ付加 (5G帯)
        # 5G帯の基準パスロスを計算 (Friis式)
        ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
        pl_p = PathLossParameters(t.distance_m, params.sync_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, ref_pl_5g)
        pl_db = calculate_path_loss(pl_p)
        total_loss_lin = 10^(-(pl_db + t.shadowing_db)/10)
        sig_rx = sig_tx_low * sqrt(total_loss_lin)
        noise = generate_awgn_noise(noise_dbm, length(sig_rx))
        rx_power = abs2.(sig_rx + noise)
        
        # B. 信号検出 (main_simulation.jl の高精度ロジック)
        nf_dbm = estimate_noise_floor_integrated(rx_power, time_rx; pre_signal_end_ms=params.noise_floor_window_ms)
        thresh_dbm = nf_dbm + params.detection_margin_db
        thresh_w = 10^(thresh_dbm/10) * 1e-3
        
        # 1. 生の閾値超過検出 (線形補間あり)
        raw_cross = detect_crossings_integrated(rx_power, time_rx, thresh_w)
        
        # 2. サンプル数フィルタ
        sampling_interval_ms = 1.0 / (params.rx_sampling_rate_mhz * 1e6) * 1000
        min_samples = params.min_samples
        
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
        final_cross = debounce_integrated(crossings_filtered, params.debounce_time_ms)
        
        # C. 起動遅延を考慮して、最初に掴むべきビーコンを決定
        startup_ms = rand() * params.max_startup_delay_ms
        
        # ★ 起動時刻をコンソール出力 ★
        println("  [Term $(t.terminal_id)] Startup Time: $(round(startup_ms, digits=4)) ms")
        
        # 起動時刻以降で最初に見つかったビーコン (同期ポイントとして end_times を使用: 立ち下がり基準)
        # main_simulation.jl では start_times が補間された正確な上抜け時刻
        valid_times = final_cross[:end_times]
        
        first_beacon_idx = findfirst(t -> t >= startup_ms, valid_times)
        
        if first_beacon_idx !== nothing
            beacon_arrival_ms = valid_times[first_beacon_idx]
            
            # 同期完了！スロット開始時刻を計算
            first_slot_start_ms = beacon_arrival_ms + params.initial_wait_ms
            
            terminal_sync_infos[t.terminal_id] = first_slot_start_ms
            
            # ログ記録
            push!(sync_log_df, (t.terminal_id, "Success", t.distance_m, 10*log10(mean(rx_power)*1000), beacon_arrival_ms, first_slot_start_ms))
            println("  [Term $(t.terminal_id)] Sync Success: Detected at $(round(beacon_arrival_ms, digits=4)) ms (Dist: $(round(t.distance_m, digits=1))m, NF: $(round(nf_dbm, digits=1)) dBm, Thresh: $(round(thresh_dbm, digits=1)) dBm)")

            # ★ 端末1の受信電力データを保存 (同期成功時: 間欠受信モード) ★
            if t.terminal_id == 1
                println("Saving Terminal 1 Received Power Data (Windowed Mode)...")
                timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
                min_len = min(length(time_rx), length(rx_power))
                
                # フィルタリングロジック
                # 1. 初期探索期間: 0 ～ 最初の同期時刻 + マージン
                # 2. 以降: 次のビーコン予定時刻 ± マージン
                window_margin_ms = 1.0 # 前後1msを開く
                
                indices_to_save = Int[]
                
                # 次のビーコン予定時刻を計算
                next_beacon_time = beacon_arrival_ms
                
                for i in 1:min_len
                    t_val = time_rx[i]
                    
                    is_save = false
                    
                    # A. 初期探索期間 (同期確定まではずっとON)
                    if t_val <= (beacon_arrival_ms + window_margin_ms)
                        is_save = true
                    else
                        # B. 間欠受信期間
                        # 現在時刻が「次のビーコン予定時刻」の窓に入っているか？
                        if t_val >= (next_beacon_time - window_margin_ms) && t_val <= (next_beacon_time + window_margin_ms)
                            is_save = true
                        end
                        
                        # 窓を過ぎたら次のターゲットへ更新
                        if t_val > (next_beacon_time + window_margin_ms)
                            next_beacon_time += params.beacon_interval_ms
                        end
                    end
                    
                    if is_save
                        push!(indices_to_save, i)
                    end
                end
                
                power_df = DataFrame(
                    time_ms = time_rx[indices_to_save],
                    power_mw = rx_power[indices_to_save] * 1000,
                    power_dbm = 10*log10.(rx_power[indices_to_save] * 1000)
                )
                mkpath(output_dir)
                out_path = joinpath(output_dir, "integrated_term1_power_$(timestamp).csv")
                CSV.write(out_path, power_df)
                println(" -> Saved to $out_path (Rows: $(length(indices_to_save)) / $min_len)")
            end

        else
            # 同期失敗 (カバレッジ外など)
            terminal_sync_infos[t.terminal_id] = nothing
            
            # ログ記録
            push!(sync_log_df, (t.terminal_id, "Failed", t.distance_m, 10*log10(mean(rx_power)*1000), missing, missing))
            println("  [Term $(t.terminal_id)] Sync Failed (Dist: $(round(t.distance_m, digits=1))m, NF: $(round(nf_dbm, digits=1)) dBm, Thresh: $(round(thresh_dbm, digits=1)) dBm)")

            # ★ 端末1の受信電力データを保存 (同期失敗時: 連続受信モード) ★
            if t.terminal_id == 1
                println("Saving Terminal 1 Received Power Data (Continuous Mode - Sync Failed)...")
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
    
    # ★ LoRa パラメータから ToA を自動計算して設定 ★
    lora_params = create_lora_params(params.spreading_factor, params.lora_payload_bytes)
    params.packet_airtime_ms = calculate_lora_airtime(lora_params)
    
    println("\nLoRa 設定:")
    println("  SF: $(params.spreading_factor)")
    println("  ペイロード: $(params.lora_payload_bytes) bytes")
    println("  計算された ToA: $(round(params.packet_airtime_ms, digits=2)) ms")
    println("  スロット長: $(params.slot_length_ms) ms")
    println()
    
    # 1. 端末配置 (5G帯でのパスロス計算)
    # 5G帯の基準パスロスを計算
    ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
    dep_p = TerminalDeploymentParameters("random_fixed", 0.0, params.num_terminals, params.area_size_m, 10.0, params.area_size_m/2, params.sync_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, ref_pl_5g)
    terminals = deploy_terminals(dep_p, params.shadowing_std_db, params.shadowing_enabled, params.gw_tx_power_dbm)
    
    # 端末情報表示（5G帯のみ）
    println("\n端末配置（5G同期信号）:")
    for t in terminals
        println("• 端末$(t.terminal_id): 位置($(round(t.x_m, digits=1)), $(round(t.y_m, digits=1))) m, 距離$(round(t.distance_m, digits=1)) m")
        println("  - パスロス: $(round(t.path_loss_db, digits=2)) dB")
        println("  - シャドウイング: $(round(t.shadowing_db, digits=1)) dB")
        println("  - 総損失: $(round(t.total_loss_db, digits=2)) dB")
        println("  - 受信電力: $(round(t.rx_power_dbm, digits=1)) dBm")
    end
    
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
        
        # ランダムな初期オフセットを追加（衝突を観察するため）
        random_offset = rand() * params.slot_length_ms
        curr_ms = start_ms + random_offset
        idx = 1
        while curr_ms < params.simulation_duration_ms
            push!(candidates, CandidateSlot(curr_ms, t, idx))
            curr_ms += params.slot_length_ms
            idx += 1
        end
    end
    
    # ★ここが重要：時刻順にソート★
    # ★時刻順にソート (DES用スタックとして使うため、降順にして末尾からpopする)★
    sort!(candidates, by = x -> x.time_global_ms, rev=true)
    println("   -> Total $(length(candidates)) slots scheduled.")
    
    # 4. Phase 3: MAC層実行 (Discrete Event Simulation)
    println("Phase 3: Running MAC Layer (DES Mode)...")
    
    # ノイズフロア計算
    noise_bw_hz = params.terminal_bw_mhz * 1e6
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    
    active_tx = TransmissionRecord[]
    finished_tx = TransmissionRecord[]
    
    # 端末ごとの次回送信可能時刻 (Duty Cycle用)
    next_available_time = Dict{Int, Float64}()
    for t in terminals
        next_available_time[t.terminal_id] = 0.0
    end
    
    # メインループ (イベントキューが空になるまで)
    while !isempty(candidates)
        cand = pop!(candidates) # 末尾(最小時刻)を取得 (O(1))
        
        curr_t = cand.time_global_ms
        me = cand.terminal_node
        
        # A. 終わった通信を掃除 (curr_t 時点で終了しているもの)
        # ※ active_tx には "現在進行中" または "未来に終了する" 送信が残る
        filter!(x -> x.end_ms > curr_t, active_tx)
        
        # B. Duty Cycle チェック
        if curr_t < next_available_time[me.terminal_id]
            continue
        end
        
        # C. キャリアセンス (LBT)
        is_busy = false
        if params.enable_carrier_sense
            # active_tx 内のすべての送信 (過去に開始し、現在まだ終わっていない) との干渉確認
            for other in active_tx
                # other.start_ms <= curr_t は常に真 (DESなので)
                # other.end_ms > curr_t も filter済なので常に真
                # よって、active_tx にあるものは全て「現在送信中」
                
                # 受信電力計算 for CS
                dist = sqrt((me.x_m - other.x)^2 + (me.y_m - other.y)^2)
                pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(dist)
                
                # Shadowing (簡易的)
                if params.shadowing_enabled
                    pl += randn() * params.shadowing_std_db
                end
                
                rssi = other.tx_power_dbm - pl 
                
                if rssi > params.cs_threshold_dbm
                    is_busy = true
                    break
                end
            end
        end
        
        # D. 送信判定
        if !is_busy
            # 送信成功（チャネルアクセス取得）
            dur = params.packet_airtime_ms
            
            # Duty Cycle 計算
            off_period = dur * (1.0 / params.duty_cycle - 1.0)
            next_available_time[me.terminal_id] = curr_t + dur + off_period
            
            # GWでの受信電力
            dist_gw = sqrt(me.x_m^2 + me.y_m^2)
            pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(dist_gw)
            if params.shadowing_enabled
                pl_gw += randn() * params.shadowing_std_db
            end
            rx_gw = params.tx_power_dbm - pl_gw
            
            rec = TransmissionRecord(me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, me.x_m, me.y_m, "Success", rx_gw)
            push!(active_tx, rec)
            push!(finished_tx, rec)
            
        else
            # ビジー検出 → バックオフして再スケジュール
            # 指数バックオフなどのロジックを入れる場合、新しい時刻で再挿入する
            
            # 簡易実装: 10～100ms ランダムバックオフ
            backoff_ms = 10.0 + rand() * 90.0
            new_time = curr_t + backoff_ms
            
            if new_time < params.simulation_duration_ms
                # 新しいイベントを作成
                # SlotIndexは変わらないが時刻が変わる
                new_cand = CandidateSlot(new_time, me, cand.slot_index)
                
                # スタックに挿入 (降順を維持)
                # searchsortedfirst(A, x, rev=true) は A[i] <= x となる最初の場所を返す
                # つまり x より大きい要素の後ろ、x 以下の要素の前
                # [500, 400, 300]. insert 350. sorted logic returns index 3 (value 300).
                # insert at 3 -> [500, 400, 350, 300]. Correct.
                idx = searchsortedfirst(candidates, new_cand, by=x->x.time_global_ms, rev=true)
                insert!(candidates, idx, new_cand)
            end
        end
    end
    
    # 5. Phase 4: 結果解析 (衝突判定)
    println("Phase 4: Analyzing Results...")
    
    # SINR ベースの衝突判定（モジュール化）
    (success, collisions) = detect_collisions_sinr(finished_tx, params.spreading_factor, noise_power_dbm)
    
    println("-"^30)
    println("Result Summary:")
    println("  Total Packets: $(length(finished_tx))")
    println("  Success:       $success")
    println("  Collisions:    $collisions")
    if !isempty(finished_tx)
        println("  PER:           $(round(collisions/length(finished_tx)*100, digits=2)) %")
    else
        println("  PER:           N/A")
    end
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