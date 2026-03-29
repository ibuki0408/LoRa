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
    
    # === 決定的ジッタ（Deterministic Jitter） ===
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
    sf = 8
    payload_bytes = 10
    
    # --- ビーコン間隔設定 ---
    beacon_interval_ms = 100.0  # ビーコン送信間隔
    
    # --- 周波数設定 ---
    sync_freq_ghz = 3.7   # 5G同期信号
    data_freq_ghz = 0.92  # LoRaデータ
    
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
        400,        # num_terminals
        2000.0,    # area_size_m
        200.0,      # slot_length_ms (ビーコン間隔に合わせる)
        0.0,      # packet_airtime_ms (自動計算)
        true,  # enable_carrier_sense (true: LBT有効, false: 純粋ALOHA)
        -80.0,   # cs_threshold_dbm
        
        # === LoRa固有パラメータ ===
        sf,
        payload_bytes,
        8,        # num_channels (AS923 Japan typical (1-16))
        
        # === シミュレーション制御 ===
        beacon_interval_ms,     # beacon_interval_ms
        3600000.0, # simulation_duration_ms (60分)
        60000.0,  # max_startup_delay_ms (最大60秒)
        0.01,     # duty_cycle (1%)
        60000.0,  # mean_event_interval_ms (平均60秒間隔 = 少し余裕を持たせる)
        
        # === 環境モデル ===
        true,     # shadowing_enabled
        8.0,      # shadowing_std_db
        2.7,      # pass_loss_exp
        
        # === Out-of-band同期 ===
        sync_freq_ghz,
        data_freq_ghz,
        20*log10(data_freq_ghz*1e9) - 147.55,  # reference_path_loss_db (データ通信用)
        
        # === 同期基地局位置 ===
        0.0,      # sync_bs_x_m (デフォルト: 原点 = GWと同じ位置)
        0.0,      # sync_bs_y_m (デフォルト: 原点 = GWと同じ位置)
        
        # === 同期検出 ===
        65000.0,  # sync_observation_duration_ms (最大起動時刻60s + マージン5s)
        40.0,     # gw_tx_power_dbm
        10.0,      # noise_floor_window_ms
        15.0,     # detection_margin_db
        3,        # min_samples
        1.0,      # debounce_time_ms
        beacon_interval_ms * 2.0,    # initial_window_duration_ms (ビーコン間隔の2倍)
        10.0,     # tx_jitter_max_ms (Pattern C: optimal low jitter)
        
        # === 間欠受信 ===
        true,      # enable_intermittent_rx
        6.0,      # intermittent_window_ms (±6ms) - ジッタに対応するため拡大
        60.0,      # initial_search_duration_ms
        
        # === 決定的ジッタ ===
        false,    # use_deterministic_jitter (デフォルト: 無効 = ランダムジッタ)
        20,       # num_jitter_offsets (20端末で分散)
        5.0,      # deterministic_jitter_random_ms (微調整用)
        
        :sinr,    # collision_model (SINR with capture effect)
        
        # === ACK/再送 ===
        true,     # enable_ack (テストのため有効化)
        3,        # max_retries
        2000.0,   # ack_timeout_ms
        1000.0,   # rx1_delay_ms
        1000.0,   # backoff_base_ms
        
        # === LBT (Listen Before Talk) ===
        5.0,      # lbt_duration_ms (ARIB STD-T108: 5ms)
        1.0,      # lbt_sample_interval_ms (1ms間隔でサンプリング)
        
        # === 出力制御 ===
        true,     # enable_file_output (デフォルト: 有効)
        true,     # enable_plot_output (デフォルト: 有効)
        true      # enable_detailed_logs (デフォルト: 有効)
    )
end

# ==========================================
# 2. 既存モジュールの読み込み
# ==========================================
include("../src/modules/signal_generation.jl")
include("../src/modules/path_loss.jl")
include("../src/modules/shadowing.jl")
include("../src/modules/noise_generation.jl")
include("../src/modules/terminal_deployment.jl")
include("../src/modules/local_clock.jl")
include("../src/modules/lora_airtime.jl")
include("../src/modules/collision_detection.jl")
include("../src/modules/packet_generation.jl")

using .PacketGeneration


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
struct SyncResult
    terminal_id::Int
    success::Bool
    sync_time_ms::Float64
    sync_error_ms::Float64
end

struct CandidateSlot
    time_global_ms::Float64  # グローバル時刻での開始
    terminal_node            # 端末オブジェクト
    slot_index::Int
    channel::Int             # チャネル番号（1-based）
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
    channel::Int  # チャネル番号（1-based）
    retry_count::Int              # 現在の再送回数
    is_retransmission::Bool       # 再送パケットか
    original_packet_id::Int       # 元のパケットID（再送追跡用）
    # === パケットロス原因分析用フィールド ===
    failure_reason::String        # 詳細な失敗理由
    sinr_db::Float64             # 実際のSINR値 (dB)
    snr_db::Float64              # 実際のSNR値 (dB)
    num_interferers::Int         # 干渉パケット数
    cs_detected::Bool            # キャリアセンスで干渉検出したか
    interferer_ids::String       # 干渉した端末IDリスト (例: "1;23;45")
    interferer_distances::String # 干渉した端末との距離 [m] リスト (例: "120.5;340.2")
end

# ==========================================
# 4. メイン処理：Phase 1 (物理層同期)

# ==========================================
# 共通ロジック
# ==========================================
"""
main_simulation.jl のロジックを使って、各端末の「最初のスロット開始時刻」を決定する

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

    # 同期ログ保存用データフレーム（後でCSVにする）
    # 大量になるので、メモリ展開せずに都度書き出し or 必要なデータだけ保持
    # ここではシンプルに結果配列を返す形に変更したので、内部で処理する
    
    # 各端末で処理
    # Ratio for Downsampling (High -> Low)
    downsample_ratio = Int(round(params.tx_sampling_rate_mhz / params.rx_sampling_rate_mhz))

    # ノイズパラメータ
    bw_hz = params.terminal_bw_mhz * 1e6
    noise_dbm = -174 + 10*log10(bw_hz) + params.noise_figure_db
    # noise_p = NoiseParameters(noise_dbm, 100.0) # 未使用
    
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
        sig_rx_segment = sig_tx_segment[1:downsample_ratio:end]
        
        # 時間軸 (RXレート)、gen_start_ms 起点
        rx_len = length(sig_rx_segment)
        time_rx = collect(0:rx_len-1) * (1.0 / (params.rx_sampling_rate_mhz * 1e6)) * 1000.0 .+ gen_start_ms
        
        # --- 受信信号処理 ---
        
        # パスロス適用 (5G基準)
        ref_pl_5g = 20*log10(params.sync_center_freq_ghz*1e9) + 20*log10(1.0) - 147.55
        pl_p = PathLossParameters(t.distance_m, params.sync_center_freq_ghz*1e9, params.pass_loss_exp, 1.0, ref_pl_5g)
        pl_db = calculate_path_loss(pl_p)
        bw_mismatch_loss_db = 10 * log10(params.signal_bw_mhz / params.terminal_bw_mhz)
        total_loss_lin = 10^(-(pl_db + t.shadowing_db + bw_mismatch_loss_db)/10)
        
        # 信号減衰
        sig_rx = sig_rx_segment * sqrt(total_loss_lin)
        
        # 信号電力を保存（ノイズ加算前）
        sig_power = abs2.(sig_rx)
        
        # ノイズ加算
        noise = generate_awgn_noise(noise_dbm, length(sig_rx))
        rx_signal_noisy = sig_rx + noise
        rx_power = abs2.(rx_signal_noisy)
        
        # ★ 間欠受信（省電力化）★
        if params.enable_intermittent_rx
            # 全区間でビーコン検出を試みる
            nf_dbm_full = estimate_noise_floor_integrated(rx_power, time_rx; pre_signal_end_ms=params.noise_floor_window_ms)
            thresh_dbm_full = nf_dbm_full + params.detection_margin_db
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
                
                if params.enable_detailed_logs
                    active_samples = sum(rx_mask)
                    total_samples = length(rx_mask)
                    reduction_pct = (1 - active_samples / total_samples) * 100
                    first_beacon_rel_time = first_beacon_time - startup_ms
                    println("  [Term $(t.terminal_id)] Intermittent RX: First beacon at $(round(first_beacon_rel_time, digits=1))ms, $(active_samples)/$(total_samples) samples ($(round(reduction_pct, digits=1))% reduction)")
                end
            end
            # 最初のビーコンが見つからない場合は、全区間を使用（マスクなし）
        end
        
        # --- 検出ロジック ---
        
        # 間欠受信が有効で、既にビーコンを検出している場合は、それを使用
        if params.enable_intermittent_rx && @isdefined(full_cross) && !isempty(full_cross[:end_times])
            # 既に検出済みのビーコンを使用
            raw_cross = full_cross
        else
            # 通常の検出ロジック
            # ノイズフロア推定
            nf_dbm = estimate_noise_floor_integrated(rx_power, time_rx; pre_signal_end_ms=params.noise_floor_window_ms)
            thresh_dbm = nf_dbm + params.detection_margin_db
            thresh_w = 10^(thresh_dbm/10) * 1e-3
            
            # クロッシング検出
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
            
            # 理想時刻との誤差
            # 最も近い理想ビーコンを探す
            # 理想ビーコン時刻は周期的に計算可能だが、jitterが入っている場合は beacon_times を参照
            # ただし beacon_times は送信開始時刻なので、信号長を考慮するか、detectが立ち下がりなら信号長を足す
            # detect_crossings_integrated は "立ち下がり(end_times)" をスロット開始としている
            # 理想的な終了時刻 = beacon_start + signal_duration
            
            signal_dur_ms = params.signal_duration_us / 1000.0
            min_diff = 1e9
            nearest_ideal_time = 0.0
            
            for bt in beacon_times
                ideal_end = bt + signal_dur_ms + (t.distance_m / 3e8 * 1000) # 伝搬遅延含む？
                # original logic compared detect time with "ideal beacon time"
                # Let's verify original logic later. For now, simple closest match.
                diff = abs(slot_start_beacon_ms - ideal_end)
                if diff < min_diff
                    min_diff = diff
                    nearest_ideal_time = ideal_end
                end
            end
            
            sync_error = min_diff # 符号なし誤差にする (or slot_start_beacon_ms - nearest_ideal_time)
            
            push!(sync_results, SyncResult(t.terminal_id, true, slot_start_beacon_ms, sync_error))
            
            # Log Success
            push!(sync_log_df, (t.terminal_id, "Success", t.distance_m, 10*log10(mean(rx_power)*1000), sync_error, slot_start_beacon_ms, missing, missing, missing))

            if params.enable_detailed_logs
                println("  [Term $(t.terminal_id)] Sync Success: Collected $(length(beacons_in_window)) beacons in $(window_duration)ms window")
            end
            
            # 端末1データ保存
            if t.terminal_id == 1
                println("Saving Terminal 1 Received Power Data (Window + 200ms)...")
                save_end_time = window_end_ms + 200.0
                indices_to_save = findall(x -> x >= startup_ms && x <= save_end_time, time_rx)
                
                # 詳細ログが有効な場合のみ保存
                if params.enable_detailed_logs && params.enable_file_output
                    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
                    out_path = joinpath(output_dir, "integrated_term1_success_$(timestamp).csv")
                    
                    # IO効率化のため直接書き込みなどを検討だが、一旟DataFrame
                    # 必要な部分だけ切り出す
                    if !isempty(indices_to_save)
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
                       println(" -> Saved to $out_path")
                    end
                end
            end
            
        else
            # 失敗
            terminal_sync_infos[t.terminal_id] = nothing
            push!(sync_results, SyncResult(t.terminal_id, false, 0.0, 0.0))
            
            # Log Failure
            push!(sync_log_df, (t.terminal_id, "Failed", t.distance_m, 10*log10(mean(rx_power)*1000), missing, missing, missing, missing, missing))
            
            if params.enable_detailed_logs
                println("  [Term $(t.terminal_id)] Sync Failed: No beacons detected")
            end
            
            # 失敗時ログ (Term 1)
            if t.terminal_id == 1
                println("Saving Terminal 1 Received Power Data (Fail)...")
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
                   println(" -> Saved to $out_path")
                end
            end
        end
        
        # メモリ解放
        sig_tx_segment = nothing
        sig_rx_segment = nothing
        rx_signal_noisy = nothing
        rx_power = nothing
        # [CRITICAL OPTIMIZATION] 手動GCを削除。非常に重い。
        # GC.gc()
    end
    
    # 同期ログをCSV保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    mkpath(output_dir)
    # 詳細ログが有効な場合のみ保存
    if params.enable_detailed_logs && params.enable_file_output
        out_path = joinpath(output_dir, "integrated_sync_log_$(timestamp).csv")
        CSV.write(out_path, sync_log_df)
        println("Sync Log saved to $out_path")
    end
    

    
    # 同期成功率を表示
    total_terminals = length(terminals)
    synced_terminals = count(x -> x !== nothing, values(terminal_sync_infos))
    sync_rate = (synced_terminals / total_terminals) * 100
    
    println("\n" * "="^60)
    println("Synchronization Summary:")
    println("  Total Terminals:    $total_terminals")
    println("  Synced Terminals:   $synced_terminals")
    println("  Failed Terminals:   $(total_terminals - synced_terminals)")
    println("  Sync Success Rate:  $(round(sync_rate, digits=2)) %")
    println("="^60)
    
    # 同期精度分析（成功した端末のみ）
    synced_data = filter(row -> row.status == "Success", sync_log_df)
    
    if nrow(synced_data) > 0
        errors = synced_data.sync_error_ms
        distances = synced_data.distance_m
        snrs = synced_data.snr_db
        
        # 統計量計算
        mean_error = mean(errors)
        std_error = std(errors)
        min_error = minimum(errors)
        max_error = maximum(errors)
        
        # 相関係数計算
        cor_distance = cor(errors, distances)
        cor_snr = cor(errors, snrs)
        
        println("\n" * "="^60)
        println("Synchronization Accuracy Analysis:")
        println("  Mean Sync Error:    $(round(mean_error, digits=4)) ms")
        println("  Std Dev:            $(round(std_error, digits=4)) ms")
        println("  Min Error:          $(round(min_error, digits=4)) ms")
        println("  Max Error:          $(round(max_error, digits=4)) ms")
        println("")
        println("  Correlation with Distance: $(round(cor_distance, digits=3))")
        println("  Correlation with SNR:      $(round(cor_snr, digits=3))")
        println("="^60)
        
        # 精度サマリーをCSV保存
        accuracy_summary = DataFrame(
            metric = ["mean_error_ms", "std_dev_ms", "min_error_ms", "max_error_ms", "correlation_distance", "correlation_snr"],
            value = [mean_error, std_error, min_error, max_error, cor_distance, cor_snr]
        )
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        # ファイル出力が有効な場合のみ保存
        if params.enable_file_output
            accuracy_path = joinpath(output_dir, "integrated_sync_accuracy_$(timestamp).csv")
            CSV.write(accuracy_path, accuracy_summary)
            println("Accuracy Summary saved to $accuracy_path")
        end
    end
    
    return terminal_sync_infos, sync_log_df
end

# ==========================================
# ヘルパー関数: 再送スケジューリング
# ==========================================
function schedule_retransmission!(queue, packet, params, current_time)
    # バックオフ時間計算（固定 + ランダム）
    backoff_ms = params.backoff_base_ms + rand() * 4000.0
  # 1000-1500ms
    
    # Duty Cycle考慮
    airtime = packet.end_ms - packet.start_ms
    dc_wait = airtime * (1.0 / params.duty_cycle - 1.0)
    next_available = packet.end_ms + dc_wait
    
    # 再送時刻決定（バックオフとDCの遅い方）
    retry_time = max(current_time + backoff_time, next_available)
    
    # 新しいパケット作成
    retry_packet = TransmissionRecord(
        packet.terminal_id,
        retry_time,
        retry_time + airtime,
        packet.tx_power_dbm,
        packet.x, packet.y,
        "Pending",
        packet.rx_power_at_gw,
        packet.channel,
        packet.retry_count + 1,
        true,  # is_retransmission
        packet.original_packet_id,
        "None",  # failure_reason
        0.0,     # sinr_db
        0.0,     # snr_db
        0,       # num_interferers
        false,   # cs_detected
        "",      # interferer_ids
        ""       # interferer_distances
    )
    
    # キューに追加（時刻順にソート）
    push!(queue, retry_packet)
    sort!(queue, by = x -> x.start_ms)
end

# ==========================================
# 5. メイン処理：Phase 2 & 3 (MAC実行)
# ==========================================
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
        params = IntegratedParameters(
            params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
            params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
            params.num_terminals, params.area_size_m, params.slot_length_ms,
            calculate_lora_airtime(lora_params),  # ← ToAを計算
            params.enable_carrier_sense, params.cs_threshold_dbm,
            params.spreading_factor, params.lora_payload_bytes, params.num_channels,
            params.beacon_interval_ms, params.simulation_duration_ms,
            params.max_startup_delay_ms,
            params.duty_cycle, params.mean_event_interval_ms,
            params.shadowing_enabled, params.shadowing_std_db, params.pass_loss_exp,
            params.sync_center_freq_ghz, params.data_center_freq_ghz, params.reference_path_loss_db,
            params.sync_bs_x_m, params.sync_bs_y_m,  # ← BS位置
            params.sync_observation_duration_ms, params.gw_tx_power_dbm,
            params.noise_floor_window_ms, params.detection_margin_db,
            params.min_samples, params.debounce_time_ms, params.initial_window_duration_ms,
            params.tx_jitter_max_ms,
            params.enable_intermittent_rx, params.intermittent_window_ms, params.initial_search_duration_ms,
            params.use_deterministic_jitter, params.num_jitter_offsets, params.deterministic_jitter_random_ms,
            params.collision_model,
            params.enable_ack, params.max_retries, params.ack_timeout_ms,
            params.rx1_delay_ms, params.backoff_base_ms,
            params.lbt_duration_ms, params.lbt_sample_interval_ms,
            params.enable_file_output, params.enable_plot_output, params.enable_detailed_logs
        )
    end
    
    # ★ スロット長が未設定の場合 (0.0) は ToA + 100ms に設定 ★
    if params.slot_length_ms == 0.0
        params = IntegratedParameters(
            params.signal_duration_us, params.signal_bw_mhz, params.terminal_bw_mhz,
            params.tx_sampling_rate_mhz, params.rx_sampling_rate_mhz, params.tx_power_dbm, params.noise_figure_db,
            params.num_terminals, params.area_size_m, 
            params.packet_airtime_ms + 100.0,  # slot_length_ms = ToA + 100ms
            params.packet_airtime_ms,
            params.enable_carrier_sense, params.cs_threshold_dbm,
            params.spreading_factor, params.lora_payload_bytes, params.num_channels,
            params.beacon_interval_ms, params.simulation_duration_ms,
            params.max_startup_delay_ms,
            params.duty_cycle, params.mean_event_interval_ms,
            params.shadowing_enabled, params.shadowing_std_db, params.pass_loss_exp,
            params.sync_center_freq_ghz, params.data_center_freq_ghz, params.reference_path_loss_db,
            params.sync_bs_x_m, params.sync_bs_y_m,
            params.sync_observation_duration_ms, params.gw_tx_power_dbm,
            params.noise_floor_window_ms, params.detection_margin_db,
            params.min_samples, params.debounce_time_ms, params.initial_window_duration_ms,
            params.tx_jitter_max_ms,
            params.enable_intermittent_rx, params.intermittent_window_ms, params.initial_search_duration_ms,
            params.use_deterministic_jitter, params.num_jitter_offsets, params.deterministic_jitter_random_ms,
            params.collision_model,
            params.enable_ack, params.max_retries, params.ack_timeout_ms,
            params.rx1_delay_ms, params.backoff_base_ms,
            params.lbt_duration_ms, params.lbt_sample_interval_ms,
            params.enable_file_output, params.enable_plot_output, params.enable_detailed_logs
        )
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
    
    # ★ クロックドリフトの影響を排除（すべての端末のドリフトを0にする） ★
    for t in terminals
        # struct TerminalInfo is immutable, so we might need to recreate if it was immutable,
        # but let's check modules/terminal_deployment.jl first.
        # Actually it's better to just swap them in the array if needed.
        # Wait, TerminalInfo is a struct (immutable by default in Julia).
    end
    # Re-creating terminals with zero drift
    terminals = [TerminalInfo(t.terminal_id, t.x_m, t.y_m, t.distance_m, t.path_loss_db, t.shadowing_db, t.total_loss_db, t.rx_power_dbm, 0.0, 1.0) for t in terminals]
    
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
    sync_results, sync_log_df = perform_hifi_synchronization(terminals, beacon_times, sig_gen_params, params, output_dir)
    
    # 3. Phase 2: スロット候補の生成とソート
    println("Phase 2: Scheduling & Sorting...")
    candidates = CandidateSlot[]
    
    for t in terminals
        start_ms = sync_results[t.terminal_id]
        
        # クロックドリフトを考慮したToAとオフ期間
        actual_airtime = params.packet_airtime_ms * t.clock_drift_factor
        min_off_period = actual_airtime * (1.0 / params.duty_cycle - 1.0)
        
        if start_ms === nothing
            # ★ 同期失敗端末: 非同期送信（LBT-ALOHAと同じ動作）★
            println("  Terminal $(t.terminal_id): Sync failed - using async transmission")
            
            # ランダムな起動時刻
            async_start = rand() * params.max_startup_delay_ms
            curr_ms = async_start
            idx = 1
            
            # 最初のパケット送信（ポアソン分布）
            initial_interval = -params.mean_event_interval_ms * log(rand()) * t.clock_drift_factor
            first_tx_time = curr_ms + initial_interval
            
            # ★ジッタなし（LBT-ALOHAと同じ動作）★
            # 公平な比較のため、同期失敗端末はジッタを追加しない
            first_tx_final = first_tx_time

            
            if first_tx_final < params.simulation_duration_ms
                channel = rand(1:params.num_channels)
                push!(candidates, CandidateSlot(first_tx_final, t, idx, channel))
                curr_ms = first_tx_final
                idx += 1
            else
                continue
            end
            
            # 次回送信可能時刻
            next_available = curr_ms + actual_airtime + min_off_period
            
            # 以降のパケット送信
            while true
                # ポアソン過程で次の送信時刻を生成
                actual_next_tx = generate_next_poisson_time(curr_ms, params.mean_event_interval_ms, t.clock_drift_factor, next_available)
                
                # ★非同期なのでスロット境界にスナップしない★
                # ★ジッタなし（LBT-ALOHAと同じ動作）★
                final_tx_time = actual_next_tx
                
                if final_tx_time >= params.simulation_duration_ms
                    break
                end
                
                channel = rand(1:params.num_channels)
                push!(candidates, CandidateSlot(final_tx_time, t, idx, channel))
                
                curr_ms = final_tx_time
                next_available = curr_ms + actual_airtime + min_off_period
                idx += 1
            end
            
            continue  # 次の端末へ
        end
        
        # ★ 同期成功端末: スロット境界に揃える ★
        curr_ms = start_ms
        idx = 1
        
        # スロット長（クロックドリフト考慮）
        slot_len = params.slot_length_ms * t.clock_drift_factor
        
        # ★ 最初のパケット送信もポアソン過程に従う ★
        # スロット境界で即座に送信するのではなく、ランダムな待機時間を設定
        initial_interval = -params.mean_event_interval_ms * log(rand()) * t.clock_drift_factor
        
        # 最初の送信時刻（ポアソン分布）
        first_tx_time = curr_ms + initial_interval
        
        # スロット境界にスナップ
        diff_from_base = first_tx_time - start_ms
        num_slots = round(max(0.0, diff_from_base) / slot_len)
        first_tx_snapped = start_ms + num_slots * slot_len
        
        # ★ジッタをスロット内に制限（スロット境界を越えないようにする）★
        # ガードバンド: パケット送信が完了する前にスロットが終わらないようにする
        guard_band_ms = 10.0
        max_jitter_in_slot = max(0.0, params.slot_length_ms - params.packet_airtime_ms - guard_band_ms)
        
        # ジッタを追加（スロット内制限 & tx_jitter_max_ms 制限）
        if params.use_deterministic_jitter
            base_offset = (t.terminal_id % params.num_jitter_offsets) * (max_jitter_in_slot / params.num_jitter_offsets)
            jitter_limit = min(params.deterministic_jitter_random_ms, max_jitter_in_slot - base_offset)
            jitter = base_offset + rand() * jitter_limit
        else
            jitter_limit = min(params.tx_jitter_max_ms, max_jitter_in_slot)
            jitter = rand() * jitter_limit
        end
        first_tx_final = first_tx_snapped + jitter
        
        if first_tx_final < params.simulation_duration_ms
            # ランダムチャネル選択
            channel = rand(1:params.num_channels)
            push!(candidates, CandidateSlot(first_tx_final, t, idx, channel))
            
            # 次回送信可能時刻 (DC明け)
            next_available = first_tx_final + actual_airtime + min_off_period
            curr_ms = first_tx_final
            idx += 1
        else
            # 最初のパケットがシミュレーション時間外なら終了
            continue
        end
        
        while true
            # 次の送信時刻（ポアソン分布）
            actual_next_tx = generate_next_poisson_time(curr_ms, params.mean_event_interval_ms, t.clock_drift_factor, next_available)
            
            # スロット境界にスナップ
            diff_from_base = actual_next_tx - start_ms
            num_slots = ceil(max(0.0, diff_from_base) / slot_len) 
            
            next_tx_snapped = start_ms + num_slots * slot_len
            
            # ★ジッタをスロット内に制限（スロット境界を越えないようにする）★
            guard_band_ms = 10.0
            max_jitter_in_slot = max(0.0, params.slot_length_ms - params.packet_airtime_ms - guard_band_ms)
            
            # スロット境界にジッタ（スロット内制限 & tx_jitter_max_ms 制限）を加算
            if params.use_deterministic_jitter
                base_offset = (t.terminal_id % params.num_jitter_offsets) * (max_jitter_in_slot / params.num_jitter_offsets)
                jitter_limit = min(params.deterministic_jitter_random_ms, max_jitter_in_slot - base_offset)
                jitter = base_offset + rand() * jitter_limit
            else
                jitter_limit = min(params.tx_jitter_max_ms, max_jitter_in_slot)
                jitter = rand() * jitter_limit
            end
            next_tx_final = next_tx_snapped + jitter
            
            if next_tx_final >= params.simulation_duration_ms
                break
            end
            
            curr_ms = next_tx_final
            idx += 1
            
            channel = rand(1:params.num_channels)
            push!(candidates, CandidateSlot(curr_ms, t, idx, channel))
            
            # 次のDC明け更新
            next_available = curr_ms + actual_airtime + min_off_period
        end
    end
    
    # ★ここが重要:時刻順にソート★
    # ★時刻順にソート (DES用スタックとして使うため、降順にして末尾からpopする)★
    sort!(candidates, by = x -> x.time_global_ms, rev=true)
    println("   -> Total $(length(candidates)) slots scheduled.")
    
    # ★同期成功端末のIDセットを作成（バックオフ計算で使用）★
    synced_terminal_ids = Set([tid for (tid, start_ms) in sync_results if start_ms !== nothing])
    
    # 4. Phase 3: MAC層実行 (Discrete Event Simulation)
    println("Phase 3: Running MAC Layer (DES Mode)...")
    
    # ノイズフロア計算
    noise_bw_hz = params.terminal_bw_mhz * 1e6
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    
    active_tx = TransmissionRecord[]
    finished_tx = TransmissionRecord[]
    
    # CS統計
    cs_stats = Tuple{Bool, Float64}[]
    
    # 端末ごとの次回送信可能時刻 (Duty Cycle用)
    next_available_time = Dict{Int, Float64}()
    for t in terminals
        next_available_time[t.terminal_id] = 0.0
    end
    
    # 再送キュー
    retransmission_queue = TransmissionRecord[]
    
    # メインループ (イベントキューが空になるまで)
    while !isempty(candidates)
        cand = pop!(candidates) # 末尾(最小時刻)を取得 (O(1))
        
        curr_t = cand.time_global_ms
        me = cand.terminal_node
        
        # A. 終わった通信を掃除 (curr_t 時点で終了しているもの)
        # ※ LBT窓 (5ms) を考慮し、curr_t - 5ms 以降に終了したパケットは一旦維持する
        filter!(x -> x.end_ms > curr_t - params.lbt_duration_ms, active_tx)
        
        # B. Duty Cycle チェック
        if curr_t < next_available_time[me.terminal_id]
            continue
        end
        
        # C. キャリアセンス (LBT) - ARIB STD-T108準拠
        is_busy = false
        max_rssi = -Inf
        
        if params.enable_carrier_sense
            # 5ms LBT: 送信予定時刻の5ms前から現在時刻までサンプリング
            cs_start_time = curr_t - params.lbt_duration_ms
            num_samples = Int(ceil(params.lbt_duration_ms / params.lbt_sample_interval_ms))
            
            for i in 0:num_samples
                sample_time = cs_start_time + i * params.lbt_sample_interval_ms
                
                # この時刻での合計受信電力を計算 (線形加算)
                total_power_w = 0.0
                
                for other in active_tx
                    # sample_time時点で送信中か？
                    if other.start_ms <= sample_time && other.end_ms > sample_time
                        # 同一チャネルのみチェック
                        if other.channel != cand.channel
                            continue
                        end
                        
                        # 受信電力計算 for CS
                        dist = sqrt((me.x_m - other.x)^2 + (me.y_m - other.y)^2)
                        pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist, 1.0))
                        
                        # Shadowing: 本来は送信端末(other)の固定シャドウイング値を使用すべき
                        # しかし、TransmissionRecordにshadowing情報がないため、
                        # 現状では簡易的にランダム値を使用
                        # TODO: TransmissionRecordにshadowing_dbを追加するか、
                        # 端末IDから端末情報を参照できるようにする
                        if params.shadowing_enabled
                            pl += randn() * params.shadowing_std_db
                        end
                        
                        rssi = other.tx_power_dbm - pl
                        power_w = 10^(rssi/10) / 1000  # dBm → W
                        total_power_w += power_w
                    end
                end
                
                # ノイズ電力を加算（現実的なCS）
                noise_power_w = 10^(noise_power_dbm/10) / 1000  # dBm → W
                total_power_w += noise_power_w
                
                if total_power_w > 0
                    total_rssi = 10 * log10(total_power_w * 1000)  # W → dBm
                    if total_rssi > max_rssi
                        max_rssi = total_rssi
                    end
                    
                    if total_rssi > params.cs_threshold_dbm
                        is_busy = true
                        break  # 一度でもビジーなら終了
                    end
                end
                if is_busy
                    break  # 一度でもビジーなら終了
                end
            end
        end
        
        # 統計収集
        push!(cs_stats, (is_busy, max_rssi))
        
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
            # シャドウイングは端末固定値を使用（me.shadowing_db）
            rx_gw = params.tx_power_dbm - (pl_gw + me.shadowing_db)
            
            # パケットID生成（グローバルカウンタ）
            packet_id = length(finished_tx) + 1
            
            rec = TransmissionRecord(
                me.terminal_id, curr_t, curr_t+dur, params.tx_power_dbm, 
                me.x_m, me.y_m, "Success", rx_gw, cand.channel,
                0,      # retry_count (初回送信)
                false,  # is_retransmission
                packet_id,  # original_packet_id
                "None",  # failure_reason
                0.0,     # sinr_db
                0.0,     # snr_db
                0,       # num_interferers
                is_busy, # cs_detected (キャリアセンス結果を記録)
                "",      # interferer_ids
                ""       # interferer_distances
            )
            push!(active_tx, rec)
            push!(finished_tx, rec)
            
        else
            # ビジー検出 → バックオフ + チャネル再選択
            
            # バックオフスロット数をランダムに選択（1～3スロット後）
            backoff_slots = rand(1:3)
            
            # ★同期端末と非同期端末で異なる処理★
            if me.terminal_id in synced_terminal_ids
                # 同期端末: スロット境界に揃える
                first_slot_start = sync_results[me.terminal_id]
                slot_len = params.slot_length_ms * me.clock_drift_factor
                
                # 現在時刻から次のスロット境界を計算
                elapsed = curr_t - first_slot_start
                current_slot = floor(elapsed / slot_len)
                next_slot = current_slot + backoff_slots
                new_time_snapped = first_slot_start + next_slot * slot_len
                
                # スロット内ジッタを追加
                guard_band_ms = 10.0
                max_jitter_in_slot = max(0.0, params.slot_length_ms - params.packet_airtime_ms - guard_band_ms)
                jitter_limit = min(params.tx_jitter_max_ms, max_jitter_in_slot)
                jitter = rand() * jitter_limit
                new_time = new_time_snapped + jitter
            else
                # 非同期端末: 通常のバックオフ（ジッタなし、LBT-ALOHAと同じ）
                backoff_ms = 10.0 + rand() * 90.0
                new_time = curr_t + backoff_ms
            end
            
            if new_time < params.simulation_duration_ms
                # 新しいチャネルをランダムに選択（チャネル再選択）
                new_channel = rand(1:params.num_channels)
                
                # 新しいイベントを作成
                new_cand = CandidateSlot(new_time, me, cand.slot_index, new_channel)
                
                # スタックに挿入 (降順を維持)
                idx = searchsortedfirst(candidates, new_cand, by=x->x.time_global_ms, rev=true)
                insert!(candidates, idx, new_cand)
            end
        end
    end
    
    # 5. Phase 4: 結果解析 (衝突判定)
    println("Phase 4: Analyzing Results...")
    
    success_count = 0
    collision_count = 0
    
    if params.collision_model == :overlap
        println("   Mode: Simple Overlap (No Capture)")
        (success_count, collision_count) = detect_collisions_simple_overlap(finished_tx)
    else
        println("   Mode: SINR-based (Capture Effect Enabled)")
        (success_count, collision_count) = detect_collisions_sinr(finished_tx, params.spreading_factor, noise_power_dbm)
    end
    
    # ★ 元のパケットのみでPER計算（再送パケットは含めない）★
    original_packets_list = filter(tx -> !tx.is_retransmission, finished_tx)
    original_total = length(original_packets_list)
    original_collisions = count(tx -> tx.status == "Collision", original_packets_list)
    original_success = count(tx -> tx.status == "Success", original_packets_list)
    per = original_total == 0 ? 0.0 : original_collisions / original_total * 100.0
    
    println("-"^30)
    println("Result Summary (Original Packets Only):")
    println("  Total Packets: $original_total")
    println("  Success:       $original_success")
    println("  Collisions:    $original_collisions")
    if original_total > 0
        println("  PER:           $(round(per, digits=2)) %")
    else
        println("  PER:           N/A")
    end
    
    # === スループット計算 ===
    norm_throughput = (original_success * params.packet_airtime_ms) / (params.simulation_duration_ms * params.num_channels)
    total_bits = original_success * params.lora_payload_bytes * 8
    throughput_bps = total_bits / (params.simulation_duration_ms / 1000.0)

    println("\nThroughput Statistics:")
    println("  Normalized:           $(round(norm_throughput, digits=5))")
    println("  Effective:            $(round(throughput_bps, digits=2)) bps")
    
    # ★ 同期成功/失敗端末の分類統計 ★
    synced_terminal_ids = Set([tid for (tid, start_ms) in sync_results if start_ms !== nothing])
    synced_packets = filter(tx -> tx.terminal_id in synced_terminal_ids, original_packets_list)
    async_packets = filter(tx -> !(tx.terminal_id in synced_terminal_ids), original_packets_list)
    
    if !isempty(synced_packets)
        synced_total = length(synced_packets)
        synced_collisions = count(tx -> tx.status == "Collision", synced_packets)
        synced_per = (synced_collisions / synced_total) * 100.0
        println("  ├─ Synced terminals:  $(synced_total) packets, PER $(round(synced_per, digits=2))%")
    end
    
    if !isempty(async_packets)
        async_total = length(async_packets)
        async_collisions = count(tx -> tx.status == "Collision", async_packets)
        async_per = (async_collisions / async_total) * 100.0
        println("  └─ Async terminals:   $(async_total) packets, PER $(round(async_per, digits=2))%")
    end
    
    # ACK/再送統計
    if params.enable_ack
        println("-"^30)
        println("ACK/Retransmission Statistics:")
        
        # 再送が必要だったパケット数
        retransmission_needed = original_collisions
        
        # 簡易的なPDR推定（再送を実際にシミュレートしていないため推定値）
        # 仮定: 各再送は独立で、同じPERで失敗する
        # 1回目失敗 -> 2回目成功の確率 = PER * (1 - PER)
        # 最大3回再送なので、最終的な配送率を推定
        
        if per > 0 && per < 100
            # 1回目で成功
            p_success_1st = (100.0 - per) / 100.0
            # 1回目失敗、2回目成功
            p_success_2nd = (per / 100.0) * ((100.0 - per) / 100.0)
            # 1回目失敗、2回目失敗、3回目成功
            p_success_3rd = (per / 100.0)^2 * ((100.0 - per) / 100.0)
            # 1回目失敗、2回目失敗、3回目失敗、4回目成功（最大再送回数3）
            p_success_4th = (per / 100.0)^3 * ((100.0 - per) / 100.0)
            
            estimated_pdr = (p_success_1st + p_success_2nd + p_success_3rd + p_success_4th) * 100.0
            estimated_failed = original_total * (1.0 - estimated_pdr / 100.0)
            estimated_retries = original_collisions * (1.0 + per/100.0 + (per/100.0)^2)  # 期待再送回数
            
            println("  Estimated PDR:    $(round(estimated_pdr, digits=2)) % (with max $(params.max_retries) retries)")
            println("  Estimated Retransmissions: $(round(estimated_retries, digits=1))")
            println("  Note: Retransmissions are estimated, not simulated")
        else
            println("  PDR:              $(round(100.0 - per, digits=2)) %")
            println("  Note: No retransmissions needed (PER = $(round(per, digits=2))%)")
        end
    end
    println("-"^30)
    
    # RSSI統計表示
    println("Carrier Sense Statistics:")
    busy_rssis = [x[2] for x in cs_stats if x[1]]
    idle_rssis = [x[2] for x in cs_stats if !x[1] && x[2] > -Inf]
    
    println("  Threshold:        $(params.cs_threshold_dbm) dBm")
    println("  Noise Floor:      $(round(noise_power_dbm, digits=2)) dBm")
    println("  Total Attempts:   $(length(cs_stats))")
    println("  Busy Count:       $(length(busy_rssis))")
    
    if !isempty(busy_rssis)
        println("  Busy RSSI (detected interference):")
        println("    Mean: $(round(mean(busy_rssis), digits=2)) dBm")
        println("    Min:  $(round(minimum(busy_rssis), digits=2)) dBm")
        println("    Max:  $(round(maximum(busy_rssis), digits=2)) dBm")
        
        # ヒストグラム的な表示
        println("    Distribution:")
        for range_start in -120:10:-60
            cnt = count(x -> x >= range_start && x < range_start+10, busy_rssis)
            println("      [$range_start, $(range_start+10)): $cnt")
        end
    end
    
    println("  Idle Count:       $(length(cs_stats) - length(busy_rssis))")
    if !isempty(idle_rssis)
        println("  Idle RSSI (detected but below threshold):")
        println("    Mean: $(round(mean(idle_rssis), digits=2)) dBm")
        println("    Max:  $(round(maximum(idle_rssis), digits=2)) dBm")
    end
    println("-"^30)
    
    # ==========================================
    # パケットロス原因分析
    # ==========================================
    println("\n" * "="^60)
    println("Packet Loss Cause Analysis")
    println("="^60)
    
    # 同期失敗端末数を計算
    total_terminals = length(terminals)
    synced_terminals = count(t -> sync_results[t.terminal_id] !== nothing, terminals)
    sync_failed_terminals = total_terminals - synced_terminals
    
    # 失敗パケットの分類
    failed_packets = filter(tx -> tx.status == "Collision", finished_tx)
    
    # 失敗理由ごとにカウント
    reason_counts = Dict{String, Int}()
    for tx in failed_packets
        reason = tx.failure_reason
        reason_counts[reason] = get(reason_counts, reason, 0) + 1
    end
    
    # 隠れ端末 vs 同時送信の分類
    hidden_terminal_count = 0
    simultaneous_count = 0
    
    for tx in failed_packets
        if tx.failure_reason == "SINR_Insufficient" || tx.failure_reason == "SINR_and_SNR_Insufficient"
            # キャリアセンスで検出できなかった場合は隠れ端末
            if !tx.cs_detected
                hidden_terminal_count += 1
            else
                # 検出したのに衝突した場合は同時送信の可能性
                # より詳細には、干渉パケットとの開始時刻差を確認
                # ここでは簡易的に cs_detected=true なら同時送信とする
                simultaneous_count += 1
            end
        end
    end
    
    # 統計計算
    total_failed = length(failed_packets)
    
    println("\nOverall Statistics:")
    println("  Total Terminals:       $total_terminals")
    println("  Synced Terminals:      $synced_terminals")
    println("  Sync Failed Terminals: $sync_failed_terminals")
    println("\n  Total Packets Sent:    $original_total")
    println("  Success:               $original_success ($(round(original_success/max(original_total,1)*100, digits=1))%)")
    println("  Failed:                $total_failed ($(round(total_failed/max(original_total,1)*100, digits=1))%)")
    
    if total_failed > 0
        println("\nFailure Breakdown:")
        
        # SINR不足の詳細
        sinr_insufficient = get(reason_counts, "SINR_Insufficient", 0)
        snr_insufficient = get(reason_counts, "SNR_Insufficient", 0)
        both_insufficient = get(reason_counts, "SINR_and_SNR_Insufficient", 0)
        
        if hidden_terminal_count > 0
            println("  • Hidden Terminal (SINR):  $hidden_terminal_count ($(round(hidden_terminal_count/total_failed*100, digits=1))% of failures)")
        end
        if simultaneous_count > 0
            println("  • Simultaneous Tx (SINR):  $simultaneous_count ($(round(simultaneous_count/total_failed*100, digits=1))% of failures)")
        end
        if snr_insufficient > 0
            println("  • SNR Insufficient:        $snr_insufficient ($(round(snr_insufficient/total_failed*100, digits=1))% of failures)")
        end
        if both_insufficient > 0
            println("  • SINR+SNR Insufficient:   $both_insufficient ($(round(both_insufficient/total_failed*100, digits=1))% of failures)")
        end
        
        # 詳細統計
        println("\nDetailed Metrics (Failed Packets):")
        
        failed_with_interferers = filter(tx -> tx.num_interferers > 0, failed_packets)
        if !isempty(failed_with_interferers)
            avg_interferers = mean([tx.num_interferers for tx in failed_with_interferers])
            println("  • Avg Interferers:     $(round(avg_interferers, digits=2))")
            
            sinr_values = [tx.sinr_db for tx in failed_with_interferers if tx.sinr_db > -Inf]
            if !isempty(sinr_values)
                println("  • Avg SINR:            $(round(mean(sinr_values), digits=2)) dB")
                println("  • Min SINR:            $(round(minimum(sinr_values), digits=2)) dB")
                println("  • Max SINR:            $(round(maximum(sinr_values), digits=2)) dB")
            end
        end
        
        snr_values = [tx.snr_db for tx in failed_packets if tx.snr_db > -Inf]
        if !isempty(snr_values)
            println("  • Avg SNR:             $(round(mean(snr_values), digits=2)) dB")
            println("  • Min SNR:             $(round(minimum(snr_values), digits=2)) dB")
            println("  • Max SNR:             $(round(maximum(snr_values), digits=2)) dB")
        end
    end
    
    if sync_failed_terminals > 0
        println("\nSynchronization Failures:")
        println("  • Terminals unable to sync: $sync_failed_terminals ($(round(sync_failed_terminals/total_terminals*100, digits=1))%)")
        println("  • These terminals sent 0 packets")
    end
    
    println("="^60)


    # ★ 送信結果のCSV保存 ★
    if !isempty(finished_tx)
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        
        # 1. パケット詳細リスト（失敗原因分析フィールド追加）
        tx_df = DataFrame(
            terminal_id = [r.terminal_id for r in finished_tx],
            start_ms = [r.start_ms for r in finished_tx],
            end_ms = [r.end_ms for r in finished_tx],
            status = [r.status for r in finished_tx],
            failure_reason = [r.failure_reason for r in finished_tx],
            tx_power_dbm = [r.tx_power_dbm for r in finished_tx],
            rx_power_gw_dbm = [r.rx_power_at_gw for r in finished_tx],
            sinr_db = [r.sinr_db for r in finished_tx],
            snr_db = [r.snr_db for r in finished_tx],
            num_interferers = [r.num_interferers for r in finished_tx],
            cs_detected = [r.cs_detected for r in finished_tx],
            interferer_ids = [r.interferer_ids for r in finished_tx],
            interferer_distances = [r.interferer_distances for r in finished_tx],
            x_m = [r.x for r in finished_tx],
            y_m = [r.y for r in finished_tx],
            channel = [r.channel for r in finished_tx]
        )
        # 詳細ログが有効な場合のみ保存
        if params.enable_detailed_logs && params.enable_file_output
            out_tx = joinpath(output_dir, "integrated_tx_log_$(timestamp).csv")
            CSV.write(out_tx, tx_df)
            println("Tx Log saved to $out_tx")
        end
        
        # 2. パケットロス原因分析サマリーCSV
        loss_analysis_df = DataFrame(
            category = String[],
            count = Int[],
            percentage = Float64[]
        )
        
        push!(loss_analysis_df, ("Total_Packets", original_total, 100.0))
        push!(loss_analysis_df, ("Success", original_success, round(original_success/max(original_total,1)*100, digits=2)))
        push!(loss_analysis_df, ("Failed", total_failed, round(total_failed/max(original_total,1)*100, digits=2)))
        push!(loss_analysis_df, ("Hidden_Terminal", hidden_terminal_count, round(hidden_terminal_count/max(original_total,1)*100, digits=2)))
        push!(loss_analysis_df, ("Simultaneous_Tx", simultaneous_count, round(simultaneous_count/max(original_total,1)*100, digits=2)))
        push!(loss_analysis_df, ("SNR_Insufficient", get(reason_counts, "SNR_Insufficient", 0), round(get(reason_counts, "SNR_Insufficient", 0)/max(original_total,1)*100, digits=2)))
        push!(loss_analysis_df, ("SINR_and_SNR_Insufficient", get(reason_counts, "SINR_and_SNR_Insufficient", 0), round(get(reason_counts, "SINR_and_SNR_Insufficient", 0)/max(original_total,1)*100, digits=2)))
        push!(loss_analysis_df, ("Sync_Failed_Terminals", sync_failed_terminals, round(sync_failed_terminals/max(total_terminals,1)*100, digits=2)))
        
        # ファイル出力が有効な場合のみ保存
        if params.enable_file_output
            out_loss = joinpath(output_dir, "packet_loss_analysis_$(timestamp).csv")
            CSV.write(out_loss, loss_analysis_df)
            println("Packet Loss Analysis saved to $out_loss")
        end

        # 3. 空間同期マップ出力 (ステータス別に分離)
        tx_coords = unique(tx_df[:, [:terminal_id, :x_m, :y_m]], :terminal_id)
        sync_map_all = leftjoin(sync_log_df, tx_coords, on=:terminal_id)
        
        if params.enable_file_output
            # 成功端末のみ
            success_map = filter(row -> row.status == "Success", sync_map_all)[:, [:x_m, :y_m, :status]]
            out_success = joinpath(output_dir, "integrated_sync_success_coords_$(timestamp).csv")
            CSV.write(out_success, success_map)
            
            # 失敗端末のみ
            failed_map = filter(row -> row.status == "Failed", sync_map_all)[:, [:x_m, :y_m, :status]]
            out_failed = joinpath(output_dir, "integrated_sync_failed_coords_$(timestamp).csv")
            CSV.write(out_failed, failed_map)
            
            println("Coord Logs saved (Success/Failed separate CSVs)")
        end
        
        # 4. 端末ごとの集計
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
        # ファイル出力が有効な場合のみ保存
        if params.enable_file_output
            out_sum = joinpath(output_dir, "integrated_summary_$(timestamp).csv")
            CSV.write(out_sum, summary_df)
            println("Summary saved to $out_sum")
        end
    end
    
    # (Section 6: 可視化はユーザー要望により削除)
    
    # 結果を返す（評価スクリプト用）
    # ★ コンソール出力と同じ値（original_packetsベース）を返す ★
    # 同期成功率を計算
    total_terminals = params.num_terminals
    synced_count = length(synced_terminal_ids)
    sync_success_rate = (synced_count / total_terminals) * 100.0
    
    # 同期成功/失敗端末のPERを計算
    synced_per = 0.0
    async_per = 0.0
    
    if !isempty(synced_packets)
        synced_total = length(synced_packets)
        synced_collisions = count(tx -> tx.status == "Collision", synced_packets)
        synced_per = (synced_collisions / synced_total) * 100.0
    end
    
    if !isempty(async_packets)
        async_total = length(async_packets)
        async_collisions = count(tx -> tx.status == "Collision", async_packets)
        async_per = (async_collisions / async_total) * 100.0
    end
    
    return Dict(
        "total_packets" => original_total,           # 元のパケット数のみ
        "success" => original_success,               # 元のパケットの成功数
        "collisions" => original_collisions,         # 元のパケットの衝突数
        "per" => per,                                # 元のパケットのPER
        "sync_success_rate" => sync_success_rate,    # 同期成功率
        "synced_per" => synced_per,                  # 同期成功端末のPER
        "async_per" => async_per,                    # 非同期端末のPER
        # === パケットロス原因の内訳 ===
        "hidden_terminal" => hidden_terminal_count,
        "simultaneous_tx" => simultaneous_count,
        "snr_insufficient" => get(reason_counts, "SNR_Insufficient", 0),
        "sinr_and_snr_insufficient" => get(reason_counts, "SINR_and_SNR_Insufficient", 0),
        "sync_failed_terminals" => sync_failed_terminals,
        "throughput_bps" => throughput_bps,
        "norm_throughput" => norm_throughput
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
