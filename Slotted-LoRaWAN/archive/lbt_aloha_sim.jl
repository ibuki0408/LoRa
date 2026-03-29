# ============================================================
# LBT-ALOHA Simulation
# Unslotted ALOHA with Listen Before Talk
# ============================================================

# 必要なモジュールの読み込み
include("../src/modules/lora_airtime.jl")  # ToA計算
include("../src/modules/collision_detection.jl") # 衝突判定
include("../src/modules/terminal_deployment.jl") # 端末配置・パスロス

# using .LoRaAirtime  <- Remove this (functions are top-level)
using Random
using Statistics
using Dates
using Printf

# ========================================
# パラメータ構造体
# ========================================
struct CSMAParameters
    # === 物理層 ===
    tx_power_dbm::Float64
    noise_figure_db::Float64
    
    # === 環境 ===
    num_terminals::Int
    area_size_m::Float64
    
    # === LoRa固有 ===
    sf::Int
    payload_bytes::Int
    num_channels::Int  # マルチチャネル数（1-16、AS923準拠）
    
    # === シミュレーション制御 ===
    simulation_duration_ms::Float64
    duty_cycle::Float64
    enable_carrier_sense::Bool  # キャリアセンス有効/無効
    mean_event_interval_ms::Float64 # ポアソン分布の平均送信間隔
    max_startup_delay_ms::Float64   # 起動遅延の上限（Integratedと同じ）
    
    # === 環境モデル ===
    pass_loss_exp::Float64   # パスロス指数
    shadowing_std_db::Float64 # シャドウイング標準偏差
    reference_path_loss_db::Float64 # 基準線路損失 (LoRa帯)
    
    # === LBT (Listen Before Talk) ===
    lbt_duration_ms::Float64             # LBT期間（ARIB: 5ms）
    lbt_sample_interval_ms::Float64      # サンプリング間隔
end

function create_csma_params()
    # IntegratedSim と同じ設定を使用
    sf = 8
    payload_bytes = 10
    data_freq_ghz = 0.92
    
    # FSPL (Free Space Path Loss) at 1m for 920MHz
    # FSPL(dB) = 20log10(f) + 20log10(d) - 147.55
    # d=1m -> 20log10(f) - 147.55
    ref_pl_db = 20*log10(data_freq_ghz*1e9) - 147.55
    
    return CSMAParameters(
        # PHY
        13.0,     # tx_power_dbm
        6.0,      # noise_figure_db
        
        # Environment
        400,       # num_terminals (デバッグ用に削減)
        2000.0,    # area_size_m
        
        # LoRa
        sf,
        payload_bytes,
        8,        # num_channels (AS923 Japan typical (1-16))
        
        # Simulation
        3600000.0, # simulation_duration_ms (10分)
        0.01,     # duty_cycle (1%)
        true,     # enable_carrier_sense (true: CSMA/CA, false: Pure ALOHA)
        30000.0,  # mean_event_interval_ms (平均60秒間隔 = 少し余裕を持たせる)
        60000.0,  # max_startup_delay_ms (最大60秒、Integratedと同じ)
        
        # Path Loss Model
        2.7,      # pass_loss_exp
        8.0,      # shadowing_std_db
        ref_pl_db, # reference_path_loss_db
        
        # === LBT (Listen Before Talk) ===
        5.0,      # lbt_duration_ms (ARIB STD-T108: 5ms)
        1.0       # lbt_sample_interval_ms (1ms間隔でサンプリング)
    )
end

# 送信記録構造体 (CollisionDetectionと互換)
# 送信記録構造体 (CollisionDetectionと互換)
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
    
    # 衝突詳細（CollisionDetectionモジュールで必要）
    sinr_db::Float64
    snr_db::Float64
    num_interferers::Int
    failure_reason::String
    cs_detected::Bool
    interferer_ids::String
    interferer_distances::String
    
    # コンストラクタ
    TransmissionRecord(id, s, e, txp, x, y, st, rx, ch) = 
        new(id, s, e, txp, x, y, st, rx, ch, -Inf, -Inf, 0, "None", false, "", "")
end

# ========================================
# メインシミュレーション
# ========================================
# イベント構造体
struct TransmissionEvent
    time_ms::Float64
    terminal
    attempt::Int  # 試行回数（バックオフ用）
    channel::Int  # 送信に使用するターゲットチャネル
end

"""
    run_csma_ca_with_params(params::CSMAParameters)

パラメータを受け取ってLBT-ALOHAシミュレーションを実行する（評価スクリプト用）
"""
function run_csma_ca_with_params(params::CSMAParameters)
    println("="^60)
    println("   LBT-ALOHA Simulation")
    println("   Unslotted ALOHA with Listen Before Talk")
    println("="^60)
    
    # ToA計算 (Corrected Usage)
    lora_params = create_lora_params(params.sf, params.payload_bytes)
    airtime_ms = calculate_lora_airtime(lora_params)
    
    println("LoRa Settings:")
    println("  SF: $(params.sf)")
    println("  Payload: $(params.payload_bytes) bytes")
    println("  ToA: $(round(airtime_ms, digits=2)) ms")
    println("  Duty Cycle: $(params.duty_cycle * 100)%")
    println("  Terminals: $(params.num_terminals)")
    println()
    
    # 1. 端末配置 (terminal_deployment.jl を使用)
    println("Phase 1: Initializing Settings...")
    
    # 周波数（Hz）
    freq_hz = 920.0 * 1e6 
    
    # 端末配置パラメータ作成
    dep_p = TerminalDeploymentParameters(
        "random_fixed",      # 固定数ランダム配置
        0.0,                 # lambda (不要)
        params.num_terminals,
        params.area_size_m,
        10.0,                # min_distance_m
        params.area_size_m/2, # max_distance_m (半径)
        freq_hz,
        params.pass_loss_exp,
        1.0,                 # reference_distance_m
        params.reference_path_loss_db,
        0.0, 0.0             # Gateway位置 (0,0)
    )
    
    # 端末配置実行
    # Note: CSMAParameters removed shadowing_enabled, so we assume true if std_db > 0
    shadowing_on = params.shadowing_std_db > 0.0
    terminals = deploy_terminals(dep_p, params.shadowing_std_db, shadowing_on, params.tx_power_dbm)
    
    # 2. CSMA/CA with Event-Driven Simulation
    println("Phase 2: Running LBT-ALOHA (Event-Driven)...")
    
    # 送信周期 (Debug用表示)
    println("  Mean Tx Interval: $(params.mean_event_interval_ms) ms (Poisson)")
    
    # CS閾値
    cs_threshold_dbm = -80.0
    max_retries = 5
    
    # イベントキュー初期化
    events = TransmissionEvent[]
    for t in terminals
        # 起動遅延: 0 ~ max_startup_delay_ms の一様分布
        startup_delay = rand() * params.max_startup_delay_ms
        # 初回送信時刻 = 起動遅延 + ポアソン間隔
        initial_interval = -params.mean_event_interval_ms * log(rand()) * t.clock_drift_factor
        initial_time = startup_delay + initial_interval
        # 初回チャネル選択
        initial_channel = rand(1:params.num_channels)
        push!(events, TransmissionEvent(initial_time, t, 0, initial_channel))
    end
    
    # 降順ソート（末尾からpop）
    sort!(events, by=x->x.time_ms, rev=true)
    
    println("  Channels: $(params.num_channels)")
    
    # 送信記録
    active_tx = TransmissionRecord[]
    finished_tx = TransmissionRecord[]
    
    # CS統計
    cs_stats = Tuple{Bool, Float64}[]
    
    # Duty Cycle管理
    next_available_time = Dict{Int, Float64}()
    for t in terminals
        next_available_time[t.terminal_id] = 0.0
    end
    
    # メインループ
    while !isempty(events)
        event = pop!(events)
        curr_t = event.time_ms
        t = event.terminal
        channel = event.channel  # 事前に割り当てられたチャネルを使用
        
        # 終了した送信を削除
        # ※ LBT窓 (5ms) を考慮し、curr_t - 5ms 以降に終了したパケットは一旦維持する
        filter!(x -> x.end_ms > curr_t - params.lbt_duration_ms, active_tx)
        
        # Duty Cycle チェック
        if curr_t < next_available_time[t.terminal_id]
            # クロックドリフト等でわずかに早くイベントが来た場合、
            # 破棄してはいけない（チェーンが切れる）。
            # 送信可能時刻まで待機して再スケジュールする。
            delayed_time = next_available_time[t.terminal_id]
            new_event = TransmissionEvent(delayed_time, t, event.attempt, channel)
            
            # 再挿入
            idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
            insert!(events, idx, new_event)
            continue
        end
        
        # ノイズフロア計算（CS用）
        noise_bw_hz = 125000.0  # 125 kHz
        noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
        
        # Carrier Sense - ARIB STD-T108準拠 (5ms LBT)
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
                    # 自分自身はスキップ
                    if other.terminal_id == t.terminal_id
                        continue
                    end
                    
                    # sample_time時点で送信中か？
                    if other.start_ms <= sample_time && other.end_ms > sample_time
                        # 同一チャネルのみチェック
                        if other.channel != channel
                            continue
                        end
                        
                        dist = sqrt((t.x_m - other.x)^2 + (t.y_m - other.y)^2)
                        pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist, 1.0))
                        
                        # Shadowing: 送信端末(other)の固定シャドウイング値を使用
                        # 注: 厳密には送受信ペアごとに異なるシャドウイングだが、
                        # 簡易的に送信端末の値を使用
                        if params.shadowing_std_db > 0.0
                            # otherはTransmissionRecordなので、元の端末情報から取得する必要がある
                            # しかし、ここではTransmissionRecordにshadowing情報がないため、
                            # 簡易的に距離ベースのパスロスのみを使用
                            # TODO: TransmissionRecordにshadowing_dbを追加するか、
                            # 端末IDから端末情報を参照できるようにする
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
                    
                    if total_rssi > cs_threshold_dbm
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

        if !is_busy
            # 送信成功
            end_time = curr_t + airtime_ms
            
            # Duty Cycle更新
            off_period = airtime_ms * (1.0 / params.duty_cycle - 1.0)
            next_available_time[t.terminal_id] = end_time + off_period
            
            # GW受信電力
            # t.distance_m is available in TerminalInfo. Apply shadowing if needed.
            # シャドウイングは端末固定値を使用（t.shadowing_db）
            pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(t.distance_m, 1.0))
            rx_gw = params.tx_power_dbm - (pl_gw + t.shadowing_db)
            
            # 事前に選んだチャネルを使用
            rec = TransmissionRecord(t.terminal_id, curr_t, end_time, params.tx_power_dbm, t.x_m, t.y_m, "Success", rx_gw, channel)
            push!(active_tx, rec)
            push!(finished_tx, rec)
            
            # 次の送信イベントをスケジュール (ポアソン分布)
            # interval = -mean * ln(rand)
            interval = -params.mean_event_interval_ms * log(rand()) * t.clock_drift_factor
            next_time = curr_t + interval
            
            # DC制約の適用 (Clamping)
            # 次回送信可能時刻 (next_available_time[t.id]) より早い場合は待たせる
            # ※ ただし、このイベントループ内では `next_available_time` はすでに `end_time + off_period` に更新済み
            # 次のイベント処理時に `if curr_t < next_available_time` でキャッチされるので、
            # ここでは単純に next_time を登録しておけばよい。
            # (ただし、極端に早い時間を登録すると無駄なループが増えるので、ここでmaxをとっても良い)
            
            if next_time < params.simulation_duration_ms
                # 次回のチャネルを新たに選択
                next_channel = rand(1:params.num_channels)
                new_event = TransmissionEvent(next_time, t, 0, next_channel)
                idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
                insert!(events, idx, new_event)
            end
        else
            # ビジー検出 → バックオフ + チャネル再選択（Integratedと同じ戦略）
            backoff_ms = 10.0 + rand() * 90.0
            new_time = curr_t + backoff_ms
            
            if new_time < params.simulation_duration_ms
                # バックオフ時は新しいチャネルをランダムに選択（チャネル再選択）
                new_channel = rand(1:params.num_channels)
                new_event = TransmissionEvent(new_time, t, event.attempt + 1, new_channel)
                idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
                insert!(events, idx, new_event)
            end
        end
    end
    
    # 3. 衝突判定 (SINRベース)
    println("Phase 3: Analysing Collisions (SINR-based)...")
    
    # ノイズフロア計算
    # BW = 125kHz
    noise_bw_hz = 125000.0
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    
    (success_count, collision_count) = detect_collisions_sinr(finished_tx, params.sf, noise_power_dbm)
    
    total_packets = length(finished_tx)
    per = total_packets == 0 ? 0.0 : collision_count / total_packets * 100
    
    println("------------------------------")
    println("Result Summary (LBT-ALOHA):")
    println("  Total Packets: $total_packets")
    println("  Success:       $success_count")
    println("  Collisions:    $collision_count")
    println("  PER:           $(round(per, digits=2)) %")
    println("-"^30)
    
    # === スループット計算 ===
    norm_throughput = (success_count * airtime_ms) / (params.simulation_duration_ms * params.num_channels)
    total_bits = success_count * params.payload_bytes * 8
    throughput_bps = total_bits / (params.simulation_duration_ms / 1000.0)

    println("Throughput Statistics:")
    println("  Normalized:           $(round(norm_throughput, digits=5))")
    println("  Effective:            $(round(throughput_bps, digits=2)) bps")
    println("-"^30)
    
    # RSSI統計表示
    println("Carrier Sense Statistics:")
    busy_rssis = [x[2] for x in cs_stats if x[1]]
    idle_rssis = [x[2] for x in cs_stats if !x[1] && x[2] > -Inf]
    
    println("  Threshold:        $cs_threshold_dbm dBm")
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
    else
        println("  Idle RSSI: None detected (Clean Channel)")
    end
    println("-"^30)
    
    # 結果保存
    output_dir = "result_csma_ca"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    summary_file = joinpath(output_dir, "csma_ca_summary_$(timestamp).csv")
    open(summary_file, "w") do io
        println(io, "Metric,Value")
        println(io, "TotalPackets,Success,Collisions,PER")
        println(io, "$total_packets,$success_count,$collision_count,$(round(per, digits=2))")
    end
    println("Summary saved to $summary_file")
    
    # 結果を返す（評価スクリプト用）
    return Dict(
        "total_packets" => total_packets,
        "success" => success_count,
        "collisions" => collision_count,
        "per" => per,
        "throughput_bps" => throughput_bps,
        "norm_throughput" => norm_throughput
    )
end

"""
    run_csma_ca()

デフォルトパラメータでLBT-ALOHAシミュレーションを実行する（既存の動作を維持）
"""
function run_csma_ca()
    params = create_csma_params()
    return run_csma_ca_with_params(params)
end

# スクリプトとして実行された場合のみ実行
if abspath(PROGRAM_FILE) == @__FILE__
    run_csma_ca()
end