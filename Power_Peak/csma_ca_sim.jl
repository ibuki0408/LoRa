# ============================================================
# CSMA/CA Simulation for Comparison with LBT
# ============================================================

# 必要なモジュールの読み込み
include("modules/lora_airtime.jl")  # ToA計算
include("modules/collision_detection.jl") # 衝突判定

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
    
    # === 環境モデル ===
    pass_loss_exp::Float64   # パスロス指数
    shadowing_std_db::Float64 # シャドウイング標準偏差
    reference_path_loss_db::Float64 # 基準線路損失 (LoRa帯)
end

function create_csma_params()
    # IntegratedSim と同じ設定を使用
    sf = 10
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
        50,       # num_terminals (デバッグ用に削減)
        500.0,    # area_size_m
        
        # LoRa
        sf,
        payload_bytes,
        8,        # num_channels (AS923 Japan typical (1-16))
        
        # Simulation
        600000.0, # simulation_duration_ms (10分)
        0.01,     # duty_cycle (1%)
        true,     # enable_carrier_sense (true: CSMA/CA, false: Pure ALOHA)
        60000.0,  # mean_event_interval_ms (平均60秒間隔 = 少し余裕を持たせる)
        
        # Path Loss Model
        2.7,      # pass_loss_exp (IntegratedSimと同じ値: 2.5)
        8.0,      # shadowing_std_db (IntegratedSimと同じ0.0に設定)
        ref_pl_db # reference_path_loss_db
    )
end

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
end

# ========================================
# メインシミュレーション
# ========================================
# イベント構造体
struct TransmissionEvent
    time_ms::Float64
    terminal
    attempt::Int  # 試行回数（バックオフ用）
end

function run_csma_ca()
    println("="^60)
    println("   CSMA/CA Simulation")
    println("   Comparison Target for LBT")
    println("="^60)
    
    params = create_csma_params()
    
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
    
    # 1. 端末配置
    println("Phase 1: Initializing Settings...")
    terminals = []
    # GWは原点 (0,0)
    
    # ランダム配置（シードなし - 毎回異なる結果）
    # Random.seed!(1234)  # コメントアウト: 毎回異なる結果を得るため
    for i in 1:params.num_terminals
        x = (rand() - 0.5) * params.area_size_m
        y = (rand() - 0.5) * params.area_size_m
        
        # クロックドリフト初期化
        drift_ppm = rand(-20.0:0.1:20.0)
        drift_factor = 1.0 + drift_ppm / 1e6
        
        push!(terminals, (id=i, x=x, y=y, drift_ppm=drift_ppm, drift_factor=drift_factor))
    end
    
    # 2. CSMA/CA with Event-Driven Simulation
    println("Phase 2: Running CSMA/CA (Event-Driven)...")
    
    # 送信周期 (Debug用表示)
    println("  Mean Tx Interval: $(params.mean_event_interval_ms) ms (Poisson)")
    
    # CS閾値
    cs_threshold_dbm = -80.0
    max_retries = 5
    
    # イベントキュー初期化
    events = TransmissionEvent[]
    for t in terminals
        # 初回はランダム (0 ～ 平均間隔 * 2 くらいで分散)
        initial_time = rand() * params.mean_event_interval_ms * 2
        push!(events, TransmissionEvent(initial_time, t, 0))
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
        next_available_time[t.id] = 0.0
    end
    
    # メインループ
    while !isempty(events)
        event = pop!(events)
        curr_t = event.time_ms
        t = event.terminal
        
        # 終了した送信を削除
        filter!(x -> x.end_ms > curr_t, active_tx)
        
        # Duty Cycle チェック
        if curr_t < next_available_time[t.id]
            # クロックドリフト等でわずかに早くイベントが来た場合、
            # 破棄してはいけない（チェーンが切れる）。
            # 送信可能時刻まで待機して再スケジュールする。
            delayed_time = next_available_time[t.id]
            new_event = TransmissionEvent(delayed_time, t, event.attempt)
            
            # 再挿入
            idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
            insert!(events, idx, new_event)
            continue
        end
        
        # Carrier Sense
        is_busy = false
        max_rssi = -Inf
        
        if params.enable_carrier_sense
            for other in active_tx
                # 自分自身はスキップ
            if other.terminal_id == t.id
                continue
            end
            
            # 同一チャネルのみチェック（異なるチャネルは干渉しない）
            if other.channel != 1  # CSMA/CAは現在1チャネル固定
                continue
            end
            
            dist = sqrt((t.x - other.x)^2 + (t.y - other.y)^2)
            pl = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist, 1.0))  # log10(0)回避
            rssi = other.tx_power_dbm - pl
            
                if rssi > max_rssi
                    max_rssi = rssi
                end
                
                if rssi > cs_threshold_dbm
                    is_busy = true
                    # break # 他の干渉源も考慮するため、breakしない
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
            next_available_time[t.id] = end_time + off_period
            
            # GW受信電力
            dist_gw = sqrt(t.x^2 + t.y^2)
            pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(max(dist_gw, 1.0))
            rx_gw = params.tx_power_dbm - pl_gw
            
            # ランダムチャネル選択（AS923準拠）
            channel = rand(1:params.num_channels)
            
            rec = TransmissionRecord(t.id, curr_t, end_time, params.tx_power_dbm, t.x, t.y, "Success", rx_gw, channel)
            push!(active_tx, rec)
            push!(finished_tx, rec)
            
            # 次の送信イベントをスケジュール (ポアソン分布)
            # interval = -mean * ln(rand)
            interval = -params.mean_event_interval_ms * log(rand()) * t.drift_factor
            next_time = curr_t + interval
            
            # DC制約の適用 (Clamping)
            # 次回送信可能時刻 (next_available_time[t.id]) より早い場合は待たせる
            # ※ ただし、このイベントループ内では `next_available_time` はすでに `end_time + off_period` に更新済み
            # 次のイベント処理時に `if curr_t < next_available_time` でキャッチされるので、
            # ここでは単純に next_time を登録しておけばよい。
            # (ただし、極端に早い時間を登録すると無駄なループが増えるので、ここでmaxをとっても良い)
            
            if next_time < params.simulation_duration_ms
                new_event = TransmissionEvent(next_time, t, 0)
                idx = searchsortedfirst(events, new_event, by=x->x.time_ms, rev=true)
                insert!(events, idx, new_event)
            end
        else
            # ビジー検出 → バックオフ（LBTと同じ：リトライ制限なし）
            backoff_ms = 10.0 + rand() * 90.0
            new_time = curr_t + backoff_ms
            
            if new_time < params.simulation_duration_ms
                new_event = TransmissionEvent(new_time, t, event.attempt + 1)
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
    
    println("-"^30)
    println("Result Summary (CSMA/CA):")
    println("  Total Packets: $total_packets")
    println("  Success:       $success_count")
    println("  Collisions:    $collision_count")
    println("  PER:           $(round(per, digits=2)) %")
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
end

run_csma_ca()