# ============================================================
# Pure ALOHA Simulation (Unslotted) for Comparison
# ============================================================

# 必要なモジュールの読み込み
include("modules/lora_airtime.jl")  # ToA計算
include("modules/collision_detection.jl") # 衝突判定

# using .LoRaAirtime  <- Remove this (functions are top-level)
using Random
using Dates
using Printf

# ========================================
# パラメータ構造体
# ========================================
struct AlohaParameters
    # === 物理層 ===
    tx_power_dbm::Float64
    noise_figure_db::Float64
    
    # === 環境 ===
    num_terminals::Int
    area_size_m::Float64
    
    # === LoRa固有 ===
    sf::Int
    payload_bytes::Int
    
    # === シミュレーション制御 ===
    simulation_duration_ms::Float64
    duty_cycle::Float64
    
    # === 環境モデル ===
    pass_loss_exp::Float64   # パスロス指数
    shadowing_std_db::Float64 # シャドウイング標準偏差
    reference_path_loss_db::Float64 # 基準線路損失 (LoRa帯)
end

function create_aloha_params()
    # IntegratedSim と同じ設定を使用
    sf = 10
    payload_bytes = 10
    data_freq_ghz = 0.92
    
    # FSPL (Free Space Path Loss) at 1m for 920MHz
    # FSPL(dB) = 20log10(f) + 20log10(d) - 147.55
    # d=1m -> 20log10(f) - 147.55
    ref_pl_db = 20*log10(data_freq_ghz*1e9) - 147.55
    
    return AlohaParameters(
        # PHY
        13.0,     # tx_power_dbm
        6.0,      # noise_figure_db
        
        # Environment
        50,       # num_terminals (IntegratedSimに合わせて変更)
        500.0,    # area_size_m
        
        # LoRa
        sf,
        payload_bytes,
        
        # Simulation
        600000.0, # simulation_duration_ms (10分)
        0.01,     # duty_cycle (1%)
        
        # Path Loss Model
        2.5,      # pass_loss_exp (IntegratedSimと同じ値: 2.5)
        0.0,      # shadowing_std_db (IntegratedSimと同じ0.0に設定)
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
end

# ========================================
# メインシミュレーション
# ========================================
function run_pure_aloha()
    println("="^60)
    println("   Pure ALOHA Simulation (Unslotted)")
    println("   Comparison Target for Integrated Sim")
    println("="^60)
    
    params = create_aloha_params()
    
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
    
    # ランダム配置（シード固定）
    Random.seed!(1234)
    for i in 1:params.num_terminals
        x = (rand() - 0.5) * params.area_size_m
        y = (rand() - 0.5) * params.area_size_m
        push!(terminals, (id=i, x=x, y=y))
    end
    
    # 2. トラフィック生成 (Pure ALOHA)
    # 各端末は DC制約を守りつつ、ランダムな初期オフセットで送信を開始
    
    println("Phase 2: Generating Traffic...")
    finished_tx = TransmissionRecord[]
    
    # 送信周期 (ToA + OffPeriod)
    # OffPeriod = ToA * (1/DC - 1)
    tx_period_ms = airtime_ms / params.duty_cycle
    println("  Tx Period per Terminal: $(round(tx_period_ms, digits=1)) ms")
    
    for t in terminals
        # ランダムな初期オフセット (0 ～ tx_period_ms)
        # これにより全端末の送信タイミングが非同期になる
        current_time = rand() * tx_period_ms
        
        while current_time < params.simulation_duration_ms
            # 送信終了時刻
            end_time = current_time + airtime_ms
            
            if end_time > params.simulation_duration_ms
                break
            end
            
            # GWでの受信電力計算
            dist_gw = sqrt(t.x^2 + t.y^2)
            # FSPL + Shadowing (Corrected Logic via params)
            pl_gw = params.reference_path_loss_db + 10 * params.pass_loss_exp * log10(dist_gw)
            
            if params.shadowing_std_db > 0
                pl_gw += randn() * params.shadowing_std_db
            end
            rx_gw = params.tx_power_dbm - pl_gw
            
            # レコード作成
            rec = TransmissionRecord(
                t.id, 
                current_time, 
                end_time, 
                params.tx_power_dbm, 
                t.x, 
                t.y, 
                "Pending", 
                rx_gw
            )
            push!(finished_tx, rec)
            
            # 次の送信時刻 (厳密なDC周期)
            current_time += tx_period_ms
        end
    end
    
    # 時間順にソート
    sort!(finished_tx, by = x -> x.start_ms)
    
    # 3. 衝突判定 (SINRベース)
    println("Phase 3: Analysing Collisions (SINR-based)...")
    
    # ノイズフロア計算
    # BW = 125kHz
    noise_bw_hz = 125000.0
    noise_power_dbm = -174 + 10*log10(noise_bw_hz) + params.noise_figure_db
    
    (success, collisions) = detect_collisions_sinr(finished_tx, params.sf, noise_power_dbm)
    
    total_sent = length(finished_tx)
    
    println("-"^30)
    println("Result Summary (Pure ALOHA):")
    println("  Total Packets: $total_sent")
    println("  Success:       $success")
    println("  Collisions:    $collisions")
    
    if total_sent > 0
        per = collisions / total_sent * 100
        println("  PER:           $(round(per, digits=2)) %")
    else
        println("  PER:           N/A")
    end
    println("-"^30)
    
    # 結果保存
    output_dir = "result_pure_aloha"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    summary_file = joinpath(output_dir, "pure_aloha_summary_$(timestamp).csv")
    open(summary_file, "w") do io
        println(io, "Metric,Value")
        println(io, "TotalPackets,$total_sent")
        println(io, "Success,$success")
        println(io, "Collisions,$collisions")
        println(io, "PER_Percent,$(total_sent > 0 ? collisions/total_sent*100 : 0)")
    end
    println("Summary saved to $summary_file")
end

run_pure_aloha()
