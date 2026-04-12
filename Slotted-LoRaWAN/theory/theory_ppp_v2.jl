# Theoretical LoRaWAN Analysis Script (Poisson-Discrete Model)
# Features:
# 1. Poisson probability P(k) for simultaneous collisions.
# 2. Capture effect probability for k interferers.
# 4. Success probability including retries (Hidden Node Model).
# 5. Distinguishes between one-shot failure and collision-only PER.

using Distributions
using Random
using Statistics
using Printf
using CSV
using DataFrames
using Dates

# Include LoRa Airtime module
include(joinpath(@__DIR__, "..", "src", "modules", "lora_airtime.jl"))
using .LoraAirtime

# ==========================================
# Parameters (Match theory.jl)
# ==========================================
R_MIN = 10.0
R_MAX = 500.0          # エリア半径 (m) (Match Prop.jl default)
TX_POWER_DBM = 13.0
PATH_LOSS_EXP = 2.7
REF_PATH_LOSS_DB = 31.65 # 20*log10(0.92e9) - 147.55
SHADOWING_STD_DB = 0.0  # Prop.jlのパラメータと合わせる

# Capture & LBT Thresholds
CS_THRESHOLD_DBM = -80.0
REQUIRED_SIR_DB = 6.0
REQUIRED_SNR_DB = -10.0
NOISE_DBM = -174 + 10*log10(125e3) + 6.0 # -117 dBm

# Traffic Parameters
SF = 8
PAYLOAD_BYTES = 10
MEAN_INTERVAL_SEC = 30.0
NUM_CHANNELS = 8
N_RANGE = 100:100:500 
SYNC_RATES = [1.0]

# --- 物理レイヤ設定 (シミュレータと同期) ---
ENABLE_CARRIER_SENSE = true
ENABLE_CAPTURE_EFFECT = true

# ==========================================
# Derived Constants
# ==========================================
# ==========================================
# Derived Constants (Effective Sensing Area)
# ==========================================
# センシング距離 R_cs の計算 (シャドウイングを考慮した積分による実効面積)
# calculate_effective_cs_area関数でシミュレーション実行前に1度だけ計算する

# ==========================================
# Analytical Functions
# ==========================================

"""
    P_capture_with_k_interferers(k, alpha, sir_db)

Probability that one signal survives when k other signals are present.
Simplified for 1D/2D random placement in distance range.
For alpha=2.7, SIR=6dB:
  k=0: 1.0 (Ignoring SNR)
  k=1: ~0.18 (Strongest of 2 survives)
  k>=2: approx 0.18 / k^1.2
"""
function P_capture_with_k_interferers(k, alpha, sir_db, enable_capture)
    if k == 0 return 1.0 end
    if !enable_capture return 0.0 end # キャプチャなし：1つでも衝突があれば失敗
    return 0.18 / (k^1.3)
end

"""
    calculate_effective_cs_area(r_max, tx_pow, cs_th, pl_ref, alpha, sigma)

シャドウイングを考慮した、実効的なキャリアセンス面積（Effective Sensing Area）を計算します。
距離 r において、平均受信電力にシャドウイング S ~ N(0, sigma^2) が加わった結果、
受信電力が cs_th を超える確率を積分します。
"""
function calculate_effective_cs_area(r_max, tx_pow, cs_th, pl_ref, alpha, sigma)
    norm_dist = Normal(0, 1)
    area = 0.0
    dr = 1.0 # 1m刻みで数値積分
    
    for r in 1.0:dr:r_max
        # 距離 r での平均受信電力
        prx_mean = tx_pow - (pl_ref + 10 * alpha * log10(r))
        
        # 閾値を超えるために必要なシャドウイング量 (dB)
        req_shadowing = cs_th - prx_mean
        
        # 必要なシャドウイング量が得られる確率 (右側確率)
        p_detect = ccdf(norm_dist, req_shadowing / sigma)
        
        # 面積要素 (2πr dr) に確率を掛けて足し合わせる
        area += 2 * pi * r * dr * p_detect
    end
    
    return area
end

"""
    calculate_average_area_ratio(r_max, sigma, cs_th, tx_pow, pl_ref, alpha)

エリア半径 r_max の円内に一様に分布する2つの点 A, B について、
点 A から見た実効的なセンス面積の平均的な割合 (AreaRatio) を、数値積分（エッジ効果補正）によって計算します。
"""
function calculate_average_area_ratio(R_max, sigma, cs_th, tx_pow, pl_ref, alpha)
    norm_dist = Normal(0, 1)
    
    total_avg_prob = 0.0
    steps = 50 # 精度を上げるため少し増やす
    
    # 規格化された積分: 2 * 2 * (1/pi * pi) ... ではない
    # \int_0^1 2r dr \int_0^1 2rho drho \int_0^pi (1/pi) dtheta = 1
    
    for r_f in range(0, 1, length=steps)
        r = r_f * R_max
        weight_r = 2 * r_f * (1.0/steps)
        
        for rho_f in range(0, 1, length=steps)
            rho = rho_f * R_max
            weight_rho = 2 * rho_f * (1.0/steps)
            
            for theta in range(0, pi, length=steps)
                # A(r, 0) と B(rho, theta) の間の距離
                dist = sqrt(r^2 + rho^2 - 2*r*rho*cos(theta))
                
                # 物理モデルに基づく検知確率
                prx_mean = tx_pow - (pl_ref + 10 * alpha * log10(max(dist, 1.0)))
                p_det = ccdf(norm_dist, (cs_th - prx_mean) / sigma)
                
                # 重みを掛けて加算 (theta方向はステップ数で割る)
                total_avg_prob += p_det * weight_r * weight_rho * (1.0/steps)
            end
        end
    end
    
    return total_avg_prob
end

function calculate_lbt_tx_prob(G_load, S_eff_cs, R_max, is_slotted=false, n_total=0, t_slot=0.1, t_interval=30.0, n_channels=8)
    S_total = pi * (R_max^2)
    area_ratio = S_eff_cs / S_total
    
    if is_slotted
        # Slotted CS Model: Based on 'number of contenders per slot'
        # M is the expected number of terminals in the sensing range targeting the same slot
        M = (n_total * area_ratio * t_slot / t_interval) / n_channels
        
        # Probability that a transmission is blocked (contention-based)
        # 1 - (expected winners / expected attempts) = 1 - ( (1-e^-M)/M )
        if M > 0
            p_block = 1.0 - (1.0 - exp(-M)) / M
            return 1.0 - p_block
        else
            return 1.0
        end
    else
        # Continuous/Async Model: Based on 'average airtime occupancy'
        return exp(-G_load * area_ratio)
    end
end

# ==========================================
# Main Processing
# ==========================================

println("="^60)
println("  Refined Theory (Poisson-Discrete Model)")
println("="^60)

lora_params = create_lora_params(SF, PAYLOAD_BYTES)
airtime_ms = calculate_lora_airtime(lora_params)
airtime_sec = airtime_ms / 1000.0

output_dir = joinpath(@__DIR__, "..", "results", "analysis")
mkpath(output_dir)
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

for sync_rate in SYNC_RATES
    println("\nAnalyzing Sync Rate: $(sync_rate * 100)%")
    df = DataFrame(num_terminals=Float64[], one_shot_per=Float64[], attempt_per=Float64[], retry_per=Float64[], cs_block_prob=Float64[])
    
    for N in N_RANGE
        # 1. Traffic Load Parameters
        SLOT_LENGTH_SEC = 0.200 # 200ms (Matches Prop.jl slot200ms)
        
        # Baseline Airtime G (for Async and LBT)
        G = (N / MEAN_INTERVAL_SEC * airtime_sec) / NUM_CHANNELS
        G_sync = G * sync_rate
        G_async = G * (1.0 - sync_rate)
        
        # MAC-level traffic density (for slotted collisions)
        # For synchronized packets, the vulnerability window is exactly the slot length
        G_sync_slot = (N * sync_rate / MEAN_INTERVAL_SEC * 100e-3) / NUM_CHANNELS
        
        # 2. LBT Suppression (Updated with Edge Effect Correction)
        # S_EFF_CS = ENABLE_CARRIER_SENSE ? calculate_effective_cs_area(R_MAX, TX_POWER_DBM, CS_THRESHOLD_DBM, REF_PATH_LOSS_DB, PATH_LOSS_EXP, SHADOWING_STD_DB) : 0.0
        # p_lbt = calculate_lbt_tx_prob(G, S_EFF_CS, R_MAX, true, N, SLOT_LENGTH_SEC, MEAN_INTERVAL_SEC, NUM_CHANNELS)
        
        # New Correction: Calculate average AreaRatio including Boundary Effects
        area_ratio = ENABLE_CARRIER_SENSE ? calculate_average_area_ratio(R_MAX, SHADOWING_STD_DB, CS_THRESHOLD_DBM, TX_POWER_DBM, REF_PATH_LOSS_DB, PATH_LOSS_EXP) : 0.0
        
        # M is now derived directly from area_ratio (probability that two random nodes sense each other)
        M_slot = (N * area_ratio * SLOT_LENGTH_SEC / MEAN_INTERVAL_SEC) / NUM_CHANNELS
        p_lbt = M_slot > 0 ? (1.0 - exp(-M_slot)) / M_slot : 1.0
        
        G_active_sync_slot = G_sync_slot * p_lbt
        G_active_sync = G_sync * p_lbt
        G_active_async = G_async * p_lbt
        
        # 3. Capture Success Rate (Poisson Sum)
        p_succ_sync = 0.0
        
        # CS有効時、送信中のノードAに干渉してくるのはAを感知できないhidden nodeのみ。
        # CSで検知できる範囲のノードはAの送信を感知して延期するため干渉=0とする。
        # Sync packets: interference from hidden sync nodes + hidden async nodes overlapping the slot
        G_hidden_sync = G_sync_slot * (1.0 - area_ratio)
        G_hidden_async_vuln = (N * (1.0 - sync_rate) / MEAN_INTERVAL_SEC * (SLOT_LENGTH_SEC + airtime_sec)) / NUM_CHANNELS * (1.0 - area_ratio)
        G_total_interf_sync = G_hidden_sync + G_hidden_async_vuln
        
        for k in 0:20
            p_k = exp(-G_total_interf_sync) * (G_total_interf_sync^k) / factorial(big(k))
            p_succ_sync += Float64(p_k) * P_capture_with_k_interferers(k, PATH_LOSS_EXP, REQUIRED_SIR_DB, ENABLE_CAPTURE_EFFECT)
        end
        
        # For Async packets, they are Pure ALOHA (No Slot) -> 2T window vulnerability.
        # Similarly, only hidden nodes contribute to interference.
        G_total_interf_async = 2 * (G_sync * (1.0 - area_ratio) + G_async * (1.0 - area_ratio))
        p_succ_async = exp(-G_total_interf_async)
        
        # --- 4. Final Aggregation ---
        p_succ_total = sync_rate * p_succ_sync + (1.0 - sync_rate) * p_succ_async

        # --- 5. PER Definitions ---
        
        # A. One-shot PER: LBT failure (= back-off) + collision loss
        one_shot_per = 1.0 - (p_lbt * p_succ_total)
        
        # B. Attempt PER (= 全送信パケット中でGWに届かなかった割合)
        # CS有効時の干渉は hidden node のみなので p_succ_total がそれを反映済み
        attempt_per = 1.0 - p_succ_total
        
        # C. retry_per は attempt_per と同義（CS有効時 hidden node が支配的）
        retry_per = attempt_per
        
        push!(df, (Float64(N), one_shot_per, attempt_per, retry_per, 1.0 - p_lbt))
    end
    
    file_name = joinpath(output_dir, "theory_poisson_SF$(SF)_sync$(round(Int,sync_rate*100))_$(timestamp).csv")
    CSV.write(file_name, df)
    
    # Report N=500
    row_500 = df[df.num_terminals .== 500, :]
    if !isempty(row_500)
        @printf("  - N=500: \n")
        @printf("    * One-shot PER:       %.4f (LBT failure = loss)\n", row_500.one_shot_per[1])
        @printf("    * Retry-aware PER:     %.4f (Hidden Nodes only, matches Prop.jl)\n", row_500.retry_per[1])
        @printf("    * Collision-only PER:   %.4f (Attempt PER)\n", row_500.attempt_per[1])
        @printf("    * CS Block Rate:      %.4f (Probability of transmission avoidance by CS)\n", row_500.cs_block_prob[1])
    end
end

println("\nDone.")