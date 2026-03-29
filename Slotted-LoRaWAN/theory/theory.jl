# Theoretical Capture + LBT Analysis Script
# Incorporates: Slotted ALOHA, Sync Success Rate, Capture Effect, and LBT Suppression.
# Outputs: Simplified CSV with major metrics.

using Distributions
using Random
using Statistics
using Plots
using Printf
using CSV
using DataFrames
using Dates

# Include LoRa Airtime module
include(joinpath(@__DIR__, "..", "src", "modules", "lora_airtime.jl"))
using .LoraAirtime

# ==========================================
# Parameters (Match Prop.jl)
# ==========================================
NUM_TRIALS_MONTE_CARLO = 10000
MAX_COLLIDING_PACKETS = 40

# Geometry
R_MIN = 10.0
R_MAX = 1000.0  # area_size / 2
TX_POWER_DBM = 13.0
PATH_LOSS_EXP = 2.7
REF_PATH_LOSS_DB = 20*log10(0.92e9) - 147.55 # 920MHz
SHADOWING_STD_DB = 8.0

# Capture & LBT Thresholds
CS_THRESHOLD_DBM = -80.0
REQUIRED_SIR_DB = 6.0
REQUIRED_SNR_DB = -10.0  # SF8
NOISE_DBM = -174 + 10*log10(125e3) + 6.0 # -117 dBm approx

# Traffic Parameters
SF = 8
PAYLOAD_BYTES = 10
MEAN_INTERVAL_SEC = 30.0 # 30.0s 
SLOT_LENGTH_MS = 200.0   # シミュレータと同じスロット長
NUM_CHANNELS = 8
N_RANGE = 100:100:500
# Sync Success Rate Scenarios
SYNC_RATES = [0.95]

# ==========================================
# Derived Parameters
# ==========================================
# センシング距離 R_cs の計算 (TX_POWER - PL(d) = CS_THRESHOLD)
R_CS = 10^((TX_POWER_DBM - CS_THRESHOLD_DBM - REF_PATH_LOSS_DB) / (10 * PATH_LOSS_EXP))

# ==========================================
# Physics Functions
# ==========================================

function get_rx_power(dist_m)
    pl = REF_PATH_LOSS_DB + 10 * PATH_LOSS_EXP * log10(max(dist_m, 1.0))
    sh = randn() * SHADOWING_STD_DB
    return TX_POWER_DBM - (pl + sh)
end

function get_dist(p1, p2)
    return sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2)
end

# ==========================================
# Monte Carlo Simulation for Capture + LBT
# ==========================================

function get_capture_lbt_probability(k_sync::Int, g_async::Float64, enable_lbt::Bool)
    if k_sync == 0 return 0.0 end
    
    success_count = 0
    noise_w = 10^(NOISE_DBM/10) * 1e-3
    sir_thresh_lin = 10^(REQUIRED_SIR_DB/10)
    snr_thresh_lin = 10^(REQUIRED_SNR_DB/10)
    min_power_w = noise_w * snr_thresh_lin

    for _ in 1:NUM_TRIALS_MONTE_CARLO
        # 1. Background Async Traffic
        k_async = rand(Poisson(2 * g_async))
        k_total = k_sync + k_async
        
        # Deploy all potential terminals and pre-calculate their properties
        # index 1:k_sync are sync, k_sync+1:k_total are async
        terminal_data = []
        for i in 1:k_total
            r = sqrt(R_MIN^2 + (R_MAX^2 - R_MIN^2) * rand())
            theta = 2π * rand()
            x = r * cos(theta)
            y = r * sin(theta)
            dist_gw = sqrt(x^2 + y^2)
            
            # Shadowing is fixed for the duration of this packet/trial
            sh = randn() * SHADOWING_STD_DB
            
            push!(terminal_data, (x=x, y=y, dist_gw=dist_gw, shadowing=sh))
        end
        # 2. LBT Suppression
        active_indices = Int[]
        if enable_lbt
            order = shuffle(1:k_total)
            for i in order
                is_busy = false
                ti = terminal_data[i]
                for active_idx in active_indices
                    ta = terminal_data[active_idx]
                    d = sqrt((ti.x - ta.x)^2 + (ti.y - ta.y)^2)
                    if d <= R_CS
                        is_busy = true
                        break
                    end
                end
                if !is_busy
                    push!(active_indices, i)
                end
            end
        else
            active_indices = collect(1:k_total)
        end
        
        # Check if any SYNC packet survived
        sync_active_indices = intersect(active_indices, 1:k_sync)
        if isempty(sync_active_indices) continue end
        
        # 3. Capture Effect at Gateway
        # Note: We need the GW received power (using GW distance and trial-fixed shadowing)
        function get_power_at_gw(t_data)
            pl = REF_PATH_LOSS_DB + 10 * PATH_LOSS_EXP * log10(max(t_data.dist_gw, 1.0))
            p_dbm = TX_POWER_DBM - (pl + t_data.shadowing)
            return 10^(p_dbm/10) * 1e-3 # Watts
        end

        active_powers_w = [get_power_at_gw(terminal_data[idx]) for idx in active_indices]
        sync_powers_w = [get_power_at_gw(terminal_data[idx]) for idx in sync_active_indices]
        
        strongest_sync = maximum(sync_powers_w)
        total_active_sum = sum(active_powers_w)
        # SINR check for the strongest SYNC packet
        others_sum = total_active_sum - strongest_sync
        sir_ok = (others_sum > 0) ? (strongest_sync / others_sum >= sir_thresh_lin) : true
        snr_ok = (strongest_sync >= min_power_w)
        
        if sir_ok && snr_ok
            success_count += 1
        end
    end
    
    # Return (expected success, expected transmissions)
    # transmissions = total sync packets that actually went to air
    return success_count / NUM_TRIALS_MONTE_CARLO, 0 # Transmission count is harder to return this way
end

# Modified Monte Carlo to return multiple values
function run_monte_carlo_detailed(k_sync::Int, g_async::Float64, enable_lbt::Bool)
    if k_sync == 0 return (0.0, 0.0) end
    
    success_count = 0
    tx_count = 0
    noise_w = 10^(NOISE_DBM/10) * 1e-3
    sir_thresh_lin = 10^(REQUIRED_SIR_DB/10)
    snr_thresh_lin = 10^(REQUIRED_SNR_DB/10)
    min_power_w = noise_w * snr_thresh_lin

    for _ in 1:NUM_TRIALS_MONTE_CARLO
        # 1. Background Async
        k_async = rand(Poisson(2 * g_async))
        k_total = k_sync + k_async
        terminal_data = []
        for i in 1:k_total
            r = sqrt(R_MIN^2 + (R_MAX^2 - R_MIN^2) * rand())
            theta = 2π * rand()
            push!(terminal_data, (x=r*cos(theta), y=r*sin(theta), dist_gw=r, shadowing=randn()*SHADOWING_STD_DB))
        end

        # 2. LBT Suppression
        active_indices = Int[]
        if enable_lbt
            # シミュレータのように同時刻イベントも内部で「順次」処理される（タイブレークが発生する）ためshuffleを使用
            order = shuffle(1:k_total)
            for i in order
                is_busy = false
                ti = terminal_data[i]
                for active_idx in active_indices
                    ta = terminal_data[active_idx]
                    d = sqrt((ti.x - ta.x)^2 + (ti.y - ta.y)^2)
                    if d <= R_CS
                        is_busy = true
                        break
                    end
                end
                if !is_busy
                    push!(active_indices, i)
                end
            end
        else
            active_indices = collect(1:k_total)
        end
        
        # Count sync transmissions
        sync_active_indices = intersect(active_indices, 1:k_sync)
        if isempty(sync_active_indices) continue end
        tx_count += length(sync_active_indices)
        
        # 3. Capture Effect
        function get_power_at_gw(t_data)
            pl = REF_PATH_LOSS_DB + 10 * PATH_LOSS_EXP * log10(max(t_data.dist_gw, 1.0))
            return 10^((TX_POWER_DBM - (pl + t_data.shadowing))/10) * 1e-3
        end
        
        active_powers_w = [get_power_at_gw(terminal_data[idx]) for idx in active_indices]
        sync_powers_w = [get_power_at_gw(terminal_data[idx]) for idx in sync_active_indices]
        
        strongest_sync = maximum(sync_powers_w)
        others_sum = sum(active_powers_w) - strongest_sync
        
        if (others_sum > 0 ? (strongest_sync/others_sum >= sir_thresh_lin) : true) && (strongest_sync >= min_power_w)
            success_count += 1
        end
    end
    return success_count / NUM_TRIALS_MONTE_CARLO, tx_count / NUM_TRIALS_MONTE_CARLO
end

# ==========================================
# Main Execution
# ==========================================

println("="^60)
println("  Advanced Theoretical Analysis: Sync + Capture + LBT")
println("="^60)

println("Calculating Throughput Curves and Saving Simplified CSVs...")
lora_params = create_lora_params(SF, PAYLOAD_BYTES)
airtime_ms = calculate_lora_airtime(lora_params)
airtime_sec = airtime_ms / 1000.0

output_dir = joinpath(@__DIR__, "..", "results", "analysis")
mkpath(output_dir)
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

    cols = ["num_terminals", "original_per", "attempt_per", "mean_throughput_bps"]
    for sync_rate in SYNC_RATES
        rate_percent = round(Int, sync_rate * 100)
        df = DataFrame()
        for col in cols
            df[!, col] = Float64[]
        end
        
        for N in N_RANGE
            # トラフィック密度（アーラン）。
            # シミュレーションの「Attempt PER」に合わせるため、衝突判定のベースには Airtime を使用。
            # (シミュレータではスロット先頭でパケットが密集しても、LBTとバックオフによって
            #  時間軸上にパケットが分散（de-clumping）されるため、実効的な衝突確率は
            #  スロット長（200ms）ではなく Airtime（50ms）ベースの密度に近くなる)
            G_sync_eff = (N / MEAN_INTERVAL_SEC * airtime_sec) / NUM_CHANNELS * sync_rate
            G_async_eff = (N / MEAN_INTERVAL_SEC * airtime_sec) / NUM_CHANNELS * (1.0 - sync_rate)
            
            s_sync_total = 0.0
            tx_sync_total = 0.0
            # 同時送信パケット数 k の計算に Airtimeベースの密度を使用
            for k in 1:MAX_COLLIDING_PACKETS
                prob_k = ( (G_sync_eff)^k * exp(-G_sync_eff) ) / factorial(big(k))
                s_sync_k, tx_sync_k = run_monte_carlo_detailed(k, G_async_eff, true)
                s_sync_total += Float64(prob_k) * s_sync_k
                tx_sync_total += Float64(prob_k) * tx_sync_k
            end
            
            # パケット送信1回あたりの成功率（Attempt Success Rate）
            p_attempt_succ_sync = tx_sync_total > 0 ? (s_sync_total / tx_sync_total) : 1.0
            p_succ_async = exp(-2 * (G_sync_eff + G_async_eff)) # 非同期ALOHAの近似
            
            # Attempt PER (衝突率)
            attempt_per = 1.0 - (p_attempt_succ_sync * sync_rate + p_succ_async * (1.0 - sync_rate))
            
            # メッセージ不達率（Original PER）
            # LBTによる一時的な送信抑制（シミュレータのBufferDropに相当）を極小と仮定
            original_per = attempt_per 
            
            thr_bps = (N / MEAN_INTERVAL_SEC) * (1.0 - original_per) * PAYLOAD_BYTES * 8
            
            push!(df, (
                num_terminals = Float64(N),
                original_per = original_per,
                attempt_per = attempt_per,
                mean_throughput_bps = thr_bps
            ))
        end
        
        file_name = joinpath(output_dir, "theory_SF$(SF)_sync$(rate_percent)_$(timestamp).csv")
        CSV.write(file_name, df)
        println("Saved detailed theory CSV: $file_name")
    end

