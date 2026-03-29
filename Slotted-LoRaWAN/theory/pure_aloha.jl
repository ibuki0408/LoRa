# Classic ALOHA Theoretical Analysis
# This script calculates the standard Pure ALOHA and Slotted ALOHA throughput/PER.
# Formula: 
#   Slotted ALOHA: S = G * exp(-G)
#   Pure ALOHA:    S = G * exp(-2G)

using CSV
using DataFrames
using Printf
using Dates

# Include LoRa Airtime module
include(joinpath(@__DIR__, "..", "src", "modules", "lora_airtime.jl"))
using .LoraAirtime

# ==========================================
# Parameters (Match Prop.jl Simulation)
# ==========================================
SF = 8
PAYLOAD_BYTES = 10
MEAN_INTERVAL_SEC = 30.0
NUM_CHANNELS = 8
N_RANGE = 100:100:500

# ==========================================
# Main Execution
# ==========================================

println("="^60)
println("  Classic ALOHA Theory Analysis")
println("="^60)

# Exact Airtime calculation
lora_params = create_lora_params(SF, PAYLOAD_BYTES)
airtime_ms = calculate_lora_airtime(lora_params)
airtime_sec = airtime_ms / 1000.0
@printf("Using accurate Airtime: %.2f ms (SF%d, %d Bytes)\n", airtime_ms, SF, PAYLOAD_BYTES)


results = DataFrame(
    num_terminals = Int[],
    total_load_G = Float64[],
    pure_aloha_per = Float64[],
    slotted_aloha_per = Float64[],
    pure_aloha_thr_bps = Float64[],
    slotted_aloha_thr_bps = Float64[]
)

for N in N_RANGE
    # Traffic Load G (Erlangs per channel)
    G = (N / MEAN_INTERVAL_SEC * airtime_sec) / NUM_CHANNELS
    
    # PER = 1 - (S/G)
    # Pure ALOHA: S/G = exp(-2G)
    # Slotted ALOHA: S/G = exp(-G)
    per_pure = 1.0 - exp(-2 * G)
    per_slotted = 1.0 - exp(-G)
    
    # Throughput (bps) = (Total Load in bps) * (Success Rate)
    # Total Load in bps = (N / Interval) * Payload * 8
    total_bps = (N / MEAN_INTERVAL_SEC) * PAYLOAD_BYTES * 8
    thr_pure = total_bps * exp(-2 * G)
    thr_slotted = total_bps * exp(-G)
    
    push!(results, (N, G, per_pure, per_slotted, thr_pure, thr_slotted))
end

# Save Results
output_dir = joinpath(@__DIR__, "..", "results", "theory")
mkpath(output_dir)
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
file_path = joinpath(output_dir, "classic_aloha_SF$(SF)_$(timestamp).csv")

CSV.write(file_path, results)
println("Theory results saved to: $file_path")

# Print Summary for N=500
row_500 = results[results.num_terminals .== 500, :]
if !isempty(row_500)
    println("\nResults for N=500 (G=$(round(row_500.total_load_G[1], digits=4))):")
    @printf("  Pure ALOHA:    PER=%.4f, Thr=%.2f bps\n", row_500.pure_aloha_per[1], row_500.pure_aloha_thr_bps[1])
    @printf("  Slotted ALOHA: PER=%.4f, Thr=%.2f bps\n", row_500.slotted_aloha_per[1], row_500.slotted_aloha_thr_bps[1])
end

println("\nDone.")
