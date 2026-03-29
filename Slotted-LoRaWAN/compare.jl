using Printf
using DataFrames
using CSV

# --- Mute some redundant prints ---
original_stdout = stdout

include("src/Prop.jl")

N_RANGE = 100:100:500
results_prop = Float64[]

println("Running Prop.jl with matching parameters...")
for N in N_RANGE
    params = create_integrated_params()
    params.num_terminals = N
    params.area_size_m = 1000.0 # radius 500m
    params.mean_event_interval_ms = 30000.0 # 30s
    params.slot_length_ms = 200.0
    params.spreading_factor = 8
    params.lora_payload_bytes = 10
    params.simulation_duration_ms = 180000.0 # 3 min (Fast execution)
    params.enable_file_output = false
    params.enable_plot_output = false
    
    # Matching thresholds in theory.jl
    params.tx_power_dbm = 13.0
    params.cs_threshold_dbm = -80.0
    params.shadowing_std_db = 8.0
    params.pass_loss_exp = 2.7
    params.num_channels = 8
    
    res = run_integrated_simulation_with_params(params)
    push!(results_prop, res["original_per"])
end

println("Running theory.jl...")
cd("analysis")
try
    run(`julia theory.jl`)
catch e
    println("Error running theory: ", e)
end
cd("..")

# Find the most recently created theory CSV
results_dir = joinpath(@__DIR__, "results", "analysis")
csv_files = filter(f -> startswith(f, "theory_SF8_sync100_"), readdir(results_dir))
sort!(csv_files)
latest_csv = joinpath(results_dir, csv_files[end])

df_theory = CSV.read(latest_csv, DataFrame)

# Now compare
println("\n" * "="^50)
println("Comparison of Original PER (Message Loss Rate)")
println("="^50)
@printf("%-10s | %-15s | %-15s\n", "N", "theory.jl", "Prop.jl")
println("-"^45)

for i in 1:length(N_RANGE)
    n = N_RANGE[i]
    th_per = df_theory[df_theory.num_terminals .== n, :original_per][1]
    pr_per = results_prop[i]
    @printf("%-10d | %-15.6f | %-15.6f\n", n, th_per, pr_per)
end
println("="^50)
