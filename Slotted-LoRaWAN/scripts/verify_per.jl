using Distributed

# 並列ワーカーの追加
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

# 全ワーカーで環境構築
@everywhere using Statistics, Printf, Dates
@everywhere include(joinpath(@__DIR__, "..", "src", "Prop.jl"))

function run_verification()
    TERMINAL_COUNTS = [100, 200, 300, 400, 500]
    NUM_TRIALS = 10
    
    output_dir = joinpath(@__DIR__, "..", "results", "verification")
    mkpath(output_dir)
    
    summary_data = []

    for num_terminals in TERMINAL_COUNTS
        println("Evaluating: $num_terminals terminals...")
        
        params = create_integrated_params()
        params.num_terminals = num_terminals
        params.enable_carrier_sense = true
        params.enable_capture_effect = true
        params.area_size_m = 1000.0 # これで半径 500m になる
        params.use_waveform_phy = false 
        params.slot_length_ms = 200.0
        params.beacon_interval_ms = 200.0
        params.simulation_duration_ms = 600000.0 # 10 minutes for faster run
        
        # Ensure ToA is calculated
        lora_params = create_lora_params(params.spreading_factor, params.lora_payload_bytes)
        params.packet_airtime_ms = calculate_lora_airtime(lora_params)
        
        results_list = pmap(1:NUM_TRIALS) do i
            # Suppress output
            original_stdout = stdout
            redirect_stdout(devnull)
            res = nothing
            try
                res = run_integrated_simulation_with_params(params)
            finally
                redirect_stdout(original_stdout)
            end
            return res
        end
        
        collision_per_list = [r["collision_per"] for r in results_list]
        loss_rate_list = [r["message_loss_rate"] for r in results_list]
        mean_collision = mean(collision_per_list)
        mean_loss = mean(loss_rate_list)
        
        push!(summary_data, (num_terminals=num_terminals, mean_collision=mean_collision, mean_loss=mean_loss))
        @printf("  Mean Collision: %.6f | Mean Loss: %.6f\n", mean_collision, mean_loss)
    end
    
    # Save summary
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    summary_file = joinpath(output_dir, "per_full_verification_$(timestamp).csv")
    open(summary_file, "w") do io
        println(io, "num_terminals,collision_per,message_loss_rate")
        for data in summary_data
            println(io, "$(data.num_terminals),$(data.mean_collision),$(data.mean_loss)")
        end
    end
    println("Summary saved to: $summary_file")
    
    println("\nVerification Results (Simulation):")
    println("Terminals | Collision PER | Msg Loss Rate")
    println("----------|---------------|---------------")
    for data in summary_data
        @printf("%9d | %.6f      | %.6f\n", data.num_terminals, data.mean_collision, data.mean_loss)
    end
end

run_verification()
