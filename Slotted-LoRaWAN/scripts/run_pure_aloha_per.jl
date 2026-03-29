using Distributed

if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)
    addprocs(num_workers)
end

@everywhere include(joinpath(@__DIR__, "..", "src", "Pure_ALOHA.jl"))

function evaluate_pure_aloha_parallel()
    TERMINAL_COUNTS = [100,200,300,400,500]
    NUM_TRIALS = 100
    SF = 8
    ENABLE_ACK = false # Pure ALOHA: No retransmissions

    out_dir = joinpath(@__DIR__, "..", "results", "parallel_evaluation", "Pure_ALOHA")
    mkpath(out_dir)
    
    summary_data = []
    
    for num_terminals in TERMINAL_COUNTS
        println("Evaluating Pure ALOHA: $num_terminals terminals")
        params = create_integrated_params()
        params.num_terminals = num_terminals
        params.spreading_factor = SF
        params.enable_ack = ENABLE_ACK
        params.enable_detailed_logs = false
        
        results_list = pmap(1:NUM_TRIALS) do i
            original_stdout = stdout
            redirect_stdout(devnull)
            res = run_integrated_simulation_with_params(params)
            redirect_stdout(original_stdout)
            return res
        end
        
        mean_per = mean([r["original_per"] for r in results_list])
        mean_thr = mean([r["throughput_bps"] for r in results_list])
        mean_norm_thr = mean([r["norm_throughput"] for r in results_list])
        mean_total_generated = mean([r["total_generated"] for r in results_list])
        mean_total_packets = mean([r["total_packets"] for r in results_list])
        
        push!(summary_data, (
            num_terminals = num_terminals,
            mean_per = mean_per,
            mean_throughput_bps = mean_thr,
            mean_norm_throughput = mean_norm_thr,
            mean_total_generated = mean_total_generated,
            mean_total_packets = mean_total_packets
        ))
    end
    
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    summary_file = joinpath(out_dir, "Pure_ALOHA_summary_SF$(SF)_$(timestamp).csv")
    open(summary_file, "w") do io
        println(io, "num_terminals,mean_per,mean_throughput_bps,mean_norm_throughput,mean_total_generated,mean_total_packets")
        for d in summary_data
            @printf(io, "%d,%.6f,%.2f,%.6f,%.2f,%.2f\n", 
                    d.num_terminals, d.mean_per, d.mean_throughput_bps, d.mean_norm_throughput,
                    d.mean_total_generated, d.mean_total_packets)
        end
    end
    println("Summary saved to: $summary_file")
end

if abspath(PROGRAM_FILE) == @__FILE__
    evaluate_pure_aloha_parallel()
end
