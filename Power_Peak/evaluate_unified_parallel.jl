# ============================================================
# Unified PER Evaluation Script (Parallel)
# Integrated LoRa と CSMA/CA の両方を評価
# ============================================================

using Distributed
using Statistics
using Random
using Dates
using Printf
using Plots

# 並列ワーカーの追加
if nprocs() == 1
    num_workers = min(Sys.CPU_THREADS, 8)
    addprocs(num_workers)
    println("Added $num_workers worker processes")
end

# 全ワーカーで両方のシミュレーションを読み込み
@everywhere include("simulation_wrapper_fast.jl")
@everywhere include("csma_ca_wrapper_fast.jl")

# ============================================================
# 評価パラメータ
# ============================================================

TERMINAL_COUNTS = [10, 20, 30, 40, 50, 75, 100]
NUM_RUNS = 30
CHANNEL_CONFIGS = [1, 8]

# ============================================================
# 結果保存用構造体
# ============================================================

struct EvaluationResult
    sim_type::String  # "Integrated" or "CSMA/CA"
    num_terminals::Int
    num_channels::Int
    run_id::Int
    total_packets::Int
    success::Int
    collisions::Int
    per::Float64
end

# ============================================================
# 並列シミュレーション実行関数
# ============================================================

@everywhere function run_simulation_worker(sim_type::String, num_terminals::Int, num_channels::Int, seed::Int)
    try
        if sim_type == "Integrated"
            results = run_simulation_with_params_fast(num_terminals, num_channels, seed)
        else  # CSMA/CA
            results = run_csma_ca_with_params_fast(num_terminals, num_channels, seed)
        end
        
        return (
            sim_type,
            num_terminals,
            num_channels,
            seed,
            results["total_packets"],
            results["success"],
            results["collisions"],
            results["per"]
        )
    catch e
        println("Error in $sim_type simulation (terminals=$num_terminals, channels=$num_channels, seed=$seed): $e")
        rethrow(e)
    end
end

# ============================================================
# 統計処理関数
# ============================================================

function calculate_statistics(results::Vector{EvaluationResult})
    pers = [r.per for r in results]
    
    mean_per = mean(pers)
    std_per = std(pers)
    min_per = minimum(pers)
    max_per = maximum(pers)
    median_per = median(pers)
    
    n = length(pers)
    ci_95 = 1.96 * std_per / sqrt(n)
    
    return (
        mean = mean_per,
        std = std_per,
        min = min_per,
        max = max_per,
        median = median_per,
        ci_95 = ci_95,
        n = n
    )
end

# ============================================================
# メイン評価ループ
# ============================================================

function run_unified_evaluation()
    println("="^60)
    println("   Unified PER Evaluation (Parallel)")
    println("   Integrated LoRa vs CSMA/CA")
    println("   Multiple Runs: $NUM_RUNS per configuration")
    println("   Workers: $(nprocs()-1)")
    println("="^60)
    println()
    
    all_results = EvaluationResult[]
    
    for sim_type in ["Integrated", "CSMA/CA"]
        println("\n" * "="^60)
        println("Evaluating $sim_type Simulation")
        println("="^60)
        
        for num_channels in CHANNEL_CONFIGS
            println("\n" * "-"^60)
            println("$sim_type with $num_channels channel(s)")
            println("-"^60)
            
            for num_terminals in TERMINAL_COUNTS
                println("\nTerminals: $num_terminals")
                print("  Running $NUM_RUNS simulations in parallel... ")
                flush(stdout)
                
                start_time = time()
                
                tasks = [(sim_type, num_terminals, num_channels, run) for run in 1:NUM_RUNS]
                raw_results = pmap(t -> run_simulation_worker(t...), tasks)
                
                elapsed = time() - start_time
                
                config_results = [EvaluationResult(r...) for r in raw_results]
                append!(all_results, config_results)
                
                println("✓ ($(round(elapsed, digits=1))s)")
                
                stats = calculate_statistics(config_results)
                
                println("  Results:")
                @printf("    Mean PER:  %.2f ± %.2f %% (95%% CI)\n", stats.mean, stats.ci_95)
                @printf("    Std Dev:   %.2f %%\n", stats.std)
                @printf("    Range:     [%.2f, %.2f] %%\n", stats.min, stats.max)
                @printf("    Median:    %.2f %%\n", stats.median)
                @printf("    Throughput: %.1f sims/sec\n", NUM_RUNS / elapsed)
            end
        end
    end
    
    # 結果保存とグラフ生成
    save_results_and_plot(all_results)
    
    println("\n" * "="^60)
    println("Evaluation Complete!")
    println("="^60)
end

# ============================================================
# 結果保存とグラフ生成
# ============================================================

function save_results_and_plot(results::Vector{EvaluationResult})
    output_dir = "result_evaluation"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    # グラフ用データ（プロット可能な形式）
    plot_file = joinpath(output_dir, "plot_data_$(timestamp).csv")
    open(plot_file, "w") do io
        println(io, "SimType,NumChannels,NumTerminals,MeanPER,CI95Lower,CI95Upper")
        
        for sim_type in ["Integrated", "CSMA/CA"]
            for num_channels in CHANNEL_CONFIGS
                for num_terminals in TERMINAL_COUNTS
                    config_results = filter(r -> r.sim_type == sim_type && 
                                                  r.num_terminals == num_terminals && 
                                                  r.num_channels == num_channels, results)
                    if !isempty(config_results)
                        stats = calculate_statistics(config_results)
                        @printf(io, "%s,%d,%d,%.2f,%.2f,%.2f\n",
                            sim_type, num_channels, num_terminals,
                            stats.mean, stats.mean - stats.ci_95, stats.mean + stats.ci_95)
                    end
                end
            end
        end
    end
    println("\nPlot data saved to: $plot_file")
    println("  Format: SimType, NumChannels, NumTerminals, MeanPER, CI95Lower, CI95Upper")
    println("  Ready for plotting: X-axis=NumTerminals, Y-axis=MeanPER")
    
    # 詳細結果
    detail_file = joinpath(output_dir, "evaluation_detail_$(timestamp).csv")
    open(detail_file, "w") do io
        println(io, "SimType,NumTerminals,NumChannels,RunID,TotalPackets,Success,Collisions,PER")
        for r in results
            println(io, "$(r.sim_type),$(r.num_terminals),$(r.num_channels),$(r.run_id),$(r.total_packets),$(r.success),$(r.collisions),$(r.per)")
        end
    end
    println("Detail results saved to: $detail_file")
    
    # グラフ生成
    try
        generate_comparison_plot(plot_file, output_dir, timestamp)
    catch e
        println("Warning: Could not generate plot: $e")
    end
end

# ============================================================
# グラフ生成
# ============================================================

function generate_comparison_plot(plot_data_file::String, output_dir::String, timestamp::String)
    # CSVからデータ読み込み
    data = []
    open(plot_data_file, "r") do io
        readline(io)  # ヘッダースキップ
        for line in eachline(io)
            parts = split(line, ",")
            push!(data, (
                sim_type = parts[1],
                num_channels = parse(Int, parts[2]),
                num_terminals = parse(Int, parts[3]),
                mean_per = parse(Float64, parts[4]),
                ci_lower = parse(Float64, parts[5]),
                ci_upper = parse(Float64, parts[6])
            ))
        end
    end
    
    # チャネル数ごとにプロット
    for num_channels in CHANNEL_CONFIGS
        p = plot(title="PER vs Number of Terminals ($num_channels Channel(s))",
                 xlabel="Number of Terminals",
                 ylabel="Packet Error Rate (%)",
                 legend=:topleft,
                 size=(800, 600),
                 grid=true)
        
        for sim_type in ["Integrated", "CSMA/CA"]
            filtered = filter(d -> d.sim_type == sim_type && d.num_channels == num_channels, data)
            sort!(filtered, by=d -> d.num_terminals)
            
            terminals = [d.num_terminals for d in filtered]
            means = [d.mean_per for d in filtered]
            ci_lower = [d.ci_lower for d in filtered]
            ci_upper = [d.ci_upper for d in filtered]
            
            plot!(p, terminals, means, 
                  label="$sim_type",
                  marker=:circle,
                  linewidth=2,
                  ribbon=(means .- ci_lower, ci_upper .- means),
                  fillalpha=0.2)
        end
        
        plot_filename = joinpath(output_dir, "per_comparison_$(num_channels)ch_$(timestamp).png")
        savefig(p, plot_filename)
        println("Plot saved to: $plot_filename")
    end
end

# ============================================================
# 実行
# ============================================================

run_unified_evaluation()
