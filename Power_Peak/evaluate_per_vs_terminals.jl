# ============================================================
# PER vs Terminal Count Evaluation Script
# 複数回実行による統計的評価
# ============================================================

using Statistics
using Random
using Dates
using Printf

# シミュレーション本体を読み込み
include("simulation_wrapper.jl")

# ============================================================
# 評価パラメータ
# ============================================================

# 端末数の範囲
TERMINAL_COUNTS = [10, 20, 30, 40, 50, 75, 100]

# 各設定での実行回数
NUM_RUNS = 30

# チャネル数の比較
CHANNEL_CONFIGS = [1, 8]  # 1チャネル vs 8チャネル

# ============================================================
# 結果保存用構造体
# ============================================================

struct EvaluationResult
    num_terminals::Int
    num_channels::Int
    run_id::Int
    total_packets::Int
    success::Int
    collisions::Int
    per::Float64
end

# ============================================================
# シミュレーション実行関数
# ============================================================

function run_single_simulation(num_terminals::Int, num_channels::Int, seed::Int)
    try
        # ラッパー関数を呼び出し
        results = run_simulation_with_params(num_terminals, num_channels, seed)
        
        return EvaluationResult(
            num_terminals,
            num_channels,
            seed,
            results["total_packets"],
            results["success"],
            results["collisions"],
            results["per"]
        )
    catch e
        println("Error in simulation (terminals=$num_terminals, channels=$num_channels, seed=$seed): $e")
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
    
    # 95% 信頼区間
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

function run_evaluation()
    println("="^60)
    println("   PER vs Terminal Count Evaluation")
    println("   Multiple Runs: $NUM_RUNS per configuration")
    println("="^60)
    println()
    
    all_results = EvaluationResult[]
    
    for num_channels in CHANNEL_CONFIGS
        println("\n" * "="^60)
        println("Evaluating with $num_channels channel(s)")
        println("="^60)
        
        for num_terminals in TERMINAL_COUNTS
            println("\nTerminals: $num_terminals")
            print("  Running: ")
            
            config_results = EvaluationResult[]
            
            for run in 1:NUM_RUNS
                # 進捗表示
                if run % 5 == 0
                    print("$run ")
                end
                
                # シミュレーション実行
                result = run_single_simulation(num_terminals, num_channels, run)
                push!(config_results, result)
                push!(all_results, result)
            end
            
            println("✓")
            
            # 統計処理
            stats = calculate_statistics(config_results)
            
            println("  Results:")
            @printf("    Mean PER:  %.2f ± %.2f %% (95%% CI)\n", stats.mean, stats.ci_95)
            @printf("    Std Dev:   %.2f %%\n", stats.std)
            @printf("    Range:     [%.2f, %.2f] %%\n", stats.min, stats.max)
            @printf("    Median:    %.2f %%\n", stats.median)
        end
    end
    
    # 結果保存
    save_results(all_results)
    
    println("\n" * "="^60)
    println("Evaluation Complete!")
    println("="^60)
end

# ============================================================
# 結果保存
# ============================================================

function save_results(results::Vector{EvaluationResult})
    output_dir = "result_evaluation"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    # 詳細結果（全実行）
    detail_file = joinpath(output_dir, "evaluation_detail_$(timestamp).csv")
    open(detail_file, "w") do io
        println(io, "NumTerminals,NumChannels,RunID,TotalPackets,Success,Collisions,PER")
        for r in results
            println(io, "$(r.num_terminals),$(r.num_channels),$(r.run_id),$(r.total_packets),$(r.success),$(r.collisions),$(r.per)")
        end
    end
    println("\nDetail results saved to: $detail_file")
    
    # 統計サマリー
    summary_file = joinpath(output_dir, "evaluation_summary_$(timestamp).csv")
    open(summary_file, "w") do io
        println(io, "NumTerminals,NumChannels,MeanPER,StdPER,CI95,MinPER,MaxPER,MedianPER,NumRuns")
        
        for num_channels in CHANNEL_CONFIGS
            for num_terminals in TERMINAL_COUNTS
                config_results = filter(r -> r.num_terminals == num_terminals && r.num_channels == num_channels, results)
                stats = calculate_statistics(config_results)
                
                @printf(io, "%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n",
                    num_terminals, num_channels,
                    stats.mean, stats.std, stats.ci_95,
                    stats.min, stats.max, stats.median, stats.n)
            end
        end
    end
    println("Summary results saved to: $summary_file")
end

# ============================================================
# 実行
# ============================================================

run_evaluation()
