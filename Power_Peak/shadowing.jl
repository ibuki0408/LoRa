# ===== シャドウイング計算関数 =====

using Random, Statistics, Distributions

# ===== シャドウイングパラメータ =====
struct ShadowingParameters
    enabled::Bool            # シャドウイング有効/無効
    std_db::Float64         # シャドウイング標準偏差（dB）
    correlation_distance_m::Float64  # 相関距離（m）
    correlation_coefficient::Float64  # 相関係数
end

# ===== シャドウイング計算 =====
function calculate_shadowing(shadowing_params::ShadowingParameters, distance_m::Float64, fixed_value::Bool=false)
    if !shadowing_params.enabled
        return 0.0
    end
    if fixed_value
        return shadowing_params.std_db
    else
        shadowing_db = randn() * shadowing_params.std_db
        return shadowing_db
    end
end

# ===== シャドウイング計算（簡易版） =====
function calculate_shadowing_simple(std_db::Float64, enabled::Bool=true)
    if !enabled
        return 0.0
    end
    shadowing_db = randn() * std_db
    return shadowing_db
end

# ===== 相関シャドウイング計算 =====
function calculate_correlated_shadowing(shadowing_params::ShadowingParameters, distance_m::Float64, 
                                       previous_shadowing_db::Float64=0.0)
    if !shadowing_params.enabled
        return 0.0
    end
    
    # 相関係数による重み付け
    correlation_weight = shadowing_params.correlation_coefficient
    new_shadowing_db = randn() * shadowing_params.std_db
    
    # 前回の値との相関を考慮
    correlated_shadowing = correlation_weight * previous_shadowing_db + 
                           sqrt(1 - correlation_weight^2) * new_shadowing_db
    
    return correlated_shadowing
end

# ===== シャドウイング統計分析 =====
function analyze_shadowing_statistics(shadowing_params::ShadowingParameters, num_samples::Int=1000)
    shadowing_values = Float64[]
    
    for i in 1:num_samples
        shadowing_db = calculate_shadowing(shadowing_params, 50.0, false)
        push!(shadowing_values, shadowing_db)
    end
    
    mean_shadowing = mean(shadowing_values)
    std_shadowing = std(shadowing_values)
    
    println("シャドウイング統計分析:")
    println("• サンプル数: $(num_samples)")
    println("• 平均: $(round(mean_shadowing, digits=3)) dB")
    println("• 標準偏差: $(round(std_shadowing, digits=3)) dB")
    println("• 理論標準偏差: $(shadowing_params.std_db) dB")
    println("• 最小値: $(round(minimum(shadowing_values), digits=3)) dB")
    println("• 最大値: $(round(maximum(shadowing_values), digits=3)) dB")
    
    return shadowing_values
end

# ===== シャドウイング計算のテスト関数 =====
function test_shadowing()
    Random.seed!(1234)  # ランダムシードを固定
    println("=== シャドウイング計算テスト ===")
    
    # パラメータ設定
    shadowing_params = ShadowingParameters(
        true,               # シャドウイング有効
        8.0,                # 標準偏差（dB）
        50.0,               # 相関距離（m）
        0.5                 # 相関係数
    )
    
    # シャドウイング計算
    shadowing_db = calculate_shadowing(shadowing_params, 50.0, false)
    
    println("シャドウイング計算結果:")
    println("• シャドウイング: $(round(shadowing_db, digits=2)) dB")
    println("• 標準偏差: $(shadowing_params.std_db) dB")
    
    # 統計分析
    println("\n統計分析:")
    analyze_shadowing_statistics(shadowing_params, 1000)
    
    return shadowing_db
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    test_shadowing()
end
