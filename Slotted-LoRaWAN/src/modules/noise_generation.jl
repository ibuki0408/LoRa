module NoiseGeneration

using Random, Statistics

export NoiseParameters, generate_awgn_noise, generate_bandwidth_noise, generate_snr_noise, analyze_noise_statistics

# ===== ノイズパラメータ =====
struct NoiseParameters
    noise_power_dbm::Float64  # ノイズ電力（dBm）
    snr_db::Float64          # 目標SNR（dB）
end

# ===== ノイズ生成（AWGN） =====
function generate_awgn_noise(noise_power_dbm::Float64, num_samples::Int)
    noise_power_linear = 10^(noise_power_dbm / 10) * 1e-3  # dBm → W
    noise = sqrt(noise_power_linear) * (randn(num_samples) + 1im * randn(num_samples)) / sqrt(2)
    return noise
end

# ===== ノイズ生成（帯域幅指定） =====
function generate_bandwidth_noise(noise_figure_db::Float64, bandwidth_hz::Float64, num_samples::Int)
    # 熱雑音基底
    thermal_noise_dbm_per_hz = -174.0  # dBm/Hz
    
    # 帯域幅を考慮したノイズ電力
    noise_power_dbm = thermal_noise_dbm_per_hz + 10 * log10(bandwidth_hz) + noise_figure_db
    
    # ノイズ生成
    noise = generate_awgn_noise(noise_power_dbm, num_samples)
    
    return noise, noise_power_dbm
end

# ===== ノイズ生成（SNR指定） =====
function generate_snr_noise(signal_power_dbm::Float64, snr_db::Float64, num_samples::Int)
    noise_power_dbm = signal_power_dbm - snr_db
    noise = generate_awgn_noise(noise_power_dbm, num_samples)
    return noise, noise_power_dbm
end

# ===== ノイズ統計分析（修正版） =====
function analyze_noise_statistics(noise::Vector{ComplexF64})
    # リニア電力（W）を計算
    noise_power_linear = abs2.(noise)
    
    # ★修正点：先にリニア値で平均を計算する★
    mean_power_linear = mean(noise_power_linear)
    
    # 平均電力をdBmに変換
    mean_power_dbm = 10 * log10(mean_power_linear * 1000)
    
    # (参考) 各サンプルの瞬時電力（dBm）の統計
    # こちらは変動の様子を見るためのもの
    noise_power_dbm_instantaneous = 10 * log10.(noise_power_linear * 1000)
    std_power_dbm = std(noise_power_dbm_instantaneous)
    
    println("ノイズ統計分析:")
    println("• サンプル数: $(length(noise))")
    println("• 平均電力 (リニア平均): $(round(mean_power_dbm, digits=2)) dBm  <- これが正しい平均電力")
    println("• 瞬時電力の標準偏差: $(round(std_power_dbm, digits=2)) dB")
    println("• 瞬時電力の最小: $(round(minimum(noise_power_dbm_instantaneous), digits=2)) dBm")
    println("• 瞬時電力の最大: $(round(maximum(noise_power_dbm_instantaneous), digits=2)) dBm")
    
    return mean_power_dbm, std_power_dbm
end

# ===== ノイズ生成のテスト関数 =====
function test_noise_generation()
    Random.seed!(1234)  # ランダムシードを固定
    println("=== ノイズ生成テスト ===")
    
    # パラメータ設定
    noise_figure_db = 5.0
    bandwidth_hz = 125e3  # 125 kHz
    num_samples = 1000
    
    # ノイズ生成
    noise, noise_power_dbm = generate_bandwidth_noise(noise_figure_db, bandwidth_hz, num_samples)
    
    println("ノイズ生成結果:")
    println("• ノイズフィギュア: $(noise_figure_db) dB")
    println("• 帯域幅: $(bandwidth_hz/1000) kHz")
    println("• ノイズ電力: $(round(noise_power_dbm, digits=2)) dBm")
    println("• サンプル数: $(num_samples)")
    
    # 統計分析
    println("\n統計分析:")
    analyze_noise_statistics(noise)
    
    return noise, noise_power_dbm
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_noise_generation()
end

end # module
