# ===== 適応的閾値設定機能 =====

using Random, Statistics, LinearAlgebra

# ===== 閾値計算方法 =====
abstract type ThresholdMethod end

struct FixedThreshold <: ThresholdMethod
    threshold_dbm::Float64
end

struct DistanceBasedThreshold <: ThresholdMethod
    base_threshold_dbm::Float64
    distance_factor::Float64  # 距離に応じた閾値調整係数
end

struct SNRBasedThreshold <: ThresholdMethod
    base_threshold_dbm::Float64
    snr_margin_db::Float64    # SNRマージン
end

struct AdaptiveThreshold <: ThresholdMethod
    base_threshold_dbm::Float64
    noise_floor_dbm::Float64
    margin_db::Float64        # ノイズフロアからのマージン
    
    # コンストラクタ
    function AdaptiveThreshold(;
        base_threshold_dbm::Float64,
        noise_floor_dbm::Float64,
        margin_db::Float64
    )
        new(base_threshold_dbm, noise_floor_dbm, margin_db)
    end
end

# ===== 閾値計算関数 =====

"""
固定閾値を計算
"""
function calculate_threshold(terminal::TerminalInfo, method::FixedThreshold)
    return 10^(method.threshold_dbm / 10) * 1e-3  # dBm → W
end

"""
距離ベース閾値を計算
"""
function calculate_threshold(terminal::TerminalInfo, method::DistanceBasedThreshold)
    # 距離に応じて閾値を調整
    distance_factor = 1.0 + method.distance_factor * (terminal.distance_m / 100.0)
    adjusted_threshold_dbm = method.base_threshold_dbm + 10 * log10(distance_factor)
    return 10^(adjusted_threshold_dbm / 10) * 1e-3  # dBm → W
end

"""
SNRベース閾値を計算
"""
function calculate_threshold(terminal::TerminalInfo, method::SNRBasedThreshold)
    # 受信電力からSNRを推定して閾値を調整
    estimated_snr_db = terminal.rx_power_dbm - (-100.0)  # ノイズフロアを-100dBmと仮定
    adjusted_threshold_dbm = method.base_threshold_dbm - method.snr_margin_db * (estimated_snr_db / 20.0)
    return 10^(adjusted_threshold_dbm / 10) * 1e-3  # dBm → W
end

"""
適応的閾値を計算（ノイズフロアベース）
"""
function calculate_threshold(terminal::TerminalInfo, method::AdaptiveThreshold)
    # ノイズフロア + マージンで閾値を設定
    threshold_dbm = method.noise_floor_dbm + method.margin_db
    return 10^(threshold_dbm / 10) * 1e-3  # dBm → W
end

# ===== 端末別閾値設定 =====

"""
端末の位置と受信電力に基づいて個別の閾値を計算
"""
function calculate_terminal_thresholds(terminals::Vector{TerminalInfo}, method::ThresholdMethod)
    thresholds = Float64[]
    
    for terminal in terminals
        threshold = calculate_threshold(terminal, method)
        push!(thresholds, threshold)
    end
    
    return thresholds
end

# ===== 閾値設定のプリセット =====

"""
近距離端末用の閾値設定（高感度）
"""
function create_near_terminal_threshold()
    return DistanceBasedThreshold(
        base_threshold_dbm = -90.0,    # 基本閾値
        distance_factor = 0.5          # 距離係数
    )
end

"""
中距離端末用の閾値設定（中感度）
"""
function create_medium_terminal_threshold()
    return DistanceBasedThreshold(
        base_threshold_dbm = -85.0,    # 基本閾値
        distance_factor = 0.3          # 距離係数
    )
end

"""
遠距離端末用の閾値設定（低感度）
"""
function create_far_terminal_threshold()
    return DistanceBasedThreshold(
        base_threshold_dbm = -80.0,    # 基本閾値
        distance_factor = 0.1          # 距離係数
    )
end

"""
適応的閾値設定（推奨）
"""
function create_adaptive_threshold()
    return AdaptiveThreshold(
        base_threshold_dbm = -90.0,    # 基本閾値
        noise_floor_dbm = -100.0,      # ノイズフロア
        margin_db = 10.0               # マージン
    )
end

# ===== 閾値最適化 =====

"""
端末の受信電力分布に基づいて最適な閾値を計算
"""
function optimize_thresholds(terminals::Vector{TerminalInfo}, target_detection_rate::Float64 = 0.8)
    # 受信電力の分布を分析
    rx_powers = [terminal.rx_power_dbm for terminal in terminals]
    mean_power = mean(rx_powers)
    std_power = std(rx_powers)
    
    # 目標検出率に基づいて閾値を調整
    # 正規分布を仮定して、目標検出率に対応する閾値を計算
    z_score = quantile(Normal(), 1.0 - target_detection_rate)
    optimal_threshold_dbm = mean_power + z_score * std_power - 5.0  # 5dBマージン
    
    # 各端末の個別閾値を計算
    thresholds = Float64[]
    for terminal in terminals
        # 端末の受信電力に基づいて閾値を調整
        terminal_threshold_dbm = optimal_threshold_dbm + (terminal.rx_power_dbm - mean_power) * 0.5
        threshold = 10^(terminal_threshold_dbm / 10) * 1e-3
        push!(thresholds, threshold)
    end
    
    return thresholds
end

# ===== 閾値設定の表示 =====

"""
閾値設定の詳細を表示
"""
function display_threshold_settings(terminals::Vector{TerminalInfo}, thresholds::Vector{Float64})
    println("=== 端末別閾値設定 ===")
    
    for (i, terminal) in enumerate(terminals)
        threshold_dbm = 10 * log10(thresholds[i] * 1000)
        margin_db = terminal.rx_power_dbm - threshold_dbm
        
        println("端末 $i:")
        println("  - 位置: ($(round(terminal.x_m, digits=1)), $(round(terminal.y_m, digits=1))) m")
        println("  - 距離: $(round(terminal.distance_m, digits=1)) m")
        println("  - 受信電力: $(round(terminal.rx_power_dbm, digits=1)) dBm")
        println("  - 閾値: $(round(threshold_dbm, digits=1)) dBm")
        println("  - マージン: $(round(margin_db, digits=1)) dB")
        println()
    end
end
