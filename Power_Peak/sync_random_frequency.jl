using Random, Statistics, Printf, FFTW, LinearAlgebra, Plots, DSP, Distributions

# ===== 構造体の定義 =====

# ビーコン信号のパラメータ
struct BeaconParameters
    duration_s::Float64       # 信号持続時間（s）
    peak_power_dBm::Float64   # ピーク電力（dBm）
    center_frequency_hz::Float64 # 中心周波数（Hz）
    bandwidth_hz::Float64     # 帯域幅（Hz）
end

# 端末
mutable struct Terminal
    id::Int
    x::Float64 # X座標（km）
    y::Float64 # Y座標（km）
    detected_peaks::Vector{Tuple{Float64, Float64, Float64}} # (検出時刻, 電力, 信頼度)
end

# ビーコン送信機
mutable struct BeaconTransmitter
    x::Float64 # X座標（km）
    y::Float64 # Y座標（km）
    parameters::BeaconParameters
end

# シミュレーション環境
mutable struct SimulationEnvironment
    terminals::Vector{Terminal}
    beacon::BeaconTransmitter
    path_loss_alpha::Float64
    path_loss_beta::Float64
    path_loss_gamma::Float64
    shadowing_std::Float64
    noise_floor_dBm::Float64
    sampling_rate_hz::Float64
end

# ===== 信号生成関数 (教授の提案に基づく) =====

"""
周波数領域で定義したランダムな広帯域信号（OFDMライクな信号）を生成する
"""
function generate_random_ofdm_signal(bandwidth_hz::Float64, duration_s::Float64, fs_hz::Float64)
    # FFTのサイズ（＝時間領域のサンプル数）を決定
    N = Int(ceil(duration_s * fs_hz))
    
    # 周波数領域の配列をゼロで初期化
    freq_domain_signal = zeros(ComplexF64, N)
    
    # 信号を配置する帯域幅を周波数ビンの数に変換
    bw_bins = Int(floor(bandwidth_hz / (fs_hz / N)))
    
    # 中心付近の指定帯域幅のビンにランダムな複素数を設定
    # 直流成分(0Hz)を中心に、左右にbw_bins/2ずつ配置する
    start_bin = div(N, 2) - div(bw_bins, 2)
    end_bin = div(N, 2) + div(bw_bins, 2)
    
    for i in start_bin:end_bin
        # rand(ComplexF64)でランダムな振幅と位相を持つ複素数を生成
        freq_domain_signal[i] = rand(ComplexF64)
    end
    
    # IFFTを適用して時間領域の信号に変換
    # fftshiftは直流成分(0Hz)を配列の中心から配列の先頭に戻すための操作
    time_domain_signal_complex = ifft(ifftshift(freq_domain_signal))
    
    return time_domain_signal_complex
end


# ===== 伝搬モデル =====

function calculate_path_loss(distance_km::Float64, alpha::Float64, beta::Float64, gamma::Float64, frequency_hz::Float64)
    if distance_km == 0
        return 0.0
    end
    return 10 * alpha * log10(distance_km) + beta + 10 * gamma * log10(frequency_hz / 1e6)
end

function calculate_position_based_shadowing(x::Float64, y::Float64, shadowing_std::Float64, rng::AbstractRNG)
    return randn(rng) * shadowing_std
end


# ===== ピーク検出アルゴリズム (電力検出) =====

function detect_peak_advanced(signal_power::Vector{Float64}, sampling_rate::Float64, noise_floor_dBm::Float64, threshold_dB::Float64)
    if isempty(signal_power)
        return nothing, 0.0, 0.0
    end
    
    valid_indices = findall(p -> p >= noise_floor_dBm + threshold_dB, signal_power)
    if isempty(valid_indices)
        return nothing, 0.0, 0.0
    end
    
    max_power, max_index_in_valid = findmax(signal_power[valid_indices])
    max_index = valid_indices[max_index_in_valid]
    
    peak_time = (max_index - 1) / sampling_rate
    confidence = min(1.0, (max_power - noise_floor_dBm) / (threshold_dB * 2.0))
    
    return peak_time, max_power, confidence
end


# ===== シミュレーション本体 =====

function simulate_broadcast_peak_detection(env::SimulationEnvironment, transmission_time::Float64, rng::AbstractRNG)
    println("\n--- ブロードキャスト送信時刻: $(round(transmission_time, digits=6))s ---")
    
    # 送信信号を生成 (OFDMライクな複素信号)
    base_signal_complex = generate_random_ofdm_signal(
        env.beacon.parameters.bandwidth_hz,
        env.beacon.parameters.duration_s,
        env.sampling_rate_hz
    )

    # 送信信号の電力を計算 (正規化のため)
    # abs2. は複素数の絶対値の2乗を計算し、電力に比例する
    base_signal_power = abs2.(base_signal_complex)
    avg_base_power = mean(base_signal_power)
    
    # 各端末での受信処理
    for terminal in env.terminals
        # 距離、パスロス、シャドウイング計算
        distance_km = sqrt((terminal.x - env.beacon.x)^2 + (terminal.y - env.beacon.y)^2)
        path_loss_dB = calculate_path_loss(distance_km, env.path_loss_alpha, env.path_loss_beta, env.path_loss_gamma, env.beacon.parameters.center_frequency_hz)
        shadowing_dB = calculate_position_based_shadowing(terminal.x, terminal.y, env.shadowing_std, rng)
        
        # 平均受信電力(dBm)を計算
        avg_received_power_dBm = env.beacon.parameters.peak_power_dBm - path_loss_dB - shadowing_dB
        
        # 線形スケールでの受信電力
        avg_received_power_linear = 10^(avg_received_power_dBm / 10.0)
        
        # 送信信号の電力にスケーリングを適用
        power_scaling_factor = avg_received_power_linear / avg_base_power
        received_signal_power = base_signal_power .* power_scaling_factor
        
        # dBmに変換
        received_signal_dBm = 10 .* log10.(received_signal_power)
        
        # 熱雑音を追加 (線形領域で加算し、再度dBmに)
        noise_power_linear = 10^(env.noise_floor_dBm / 10.0)
        total_signal_power_linear = received_signal_power .+ noise_power_linear
        total_signal_dBm = 10 .* log10.(total_signal_power_linear)

        # ピーク検出 (電力信号を渡す)
        peak_time, peak_power, confidence = detect_peak_advanced(
            total_signal_dBm, 
            env.sampling_rate_hz,
            env.noise_floor_dBm, 
            3.0 # 検出閾値 3dB
        )
        
        if peak_time !== nothing
            detection_time = transmission_time + peak_time
            push!(terminal.detected_peaks, (detection_time, peak_power, confidence))
            @printf("  端末%d: ピーク検出! 時刻=%.6fs, 電力=%.1fdBm (距離=%.1fm)\n",
                    terminal.id, detection_time, peak_power, distance_km * 1000)
        else
            @printf("  端末%d: ピーク未検出 (距離=%.1fm, 平均受信電力=%.1fdBm)\n",
                    terminal.id, distance_km * 1000, avg_received_power_dBm)
        end
    end
end


# ===== メイン実行関数 =====

function main()
    println("="^60)
    println("同期ビーコン信号シミュレーション (OFDMライク信号版)")
    println("="^60)

    # シミュレーションパラメータ
    NUM_TERMINALS = 10
    SIMULATION_DURATION_S = 0.5
    BEACON_INTERVAL_S = 0.1
    COVERAGE_RADIUS_KM = 1.0
    
    rng = MersenneTwister(1234) # 乱数シードを固定

    # 環境設定
    env = SimulationEnvironment(
        [], # 端末リスト
        BeaconTransmitter(
            0.0, 0.0, # 送信機位置 (km)
            BeaconParameters(
                0.001,      # 信号持続時間 (1ms)
                20.0,       # 送信電力 (dBm)
                4.7e9,      # 中心周波数 (4.7GHz, ローカル5G帯)
                10e6        # 帯域幅 (10MHz)
            )
        ),
        3.5, 9.0, 2.0,  # パスロス係数
        5.0,            # シャドウイング標準偏差 (dB)
        -95.0,          # ノイズフロア (dBm)
        20e6            # サンプリングレート (20MHz)
    )

    # 端末をランダムに配置
    for i in 1:NUM_TERMINALS
        r = COVERAGE_RADIUS_KM * sqrt(rand(rng))
        theta = 2π * rand(rng)
        terminal = Terminal(i, r * cos(theta), r * sin(theta), [])
        push!(env.terminals, terminal)
    end
    println("端末を$(length(env.terminals))台配置しました。")

    # シミュレーションループ
    beacon_times = 0.0:BEACON_INTERVAL_S:SIMULATION_DURATION_S
    for t in beacon_times
        simulate_broadcast_peak_detection(env, t, rng)
    end

    # 結果の表示
    println("\n" * "="^60)
    println("シミュレーション結果")
    println("="^60)
    total_detections = 0
    for terminal in env.terminals
        detections = length(terminal.detected_peaks)
        total_detections += detections
        println("端末$(terminal.id): ピーク検出回数 = $(detections)")
    end
    println("-"^60)
    println("総検出回数: $(total_detections)")
    println("シミュレーション完了。")
    
    # 可視化 (オプション)
    # plot_terminal_layout(env)
end


# ===== 可視化関数 (オプション) =====
function plot_terminal_layout(env::SimulationEnvironment)
    p = plot(xlabel="X (km)", ylabel="Y (km)", title="端末配置", aspect_ratio=:equal, legend=:outertopright)
    
    # 端末プロット
    scatter!(p, [t.x for t in env.terminals], [t.y for t in env.terminals], label="端末")
    
    # 送信機プロット
    scatter!(p, [env.beacon.x], [env.beacon.y], marker=:star, markersize=10, label="ビーコン送信機")
    
    # カバレッジ円
    theta = 0:0.01:2π
    plot!(p, cos.(theta), sin.(theta), linestyle=:dash, label="カバレッジ (1km)")
    
    savefig(p, "terminal_layout.png")
    println("端末配置図を'terminal_layout.png'に保存しました。")
    return p
end


# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end