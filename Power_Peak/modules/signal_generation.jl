# ===== 信号生成関数 =====

using Random, Statistics, FFTW, LinearAlgebra

# ===== 信号パラメータ =====
struct SignalParameters
    duration_s::Float64       # 信号持続時間（s）
    center_frequency_hz::Float64  # 中心周波数（Hz）
    bandwidth_hz::Float64    # 帯域幅（Hz）
    sampling_rate_hz::Float64    # サンプリングレート（Hz）
    tx_power_dbm::Float64    # 送信電力（dBm）
end

# ===== OFDM信号生成（QPSK変調） =====
function generate_ofdm_signal_with_qpsk(bandwidth_hz::Float64, duration_s::Float64, fs_hz::Float64, tx_power_dbm::Float64)
    
    # サンプル数とFFTサイズを計算
    N = Int(ceil(duration_s * fs_hz))
    N = max(N, 1)  # 最低1サンプル
    
    # 帯域幅に対応する周波数ビン数を計算
    bw_bins = Int(ceil(bandwidth_hz / fs_hz * N))
    bw_bins = max(bw_bins, 1)  # 最低1ビン
    
    # 周波数ドメインでQPSKシンボルを生成
    freq_domain = zeros(ComplexF64, N)
    
    # 正の周波数成分（中央から帯域幅分）
    start_bin = N ÷ 2 - bw_bins ÷ 2
    end_bin = N ÷ 2 + bw_bins ÷ 2
    
    for i in start_bin:end_bin
        if i > 0 && i <= N
            # QPSKシンボル生成（ランダム）
            real_part = (rand() > 0.5 ? 1.0 : -1.0) / sqrt(2)
            imag_part = (rand() > 0.5 ? 1.0 : -1.0) / sqrt(2)
            freq_domain[i] = real_part + 1im * imag_part
        end
    end
    
    # 負の周波数成分（共役対称）
    for i in 1:(N÷2)
        if N - i + 1 <= N && N - i + 1 > N ÷ 2
            freq_domain[N - i + 1] = conj(freq_domain[i])
        end
    end
    
    # IFFTで時間ドメイン信号を生成
    time_domain = ifft(freq_domain)
    
    # 送信電力に合わせてスケーリング（正規化の代わり）
    linear_power = 10^(tx_power_dbm / 10) * 1e-3  # dBm → W
    time_domain = time_domain * sqrt(linear_power)
    
    return time_domain
end

# ===== 周期的同期信号生成 =====
function generate_periodic_sync_signals(params::SignalParameters, interval_ms::Float64, total_duration_ms::Float64)
    signal_duration_ms = params.duration_s * 1000
    sampling_rate_hz = params.sampling_rate_hz
    total_samples = Int(ceil(total_duration_ms * sampling_rate_hz / 1000))
    
    # 信号送信タイミングを計算
    signal_times = Float64[]
    start_delay_ms = 10.0 +rand() * 20.0  # 開始遅延（ms）
    current_time = start_delay_ms
    while current_time + signal_duration_ms <= total_duration_ms
        push!(signal_times, current_time)
        current_time += interval_ms
    end
    
    # 送信信号配列を初期化
    tx_signal = zeros(ComplexF64, total_samples)
    time_axis = collect((0:total_samples-1) / sampling_rate_hz * 1000)  # ms単位
    
    # 各送信タイミングで信号を生成・配置
    signal_count = 0
    for signal_time in signal_times
        # 信号生成（各送信で異なるランダムシードを使用）
        signal_samples = generate_ofdm_signal_with_qpsk(params.bandwidth_hz, params.duration_s, params.sampling_rate_hz, params.tx_power_dbm)
        
        # 信号配置
        start_sample = Int(ceil(signal_time * sampling_rate_hz / 1000)) + 1
        end_sample = min(start_sample + length(signal_samples) - 1, total_samples)
        
        if start_sample <= total_samples
            signal_length = end_sample - start_sample + 1
            tx_signal[start_sample:end_sample] = signal_samples[1:signal_length]
            signal_count += 1
            println("送信信号 $(signal_count): 時刻 $(signal_time) ms ～ $(signal_time + signal_duration_ms) ms")
        end
    end
    
    println("送信側パラメータ:")
    println("• 開始遅延: $(start_delay_ms) ms")
    println("• 信号間隔: $(interval_ms) ms")
    println("• 信号持続時間: $(signal_duration_ms) ms")
    println("• 総持続時間: $(total_duration_ms) ms")
    println("• 送信信号数: $(signal_count)")
    println("• 送信電力: $(params.tx_power_dbm) dBm")
    println()
    
    # Return ideal beacon times (center of each transmission) for accuracy analysis
    return time_axis, tx_signal, signal_count, signal_times
end

# ===== 信号生成のテスト関数 =====
function test_signal_generation()
    println("=== 信号生成テスト ===")
    
    # パラメータ設定
    params = SignalParameters(
        66.67e-6,          # 信号持続時間（s）- QPSK_resampling.jlに合わせる
        4.7e9,             # 中心周波数（Hz）
        3.6e6,             # 帯域幅（Hz）- QPSK_resampling.jlに合わせる
        7.68e6,            # サンプリングレート（Hz）- QPSK_resampling.jlに合わせる
        43.0               # 送信電力（dBm）
    )
    
    # 信号生成
    time_axis, tx_signal, signal_count, ideal_beacon_times = generate_periodic_sync_signals(params, 20.0, 110.0)
    
    println("生成された信号:")
    println("• サンプル数: $(length(tx_signal))")
    println("• 信号数: $(signal_count)")
    println("• 平均電力: $(round(mean(abs2.(tx_signal)) * 1000, digits=3)) mW")
    println("• 最大電力: $(round(maximum(abs2.(tx_signal)) * 1000, digits=3)) mW")
    
    return time_axis, tx_signal, signal_count
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    test_signal_generation()
end