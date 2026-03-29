# ===== 信号生成関数 =====
module SignalGeneration

using Random, Statistics, FFTW, LinearAlgebra

export SignalParameters, generate_ofdm_signal_with_qpsk, generate_periodic_sync_signals, calculate_beacon_times, generate_sync_signal_segment, test_signal_generation

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
    
    # --- OLD: 帯域の端（ナイキスト周波数付近）に配置していたロジック ---
    # start_bin = N ÷ 2 - bw_bins ÷ 2
    # end_bin = N ÷ 2 + bw_bins ÷ 2
    # for i in start_bin:end_bin
    #     if i > 0 && i <= N
    #         # QPSKシンボル生成（ランダム）
    #         real_part = (rand() > 0.5 ? 1.0 : -1.0) / sqrt(2)
    #         imag_part = (rand() > 0.5 ? 1.0 : -1.0) / sqrt(2)
    #         freq_domain[i] = real_part + 1im * imag_part
    #     end
    # end
    # # 負の周波数成分（共役対称）
    # for i in 1:(N÷2)
    #     if N - i + 1 <= N && N - i + 1 > N ÷ 2
    #         freq_domain[N - i + 1] = conj(freq_domain[i])
    #     end
    # end

    # --- NEW: 0 Hz (DC) を中心とした正しいベースバンド配置 ---
    # 0 Hz (index 1) から正負に広げる
    half_bw = bw_bins ÷ 2
    
    # DCおよび正の周波数 (1から順に)
    for i in 1:(half_bw + 1)
        if i <= N
            freq_domain[i] = ((rand() > 0.5 ? 1.0 : -1.0) + 1im * (rand() > 0.5 ? 1.0 : -1.0)) / sqrt(2)
        end
    end
    
    # 負の周波数 (末尾から戻る)
    for i in 1:half_bw
        idx = N - i + 1
        if idx > half_bw + 1 # 重複回避
            freq_domain[idx] = ((rand() > 0.5 ? 1.0 : -1.0) + 1im * (rand() > 0.5 ? 1.0 : -1.0)) / sqrt(2)
        end
    end
    
    # IFFTで時間ドメイン信号を生成
    time_domain = ifft(freq_domain)
    
    # 正規化（IFFTによるパワー変動を補正し、平均パワーを1にする）
    current_power = mean(abs2.(time_domain))
    if current_power > 0
        time_domain = time_domain ./ sqrt(current_power)
    end
    
    # 送信電力に合わせてスケーリング
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
            # println("送信信号 $(signal_count): 時刻 $(signal_time) ms ～ $(signal_time + signal_duration_ms) ms")
        end
    end
    
    # サマリーのみ表示
    println("同期信号生成完了:")
    println("  送信信号数: $(signal_count)")
    println("  信号間隔: $(interval_ms) ms")
    println("  総持続時間: $(total_duration_ms) ms")
    println()
    
    # Return ideal beacon times (center of each transmission) for accuracy analysis
    return time_axis, tx_signal, signal_count, signal_times
end

# ===== ビーコン送信時刻の事前計算 =====
function calculate_beacon_times(interval_ms::Float64, total_duration_ms::Float64, start_delay_ms::Float64)
    beacon_times = Float64[]
    current_time = start_delay_ms
    # 信号長は仮で計算（パラメータから本当は取るべきだが、ここではタイミングリスト作成が主目的）
    # 実際は呼び出し元で管理するか、パラメータとして渡す
    
    while current_time < total_duration_ms
        push!(beacon_times, current_time)
        current_time += interval_ms
    end
    return beacon_times
end

# ===== 指定区間の同期信号生成 (メモリ節約版) =====
function generate_sync_signal_segment(params::SignalParameters, beacon_times::Vector{Float64}, segment_start_ms::Float64, segment_end_ms::Float64)
    sampling_rate_hz = params.sampling_rate_hz
    signal_duration_ms = params.duration_s * 1000.0
    
    # セグメントの長さ（サンプル数）
    segment_duration_ms = segment_end_ms - segment_start_ms
    num_samples = Int(ceil(segment_duration_ms * sampling_rate_hz / 1000.0))
    
    segment_signal = zeros(ComplexF64, num_samples)
    
    # このセグメントに関係するビーコンを探す
    # ビーコンの開始時刻 〜 終了時刻 が セグメントの開始 〜 終了 と重なればOK
    
    # 1つのビーコン信号波形（共通で使い回して生成）
    # パフォーマンスのため、ここで1つ作ってしまう
    # ※ 本来はビーコンごとにランダムデータ変調すべきなら都度生成だが、
    #    ここでは簡略化のため、または都度生成でも可。
    #    元の実装では "各送信で異なるランダムシードを使用" とあるので都度生成する。
    
    for beacon_start_ms in beacon_times
        beacon_end_ms = beacon_start_ms + signal_duration_ms
        
        # 重なり判定
        if beacon_end_ms > segment_start_ms && beacon_start_ms < segment_end_ms
            # 重なっている
            
            # 今回のビーコン用の信号生成
            beacon_waveform = generate_ofdm_signal_with_qpsk(params.bandwidth_hz, params.duration_s, params.sampling_rate_hz, params.tx_power_dbm)
            beacon_len = length(beacon_waveform)
            
            # 配置位置の計算
            # セグメント内でのビーコン開始位置（サンプルインデックス、1-based）
            # start_idx_in_segment = (beacon_start_ms - segment_start_ms) * fs + 1
            
            start_sample_global = Int(floor(beacon_start_ms * sampling_rate_hz / 1000.0)) + 1
            start_segment_global = Int(floor(segment_start_ms * sampling_rate_hz / 1000.0)) + 1
            
            offset = start_sample_global - start_segment_global + 1
            
            # コピー範囲
            src_start = 1
            src_end = beacon_len
            dst_start = offset
            dst_end = offset + beacon_len - 1
            
            # クリップ（セグメントからはみ出す部分をカット）
            if dst_start < 1
                cut = 1 - dst_start
                src_start += cut
                dst_start = 1
            end
            if dst_end > num_samples
                cut = dst_end - num_samples
                src_end -= cut
                dst_end = num_samples
            end
            
            if src_start <= src_end && dst_start <= dst_end
                # 重畳（加算）
                segment_signal[dst_start:dst_end] .+= beacon_waveform[src_start:src_end]
            end
        end
    end
    
    return segment_signal, num_samples
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

end # module