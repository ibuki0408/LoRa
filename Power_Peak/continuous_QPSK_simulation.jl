using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots

# ===== 信号パラメータ =====
struct SignalParameters
    duration_s::Float64       # 信号持続時間（s）
    center_frequency_hz::Float64 # 中心周波数（Hz）
    bandwidth_hz::Float64     # 帯域幅（Hz）
    sampling_rate_hz::Float64 # サンプリングレート（Hz）
end

# ===== QPSKシンボル生成関数 =====

"""
QPSKシンボル（4つの候補の中からランダムに1つ）を生成する。
平均電力が1になるように正規化済み。
"""
function generate_qpsk_symbol()
    # 実部と虚部がそれぞれ-1か+1のどちらかをランダムに選ぶ
    real_part = rand([-1, 1])
    imag_part = rand([-1, 1])
    
    # 複素数を生成し、平均電力が1になるようにsqrt(2)で割る
    return (real_part + im * imag_part) / sqrt(2)
end

# ===== 信号生成関数 (QPSK専用版) =====

"""
周波数領域で、各サブキャリアにランダムなQPSKシンボルを配置することで
広帯域信号を生成する。
"""
function generate_ofdm_signal_with_qpsk(bandwidth_hz::Float64, duration_s::Float64, fs_hz::Float64)
    # FFTのサイズ（＝時間領域のサンプル数）を決定
    N = Int(ceil(duration_s * fs_hz))
    
    # 周波数領域の配列をゼロで初期化
    freq_domain_signal = zeros(ComplexF64, N)
    
    # 信号を配置する帯域幅を周波数ビンの数に変換
    bw_bins = Int(floor(bandwidth_hz / (fs_hz / N)))
    
    # 中心付近の指定帯域幅のビンにランダムなシンボルを設定
    start_bin = div(N, 2) - div(bw_bins, 2)
    end_bin = div(N, 2) + div(bw_bins, 2)
    
    for i in start_bin:end_bin
        # 各サブキャリアに、ランダムなQPSKシンボルを1つ配置する
        freq_domain_signal[i] = generate_qpsk_symbol()
    end
    
    # IFFTを適用して時間領域の信号に変換
    time_domain_signal_complex = ifft(ifftshift(freq_domain_signal))
    
    return time_domain_signal_complex
end

# ===== 継続的信号生成関数 =====

"""
20ms間隔でQPSK信号を送信する継続的なシミュレーション
"""
function generate_continuous_qpsk_signals(params::SignalParameters, interval_ms::Float64, total_duration_ms::Float64)
    # パラメータ設定
    signal_duration_ms = params.duration_s * 1000  # 0.1428 ms
    sampling_rate_hz = params.sampling_rate_hz
    total_samples = Int(ceil(total_duration_ms * sampling_rate_hz / 1000))
    
    # 信号送信タイミングを計算
    signal_times = Float64[]
    current_time = 0.0
    while current_time + signal_duration_ms <= total_duration_ms
        push!(signal_times, current_time)
        current_time += interval_ms
    end
    
    println("継続的信号生成パラメータ:")
    println("• 信号間隔: $(interval_ms) ms")
    println("• 信号持続時間: $(signal_duration_ms) ms")
    println("• 総持続時間: $(total_duration_ms) ms")
    println("• 送信信号数: $(length(signal_times))")
    println("• 総サンプル数: $(total_samples)")
    println()
    
    # 継続的信号配列を初期化
    continuous_signal = zeros(ComplexF64, total_samples)
    time_axis = collect((0:total_samples-1) / sampling_rate_hz * 1000)  # ms単位
    
    # 各送信タイミングで信号を生成・配置
    signal_count = 0
    for (i, signal_time) in enumerate(signal_times)
        # 個別のQPSK信号を生成
        individual_signal = generate_ofdm_signal_with_qpsk(
            params.bandwidth_hz,
            params.duration_s,
            sampling_rate_hz
        )
        
        # 信号を継続的配列に配置
        signal_samples = length(individual_signal)
        start_sample = Int(ceil(signal_time * sampling_rate_hz / 1000)) + 1
        end_sample = min(start_sample + signal_samples - 1, total_samples)
        
        if start_sample <= total_samples
            actual_samples = end_sample - start_sample + 1
            continuous_signal[start_sample:end_sample] = individual_signal[1:actual_samples]
            signal_count += 1
            println("信号 $(signal_count): 時刻 $(round(signal_time, digits=3)) ms ～ $(round(signal_time + signal_duration_ms, digits=3)) ms")
        end
    end
    
    return time_axis, continuous_signal, signal_count
end

# ===== 継続的信号用CSV出力関数 =====

function save_continuous_signal_to_csv(time_axis::Vector{Float64}, continuous_signal::Vector{ComplexF64}, 
                                      sampling_rate::Float64, params::SignalParameters, 
                                      signal_count::Int, output_dir::String)
    N = length(continuous_signal)
    freq_axis = fftfreq(N, sampling_rate) / 1e6  # MHz単位
    
    # 時間領域データ
    time_data = DataFrame(
        time_ms = time_axis,
        real_part = real.(continuous_signal),
        imaginary_part = imag.(continuous_signal),
        amplitude = abs.(continuous_signal),
        power = abs2.(continuous_signal)
    )
    
    # 周波数領域データ（継続的信号全体のFFT）
    freq_domain = fft(continuous_signal)
    freq_data = DataFrame(
        frequency_mhz = fftshift(freq_axis),
        amplitude_spectrum = fftshift(abs.(freq_domain)),
        power_spectrum = fftshift(abs2.(freq_domain))
    )
    
    # パラメータ情報
    params_data = DataFrame(
        parameter = ["Duration (ms)", "Bandwidth (MHz)", "Center Frequency (GHz)", "Sampling Rate (MHz)", "Sample Count", "Modulation", "Signal Interval (ms)", "Total Signals", "Signal Duration (ms)"],
        value = [params.duration_s*1000, params.bandwidth_hz/1e6, params.center_frequency_hz/1e9, params.sampling_rate_hz/1e6, N, "QPSK", 20.0, signal_count, params.duration_s*1000]
    )
    
    # CSVファイルに保存
    CSV.write("$(output_dir)/continuous_time_domain_data.csv", time_data)
    CSV.write("$(output_dir)/continuous_frequency_domain_data.csv", freq_data)
    CSV.write("$(output_dir)/continuous_signal_parameters.csv", params_data)
    
    return time_data, freq_data, params_data
end

# ===== 継続的信号用PNG出力関数 =====

function plot_continuous_signal_time_domain(time_axis::Vector{Float64}, continuous_signal::Vector{ComplexF64}, title::String)
    # データを間引いてプロット（最大10000点）
    max_points = 10000
    step = max(1, div(length(time_axis), max_points))
    time_subset = time_axis[1:step:end]
    signal_subset = continuous_signal[1:step:end]
    
    p1 = plot(time_subset, real.(signal_subset), label="Real Part", linewidth=1, color=:blue,
              xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Real Part")
    
    p2 = plot(time_subset, imag.(signal_subset), label="Imaginary Part", linewidth=1, color=:red,
              xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Imaginary Part")
    
    p3 = plot(time_subset, abs.(signal_subset), label="Amplitude", linewidth=1, color=:green,
              xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Amplitude")
    
    p4 = plot(time_subset, abs2.(signal_subset), label="Power", linewidth=1, color=:purple,
              xlabel="Time (ms)", ylabel="Power", title="$(title) - Power")
    
    return plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 800))
end

function plot_continuous_signal_frequency_domain(continuous_signal::Vector{ComplexF64}, sampling_rate::Float64, title::String)
    N = length(continuous_signal)
    freq_domain = fft(continuous_signal)
    freq_axis = fftfreq(N, sampling_rate) / 1e6  # MHz単位
    
    # fftshiftを適用
    freq_axis_shifted = fftshift(freq_axis)
    amplitude_spectrum_shifted = fftshift(abs.(freq_domain))
    power_spectrum_shifted = fftshift(abs2.(freq_domain))
    
    p1 = plot(freq_axis_shifted, amplitude_spectrum_shifted, label="Amplitude Spectrum", linewidth=1, color=:blue,
              xlabel="Frequency (MHz)", ylabel="Amplitude", title="$(title) - Amplitude Spectrum")
    
    p2 = plot(freq_axis_shifted, power_spectrum_shifted, label="Power Spectrum", linewidth=1, color=:red,
              xlabel="Frequency (MHz)", ylabel="Power", title="$(title) - Power Spectrum")
    
    return plot(p1, p2, layout=(1,2), size=(1200, 400))
end

function save_continuous_signal_to_png(time_axis::Vector{Float64}, continuous_signal::Vector{ComplexF64}, 
                                     sampling_rate::Float64, output_dir::String)
    gr()
    default(fontfamily="Arial")
    
    title_prefix = "Continuous QPSK Signals (20ms interval)"
    time_plot = plot_continuous_signal_time_domain(time_axis, continuous_signal, title_prefix)
    freq_plot = plot_continuous_signal_frequency_domain(continuous_signal, sampling_rate, title_prefix)
    
    combined_plot = plot(time_plot, freq_plot, layout=grid(2,1,heights=[0.7, 0.3]), size=(1200, 1200))
    savefig(combined_plot, "$(output_dir)/continuous_signal_visualization.png")
    
    return combined_plot
end

# ===== メイン実行関数 =====

function main(output_mode::String="csv")
    println("="^60)
    title_str = "継続的QPSK信号生成プログラム（20ms間隔）"
    if output_mode == "csv"
        println(title_str * " - CSV出力")
    elseif output_mode == "png"
        println(title_str * " - PNG出力")
    else
        println(title_str * " - CSV+PNG出力")
    end
    println("="^60)

    # 信号パラメータ設定
    params = SignalParameters(
        142.8e-6,   # 持続時間: 142.8 µs
        4.7e9,      # 中心周波数
        7.2e6,      # 帯域幅: 7.2 MHz
        15.36e6     # サンプリングレート: 15.36 MHz
    )
    
    println("信号パラメータ:")
    println("• 変調方式: QPSK")
    println("• 持続時間: $(params.duration_s*1000) ms")
    println("• 帯域幅: $(params.bandwidth_hz/1e6) MHz")
    println("• サンプリングレート: $(params.sampling_rate_hz/1e6) MHz")
    println()

    # 乱数シードを固定（再現性のため）
    Random.seed!(1234)
    
    # 継続的信号生成パラメータ
    interval_ms = 20.0  # 20ms間隔
    total_duration_ms = 200.0  # 200ms間のシミュレーション
    
    # 継続的信号生成
    println("継続的QPSK信号を生成中（20ms間隔）...")
    time_axis, continuous_signal, signal_count = generate_continuous_qpsk_signals(
        params, interval_ms, total_duration_ms
    )
    
    println("継続的信号生成完了！")
    println("• 総サンプル数: $(length(continuous_signal))")
    println("• 送信信号数: $(signal_count)")
    println()

    # 信号の統計情報
    signal_power = abs2.(continuous_signal)
    println("継続的信号統計:")
    println("• 平均電力: $(mean(signal_power))")
    println("• 最大電力: $(maximum(signal_power))")
    println("• ピーク対平均電力比 (PAPR): $(maximum(signal_power)/mean(signal_power))")
    println()

    # 出力ディレクトリを作成
    output_dir = "results_continuous_QPSK"
    mkpath(output_dir)
    
    # 出力処理
    if output_mode == "csv" || output_mode == "both"
        println("CSVファイルを生成中...")
        save_continuous_signal_to_csv(time_axis, continuous_signal, params.sampling_rate_hz, params, signal_count, output_dir)
        println("CSVファイルを '$(output_dir)' に保存しました。")
        println("• continuous_time_domain_data.csv - 時間領域データ")
        println("• continuous_frequency_domain_data.csv - 周波数領域データ")
        println("• continuous_signal_parameters.csv - パラメータ情報")
        println()
    end
    
    if output_mode == "png" || output_mode == "both"
        println("PNGファイルを生成中...")
        save_continuous_signal_to_png(time_axis, continuous_signal, params.sampling_rate_hz, output_dir)
        println("PNGファイルを '$(output_dir)' に保存しました。")
        println("• continuous_signal_visualization.png - 継続的信号可視化")
        println()
    end
    
    println("\n" * "="^60)
    println("プログラム完了")
    println("="^60)
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) > 0
        output_mode = ARGS[1]
        if !(output_mode in ["csv", "png", "both"])
            println("エラー: 不正な引数です。'csv', 'png', または 'both' を指定してください。")
            println("デフォルトの 'csv' モードで実行します。")
            output_mode = "csv"
        end
        main(output_mode)
    else
        main("csv")  # デフォルトはCSV出力
    end
end
