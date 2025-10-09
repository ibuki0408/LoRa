using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots, Dates

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
    
    # バウンドチェック
    start_bin = max(1, start_bin)
    end_bin = min(N, end_bin)
    
    for i in start_bin:end_bin
        # 各サブキャリアに、ランダムなQPSKシンボルを1つ配置する
        freq_domain_signal[i] = generate_qpsk_symbol()
    end
    
    # IFFTを適用して時間領域の信号に変換
    time_domain_signal_complex = ifft(ifftshift(freq_domain_signal))
    
    return time_domain_signal_complex
end

# ===== CSV出力関数 =====

function save_signal_to_csv(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, params::SignalParameters, output_dir::String)
    N = length(signal_complex)
    time_axis = (0:N-1) / sampling_rate * 1000  # ms単位
    freq_axis = fftfreq(N, sampling_rate) / 1e6  # MHz単位
    
    # 実行時刻を取得
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    # 時間領域データ
    time_data = DataFrame(
        time_ms = time_axis,
        real_part = real.(signal_complex),
        imaginary_part = imag.(signal_complex),
        amplitude = abs.(signal_complex),
        power = abs2.(signal_complex)
    )
    
    # 周波数領域データ
    freq_domain = fft(signal_complex)
    freq_data = DataFrame(
        frequency_mhz = fftshift(freq_axis), # fftshiftを適用して中央揃え
        amplitude_spectrum = fftshift(abs.(freq_domain)),
        power_spectrum = fftshift(abs2.(freq_domain))
    )
    
    # パラメータ情報
    params_data = DataFrame(
        parameter = ["Duration (ms)", "Bandwidth (MHz)", "Center Frequency (GHz)", "Sampling Rate (MHz)", "Sample Count", "Modulation", "Execution Time"],
        value = [params.duration_s*1000, params.bandwidth_hz/1e6, params.center_frequency_hz/1e9, params.sampling_rate_hz/1e6, N, "QPSK", timestamp]
    )
    
    # CSVファイルに保存（タイムスタンプ付き）
    CSV.write("$(output_dir)/time_domain_data_qpsk_$(timestamp).csv", time_data)
    CSV.write("$(output_dir)/frequency_domain_data_qpsk_$(timestamp).csv", freq_data)
    CSV.write("$(output_dir)/signal_parameters_qpsk_$(timestamp).csv", params_data)
    
    return time_data, freq_data, params_data
end

# ===== PNG出力関数 =====

function plot_signal_time_domain(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, title::String)
    N = length(signal_complex)
    time_axis = (0:N-1) / sampling_rate * 1000  # ms単位
    
    p1 = plot(time_axis, real.(signal_complex), label="Real Part", xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Real Part")
    p2 = plot(time_axis, imag.(signal_complex), label="Imaginary Part", xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Imaginary Part")
    p3 = plot(time_axis, abs.(signal_complex), label="Amplitude", xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Amplitude")
    p4 = plot(time_axis, abs2.(signal_complex), label="Power", xlabel="Time (ms)", ylabel="Power", title="$(title) - Power")
    
    return plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 750))
end

function plot_signal_frequency_domain(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, title::String)
    N = length(signal_complex)
    freq_domain = fft(signal_complex)
    freq_axis = fftfreq(N, sampling_rate) / 1e6  # MHz単位
    
    # fftshiftを適用してデータを中央揃えにする
    freq_axis_shifted = fftshift(freq_axis)
    amplitude_spectrum_shifted = fftshift(abs.(freq_domain))
    power_spectrum_shifted = fftshift(abs2.(freq_domain))
    
    p1 = plot(freq_axis_shifted, amplitude_spectrum_shifted, label="Amplitude Spectrum", xlabel="Frequency (MHz)", ylabel="Amplitude", title="$(title) - Amplitude Spectrum")
    p2 = plot(freq_axis_shifted, power_spectrum_shifted, label="Power Spectrum", xlabel="Frequency (MHz)", ylabel="Power", title="$(title) - Power Spectrum")
    
    return plot(p1, p2, layout=(1,2), size=(1000, 375))
end

function plot_signal_comprehensive(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, title_prefix::String)
    time_plot = plot_signal_time_domain(signal_complex, sampling_rate, title_prefix)
    freq_plot = plot_signal_frequency_domain(signal_complex, sampling_rate, title_prefix)
    
    return plot(time_plot, freq_plot, layout=grid(2,1,heights=[0.66, 0.34]), size=(1000, 1125))
end

function save_signal_to_png(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, output_dir::String)
    gr()
    default(fontfamily="Arial")
    
    # 実行時刻を取得
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    
    title_prefix = "OFDM-like Signal (QPSK Modulation)"
    combined_plot = plot_signal_comprehensive(signal_complex, sampling_rate, title_prefix)
    
    savefig(combined_plot, "$(output_dir)/signal_visualization_qpsk_$(timestamp).png")
    
    return combined_plot
end


# ===== メイン実行関数 =====

function main(output_mode::String="csv")
    println("="^60)
    title_str = "OFDMライク信号生成プログラム（QPSK変調）"
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
        1.0e6,      # 帯域幅: 3.6 MHz（7.2MHz→3.6MHzに半減）
        2.0e6      # サンプリングレート: 10 MHz（帯域幅の約2.8倍）
    )
    
    println("信号パラメータ:")
    println("• 変調方式: QPSK")
    println("• 持続時間: $(params.duration_s*1000) ms")
    println("• 帯域幅: $(params.bandwidth_hz/1e6) MHz")
    println("• サンプリングレート: $(params.sampling_rate_hz/1e6) MHz")
    println()

    # 乱数シードを固定（再現性のため）
    Random.seed!(1234)
    
    # 信号生成
    println("OFDMライク信号を生成中...")
    signal_complex = generate_ofdm_signal_with_qpsk(
        params.bandwidth_hz,
        params.duration_s,
        params.sampling_rate_hz
    )
    
    println("信号生成完了！")
    println("• サンプル数: $(length(signal_complex))")
    println()

    # 信号の統計情報
    signal_power = abs2.(signal_complex)
    println("信号統計:")
    println("• 平均電力: $(mean(signal_power))")
    println("• ピーク対平均電力比 (PAPR): $(maximum(signal_power)/mean(signal_power))")
    println()

    # 出力ディレクトリを作成
    output_dir = "results_QPSK"
    mkpath(output_dir)
    
    # 実行時刻を取得
    execution_timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    println("実行時刻: $(execution_timestamp)")
    println()
    
    # 出力処理
    if output_mode == "csv" || output_mode == "both"
        println("CSVファイルを生成中...")
        save_signal_to_csv(signal_complex, params.sampling_rate_hz, params, output_dir)
        println("CSVファイルを '$(output_dir)' に保存しました。")
        println("• time_domain_data_qpsk_$(execution_timestamp).csv")
        println("• frequency_domain_data_qpsk_$(execution_timestamp).csv")
        println("• signal_parameters_qpsk_$(execution_timestamp).csv")
        println()
    end
    
    if output_mode == "png" || output_mode == "both"
        println("PNGファイルを生成中...")
        save_signal_to_png(signal_complex, params.sampling_rate_hz, output_dir)
        println("PNGファイルを '$(output_dir)' に保存しました。")
        println("• signal_visualization_qpsk_$(execution_timestamp).png")
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