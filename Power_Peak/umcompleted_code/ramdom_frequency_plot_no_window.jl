using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots

# ===== 信号パラメータ =====
struct SignalParameters
    duration_s::Float64       # 信号持続時間（s）
    center_frequency_hz::Float64 # 中心周波数（Hz）
    bandwidth_hz::Float64     # 帯域幅（Hz）
    sampling_rate_hz::Float64 # サンプリングレート（Hz）
end

# ===== 信号生成関数 (ハニング窓なし) =====

"""
周波数領域で定義したランダムな広帯域信号（OFDMライクな信号）を生成する
ハニング窓を適用しないバージョン
"""
function generate_random_ofdm_signal_no_window(bandwidth_hz::Float64, duration_s::Float64, fs_hz::Float64)
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
        # 窓関数を適用せずにランダムな複素数を生成
        freq_domain_signal[i] = rand(ComplexF64)
    end
    
    # IFFTを適用して時間領域の信号に変換
    # fftshiftは直流成分(0Hz)を配列の中心から配列の先頭に戻すための操作
    time_domain_signal_complex = ifft(ifftshift(freq_domain_signal))
    
    return time_domain_signal_complex
end


# ===== CSV出力関数 =====

function save_signal_to_csv(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, params::SignalParameters, output_dir::String)
    N = length(signal_complex)
    time_axis = (0:N-1) / sampling_rate * 1000  # ms単位
    freq_axis = fftfreq(N, sampling_rate) / 1e6  # MHz単位
    
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
        frequency_mhz = freq_axis,
        amplitude_spectrum = abs.(freq_domain),
        power_spectrum = abs2.(freq_domain)
    )
    
    # パラメータ情報
    params_data = DataFrame(
        parameter = ["Duration (ms)", "Bandwidth (MHz)", "Center Frequency (GHz)", "Sampling Rate (MHz)", "Sample Count"],
        value = [params.duration_s*1000, params.bandwidth_hz/1e6, params.center_frequency_hz/1e9, params.sampling_rate_hz/1e6, N]
    )
    
    # CSVファイルに保存
    CSV.write("$(output_dir)/time_domain_data.csv", time_data)
    CSV.write("$(output_dir)/frequency_domain_data.csv", freq_data)
    CSV.write("$(output_dir)/signal_parameters.csv", params_data)
    
    return time_data, freq_data, params_data
end

# ===== PNG出力関数 =====

function plot_signal_time_domain(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, title::String)
    N = length(signal_complex)
    time_axis = (0:N-1) / sampling_rate * 1000  # ms単位
    
    # 実部と虚部をプロット
    p1 = plot(time_axis, real.(signal_complex), 
              label="Real Part", linewidth=2, color=:blue,
              xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Real Part")
    
    p2 = plot(time_axis, imag.(signal_complex), 
              label="Imaginary Part", linewidth=2, color=:red,
              xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Imaginary Part")
    
    # 振幅（絶対値）をプロット
    p3 = plot(time_axis, abs.(signal_complex), 
              label="Amplitude", linewidth=2, color=:green,
              xlabel="Time (ms)", ylabel="Amplitude", title="$(title) - Amplitude")
    
    # 電力（絶対値の2乗）をプロット
    p4 = plot(time_axis, abs2.(signal_complex), 
              label="Power", linewidth=2, color=:purple,
              xlabel="Time (ms)", ylabel="Power", title="$(title) - Power")
    
    return plot(p1, p2, p3, p4, layout=(2,2), size=(800, 600))
end

function plot_signal_frequency_domain(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, title::String)
    N = length(signal_complex)
    
    # FFTを計算
    freq_domain = fft(signal_complex)
    freq_axis = fftfreq(N, sampling_rate) / 1e6  # MHz単位
    
    # 周波数スペクトラム（振幅）
    p1 = plot(freq_axis, abs.(freq_domain), 
              label="Amplitude Spectrum", linewidth=2, color=:blue,
              xlabel="Frequency (MHz)", ylabel="Amplitude", title="$(title) - Amplitude Spectrum")
    
    # 周波数スペクトラム（電力）
    p2 = plot(freq_axis, abs2.(freq_domain), 
              label="Power Spectrum", linewidth=2, color=:red,
              xlabel="Frequency (MHz)", ylabel="Power", title="$(title) - Power Spectrum")
    
    return plot(p1, p2, layout=(1,2), size=(800, 300))
end

function plot_signal_comprehensive(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, params::SignalParameters)
    # 時間領域プロット
    time_plot = plot_signal_time_domain(signal_complex, sampling_rate, "OFDM-like Signal (No Window)")
    
    # 周波数領域プロット
    freq_plot = plot_signal_frequency_domain(signal_complex, sampling_rate, "OFDM-like Signal (No Window)")
    
    # 全体を結合
    combined_plot = plot(time_plot, freq_plot, layout=(2,1), size=(800, 900))
    
    return combined_plot
end

function save_signal_to_png(signal_complex::Vector{ComplexF64}, sampling_rate::Float64, params::SignalParameters, output_dir::String)
    # フォント設定
    gr()  # GRバックエンドを使用
    default(fontfamily="Arial")
    
    # 包括的可視化
    combined_plot = plot_signal_comprehensive(signal_complex, sampling_rate, params)
    
    # PNGファイルに保存
    savefig(combined_plot, "$(output_dir)/signal_visualization_no_window.png")
    
    return combined_plot
end


# ===== メイン実行関数 =====

function main(output_mode::String="csv")
    println("="^60)
    if output_mode == "csv"
        println("OFDMライク信号生成・CSV出力プログラム（ハニング窓なし）")
    elseif output_mode == "png"
        println("OFDMライク信号生成・PNG出力プログラム（ハニング窓なし）")
    else
        println("OFDMライク信号生成・CSV+PNG出力プログラム（ハニング窓なし）")
    end
    println("="^60)

    # 信号パラメータ設定
    params = SignalParameters(
        142.8e-6,   # 持続時間: 142.8 µs (4シンボル分)
        4.7e9,      # 中心周波数 (変更なし)
        7.2e6,      # 帯域幅: 7.2 MHz (240サブキャリア × 30kHz)
        15.36e6     # サンプリングレート: 15.36 MHz (5Gの標準値)
    )
    
    println("信号パラメータ:")
    println("• 持続時間: $(params.duration_s*1000) ms")
    println("• 帯域幅: $(params.bandwidth_hz/1e6) MHz")
    println("• 中心周波数: $(params.center_frequency_hz/1e9) GHz")
    println("• サンプリングレート: $(params.sampling_rate_hz/1e6) MHz")
    println()

    # 乱数シードを固定（再現性のため）
    rng = MersenneTwister(1234)
    
    # 信号生成
    println("OFDMライク信号を生成中（ハニング窓なし）...")
    signal_complex = generate_random_ofdm_signal_no_window(
        params.bandwidth_hz,
        params.duration_s,
        params.sampling_rate_hz
    )
    
    println("信号生成完了！")
    println("• サンプル数: $(length(signal_complex))")
    println("• 信号長: $(length(signal_complex)/params.sampling_rate_hz*1000) ms")
    println()

    # 信号の統計情報
    signal_power = abs2.(signal_complex)
    signal_amplitude = abs.(signal_complex)
    
    println("信号統計:")
    println("• 平均電力: $(mean(signal_power))")
    println("• 最大振幅: $(maximum(signal_amplitude))")
    println("• 平均振幅: $(mean(signal_amplitude))")
    println("• ピーク対平均比: $(maximum(signal_power)/mean(signal_power))")
    println("• 開始時電力: $(signal_power[1])")
    println("• 終了時電力: $(signal_power[end])")
    println()

    # 出力ディレクトリを作成
    output_dir = "results_no_window"
    mkpath(output_dir)
    
    # 出力処理
    if output_mode == "csv" || output_mode == "both"
        println("CSVファイルを生成中...")
        time_data, freq_data, params_data = save_signal_to_csv(signal_complex, params.sampling_rate_hz, params, output_dir)
        
        println("CSVファイルを保存しました:")
        println("• $(output_dir)/time_domain_data.csv - 時間領域データ")
        println("• $(output_dir)/frequency_domain_data.csv - 周波数領域データ")
        println("• $(output_dir)/signal_parameters.csv - パラメータ情報")
        println()
        
        # データの概要を表示
        println("データ概要:")
        println("• 時間領域データ: $(nrow(time_data)) サンプル")
        println("• 周波数領域データ: $(nrow(freq_data)) サンプル")
        println("• パラメータ数: $(nrow(params_data)) 項目")
        println()
    end
    
    if output_mode == "png" || output_mode == "both"
        println("PNGファイルを生成中...")
        combined_plot = save_signal_to_png(signal_complex, params.sampling_rate_hz, params, output_dir)
        
        println("PNGファイルを保存しました:")
        println("• $(output_dir)/signal_visualization_no_window.png - 信号可視化")
        println()
    end
    
    println("\n" * "="^60)
    println("プログラム完了")
    println("="^60)
end




# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    # コマンドライン引数をチェック
    if length(ARGS) > 0
        output_mode = ARGS[1]
        if output_mode in ["csv", "png", "both"]
            main(output_mode)
        else
            println("使用方法: julia ramdom_frequency_plot_no_window.jl [csv|png|both]")
            println("• csv: CSVファイルのみ出力")
            println("• png: PNGファイルのみ出力") 
            println("• both: CSVとPNG両方出力")
            println("引数なし: CSV出力（デフォルト）")
            main("csv")
        end
    else
        main("csv")  # デフォルトはCSV出力
    end
end
