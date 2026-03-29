# ===== main_simulation2.jl (最終決定版) =====
# 複数端末・非同期・B案ロジック・衝突判定機能付き
# 評価指標: パケット衝突率 (PER) と 到達率 (PDR)

using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, Distributions, CSV, DataFrames, Plots, Dates

# 各機能モジュール
include("signal_generation.jl")
include("path_loss.jl")
include("shadowing.jl")
include("noise_generation.jl")
include("terminal_deployment.jl")
include("adaptive_threshold.jl")
include("local_clock.jl")
include("sync_observation.jl")
# include("signal_visualization.jl")

# Struct定義
struct SimulationParameters
    signal_duration_us::Float64
    center_frequency_ghz::Float64
    signal_bandwidth_mhz::Float64
    terminal_bandwidth_mhz::Float64
    tx_sampling_rate_mhz::Float64
    rx_sampling_rate_mhz::Float64
    tx_power_dbm::Float64
    shadowing_enabled::Bool
    shadowing_std_db::Float64
    deployment_mode::String
    num_terminals::Int
    area_size_m::Float64
    min_distance_m::Float64
    max_distance_m::Float64
    path_loss_exponent::Float64
    reference_distance_m::Float64
    reference_path_loss_db::Float64
    signal_interval_ms::Float64
    total_duration_ms::Float64
    threshold_margin_db::Float64
    clock_precision::String
    temperature_variation_c::Float64
    sync_accuracy_ns::Float64
    representative_selection::String
    known_beacon_period_ms::Float64
    gate_width_ms::Float64
    slot_length_ms::Float64
    lora_toa_ms::Float64
    transmission_probability::Float64
    num_channels::Int
    initial_window_duration_ms::Float64
    max_startup_delay_ms::Float64
end

# Slot simulationの読み込み
include("slot_simulation.jl") 

# ===== ヘルパー関数群 =====
function estimate_noise_floor_dbm(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}; pre_signal_end_ms::Float64=9.0)
    idxs = findall(t -> t < pre_signal_end_ms, time_axis_ms)
    if isempty(idxs)
        noise_power_w = median(rx_power)
        return 10 * log10(noise_power_w * 1000)
    end
    segment = rx_power[idxs]
    noise_power_w = median(segment)
    return 10 * log10(noise_power_w * 1000)
end

function detect_threshold_crossings(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}, threshold_w::Float64)
    crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[], :end_times => Float64[],
        :peak_times => Float64[], :peak_powers_w => Float64[], :peak_powers_dbm => Float64[]
    )
    is_above_threshold = false
    segment_peak_idx = 0
    segment_peak_power = 0.0
    
    for i in 2:length(rx_power)
        p_prev = rx_power[i-1]; p_curr = rx_power[i]
        t_prev = time_axis_ms[i-1]; t_curr = time_axis_ms[i]
        
        if !is_above_threshold && p_prev <= threshold_w && p_curr > threshold_w
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            push!(crossings[:start_times], crossing_time)
            is_above_threshold = true
            segment_peak_idx = i
            segment_peak_power = p_curr
        end
        
        if is_above_threshold
            if p_curr > segment_peak_power
                segment_peak_idx = i
                segment_peak_power = p_curr
            end
        end
        
        if is_above_threshold && p_prev > threshold_w && p_curr <= threshold_w
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            push!(crossings[:end_times], crossing_time)
            peak_time = time_axis_ms[segment_peak_idx]
            push!(crossings[:peak_times], peak_time)
            push!(crossings[:peak_powers_w], segment_peak_power)
            push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
            is_above_threshold = false
        end
    end
    if is_above_threshold && !isempty(time_axis_ms)
        push!(crossings[:end_times], time_axis_ms[end])
        peak_time = time_axis_ms[segment_peak_idx]
        push!(crossings[:peak_times], peak_time)
        push!(crossings[:peak_powers_w], segment_peak_power)
        push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
    end
    return crossings
end

function filter_by_debounce(crossings::Dict{Symbol, Vector{Float64}}, debounce_time_ms::Float64)
    if isempty(crossings[:peak_times]) return crossings end
    filtered_crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[], :end_times => Float64[],
        :peak_times => Float64[], :peak_powers_w => Float64[], :peak_powers_dbm => Float64[]
    )
    last_accepted_peak_time_ms = -Inf
    for i in 1:length(crossings[:peak_times])
        current_peak_time = crossings[:peak_times][i]
        if current_peak_time >= last_accepted_peak_time_ms + debounce_time_ms
            push!(filtered_crossings[:start_times], crossings[:start_times][i])
            push!(filtered_crossings[:end_times], crossings[:end_times][i])
            push!(filtered_crossings[:peak_times], current_peak_time)
            push!(filtered_crossings[:peak_powers_w], crossings[:peak_powers_w][i])
            push!(filtered_crossings[:peak_powers_dbm], crossings[:peak_powers_dbm][i])
            last_accepted_peak_time_ms = current_peak_time
        end
    end
    return filtered_crossings
end

function add_path_loss_and_noise(tx_signal::Vector{ComplexF64}, path_loss_params::PathLossParameters, 
                                 noise_params::NoiseParameters, shadowing_params::ShadowingParameters, sampling_rate_hz::Float64, 
                                terminal_shadowing_db::Float64)
    path_loss_db = calculate_path_loss(path_loss_params)
    path_loss_linear = 10^(-path_loss_db / 10)
    shadowing_db = terminal_shadowing_db
    shadowing_linear = 10^(-shadowing_db / 10)
    total_loss_linear = path_loss_linear * shadowing_linear
    rx_signal_with_path_loss = tx_signal * sqrt(total_loss_linear)
    signal_power_linear = mean(abs2.(rx_signal_with_path_loss))
    noise_power_dbm = noise_params.noise_power_dbm
    noise = generate_awgn_noise(noise_power_dbm, length(tx_signal))
    rx_signal_with_noise = rx_signal_with_path_loss + noise
    actual_noise_power = mean(abs2.(noise))
    actual_snr_db = 10 * log10(signal_power_linear / actual_noise_power)
    return rx_signal_with_noise, actual_snr_db, path_loss_db
end

function filter_by_sample_count(crossings_raw::Dict{Symbol, Vector{Float64}}, rx_sampling_rate_mhz::Float64, min_samples::Int)
    num_crossings_raw = length(crossings_raw[:start_times])
    if num_crossings_raw == 0 return crossings_raw end
    sampling_interval_ms = 1.0 / (rx_sampling_rate_mhz * 1e6) * 1000
    filtered_indices = Int[]
    for i in 1:num_crossings_raw
        duration_ms = crossings_raw[:end_times][i] - crossings_raw[:start_times][i]
        num_samples = ceil(duration_ms / sampling_interval_ms)
        if num_samples >= min_samples push!(filtered_indices, i) end
    end
    if !isempty(filtered_indices)
        return Dict{Symbol, Vector{Float64}}(
            :start_times => crossings_raw[:start_times][filtered_indices],
            :end_times => crossings_raw[:end_times][filtered_indices],
            :peak_times => crossings_raw[:peak_times][filtered_indices],
            :peak_powers_w => crossings_raw[:peak_powers_w][filtered_indices],
            :peak_powers_dbm => crossings_raw[:peak_powers_dbm][filtered_indices]
        )
    else
        return Dict{Symbol, Vector{Float64}}(:start_times => Float64[], :end_times => Float64[], :peak_times => Float64[], :peak_powers_w => Float64[], :peak_powers_dbm => Float64[])
    end
end

# ==================================================
# ★★★ 衝突解析・可視化関数 (強化版) ★★★
# ==================================================
function analyze_results(df::DataFrame, total_duration_ms::Float64, num_terminals::Int, output_dir::String)
    println("\n" * "="^60)
    println("シミュレーション結果解析 (判定・詳細付与・可視化)")
    println("="^60)
    
    total_packets = nrow(df)
    if total_packets == 0
        println("送信パケットがありません。")
        return df
    end

    # 1. 結果を記録する新しい列を追加
    df[!, :status] .= "Success"           # デフォルトは成功
    df[!, :collision_partners] .= ""      # 衝突相手のID (カンマ区切り文字列)
    
    # 判定用フラグ配列 (trueなら衝突)
    is_collided = falses(total_packets)
    
    # 時間順にソート (解析の前提)
    sort!(df, :actual_tx_start_global_ms)

    # === 総当たり衝突判定 ===
    for i in 1:total_packets
        for j in (i+1):total_packets
            # 時間最適化: Jの開始がIの終了を超えていれば、それ以降は衝突しない
            if df[j, :actual_tx_start_global_ms] >= df[i, :actual_tx_end_global_ms]
                break
            end

            # チャネルチェック
            if df[i, :channel_id] != df[j, :channel_id]
                continue
            end
            
            # 時間重複チェック
            start_i = df[i, :actual_tx_start_global_ms]
            end_i   = df[i, :actual_tx_end_global_ms]
            start_j = df[j, :actual_tx_start_global_ms]
            end_j   = df[j, :actual_tx_end_global_ms]
            
            if max(start_i, start_j) < min(end_i, end_j)
                # 衝突確定
                is_collided[i] = true
                is_collided[j] = true
                
                # ステータス更新
                df[i, :status] = "Collision"
                df[j, :status] = "Collision"
                
                # 相手のIDを記録 (既に書かれていれば追記)
                partner_i = string(df[i, :terminal_id])
                partner_j = string(df[j, :terminal_id])
                
                df[i, :collision_partners] = isempty(df[i, :collision_partners]) ? partner_j : df[i, :collision_partners] * "," * partner_j
                df[j, :collision_partners] = isempty(df[j, :collision_partners]) ? partner_i : df[j, :collision_partners] * "," * partner_i
            end
        end
    end
    
    # === 数値集計 ===
    collision_count = count(is_collided)
    success_count = total_packets - collision_count
    per = collision_count / total_packets
    #pdr = success_count / total_packets
    
    println("解析結果:")
    println("• 総パケット数:     $total_packets")
    println("• 成功パケット数:   $success_count")
    println("• 衝突パケット数:   $collision_count")
    println("-"^30)
    println("• パケット衝突率 (PER):     $(round(per * 100, digits=2)) %")
    #println("• パケット到達率 (PDR):     $(round(pdr * 100, digits=2)) %")
    
    # === ★★★ 可視化 (チャネル使用マップ) ★★★ ===
    println("\nグラフを生成中...")
    
    # プロットの準備
    p = plot(
        title = "Channel Occupancy & Collisions",
        xlabel = "Time (ms)",
        ylabel = "Channel ID",
        legend = :topright,
        size = (1000, 600),
        xlims = (0, total_duration_ms),
        ylims = (0, 9), # チャネル1-8を表示
        yticks = 1:8
    )
    
    # 成功パケットと衝突パケットを分けて描画
    success_df = filter(row -> row.status == "Success", df)
    collision_df = filter(row -> row.status == "Collision", df)
    
    # 成功パケット (青色)
    for row in eachrow(success_df)
        # 矩形を描画: shape(x座標配列, y座標配列)
        # チャネルIDを中心に ±0.4 の幅を持たせる
        t_start = row.actual_tx_start_global_ms
        t_end = row.actual_tx_end_global_ms
        ch = row.channel_id
        plot!(p, [t_start, t_end, t_end, t_start], [ch-0.4, ch-0.4, ch+0.4, ch+0.4], 
              seriestype=:shape, color=:blue, opacity=0.5, linecolor=:white, label="")
    end
    
    # 衝突パケット (赤色)
    for row in eachrow(collision_df)
        t_start = row.actual_tx_start_global_ms
        t_end = row.actual_tx_end_global_ms
        ch = row.channel_id
        plot!(p, [t_start, t_end, t_end, t_start], [ch-0.4, ch-0.4, ch+0.4, ch+0.4], 
              seriestype=:shape, color=:red, opacity=0.6, linecolor=:white, label="")
    end
    
    # 凡例用のダミープロット
    plot!(p, [], [], seriestype=:shape, color=:blue, label="Success")
    plot!(p, [], [], seriestype=:shape, color=:red, label="Collision")

    # 保存
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    png_filename = "$(output_dir)/channel_map_$(timestamp).png"
    savefig(p, png_filename)
    println("• グラフ保存完了: $png_filename")
    
    println("="^60)
    
    return df # 情報を付与したDataFrameを返す
end

# ===== シミュレーションパラメータ作成 =====
function create_simulation_parameters()
    return SimulationParameters(
        66.67, 0.92, 3.6, 0.125, 7.68, 0.25, 43.0,
        true, 0.0,
        "random_fixed",
        2,                 # 端末数
        100.0, 10.0, 500.0, 3.0, 1.0, 0.0,
        20.0,
        5000.0,             # ★ 総持続時間
        15.0,
        "low", 5.0, 100.0, "median",
        20.0, 2.0,
        500.0, 288.8,
        0.5,                # 送信確率
        1,                  # チャネル数
        110.0,              # 初期同期窓
        2000.0              # 最大起動遅延
    )
end


# ===== メイン実行関数 =====
function main_simulation2()
    Random.seed!(time_ns())
    
    println("="^60)
    println("複数端末シミュレーション実行 (完全統合版)")
    println("="^60)

    # 1. パラメータ
    sim_params = create_simulation_parameters()
    
    # 2. ノイズ・シャドウイング設定
    noise_figure_db = 5.0
    terminal_bandwidth_hz = sim_params.terminal_bandwidth_mhz * 1e6
    fixed_noise_power_dbm = -174 + 10 * log10(terminal_bandwidth_hz) + noise_figure_db
    noise_params = NoiseParameters(fixed_noise_power_dbm, 100.0)
    shadowing_params = ShadowingParameters(sim_params.shadowing_enabled, sim_params.shadowing_std_db, 50.0, 0.5)
    
    # 3. 端末配置
    println("端末配置を生成中...")
    deployment_params = TerminalDeploymentParameters(
        sim_params.deployment_mode, 0.001, sim_params.num_terminals, sim_params.area_size_m,
        sim_params.min_distance_m, sim_params.max_distance_m, sim_params.center_frequency_ghz * 1e9,
        sim_params.path_loss_exponent, sim_params.reference_distance_m, sim_params.reference_path_loss_db
    )
    terminals = deploy_terminals(deployment_params, shadowing_params.std_db, shadowing_params.enabled, sim_params.tx_power_dbm)
    
    # 4. 信号生成・ダウンサンプリング
    params_tx = SignalParameters(
        sim_params.signal_duration_us * 1e-6, sim_params.center_frequency_ghz * 1e9,
        sim_params.signal_bandwidth_mhz * 1e6, sim_params.tx_sampling_rate_mhz * 1e6, sim_params.tx_power_dbm
    )
    println("送信信号生成・ダウンサンプリング中...")
    time_axis_tx, tx_signal_high_rate, signal_count = generate_periodic_sync_signals(
        params_tx, sim_params.signal_interval_ms, sim_params.total_duration_ms
    )
    rate_ratio = sim_params.rx_sampling_rate_mhz / sim_params.tx_sampling_rate_mhz
    tx_signal_low_rate = resample(tx_signal_high_rate, rate_ratio)
    time_axis_rx = collect((0:length(tx_signal_low_rate)-1) / (sim_params.rx_sampling_rate_mhz * 1e6) * 1000)

    # 集約用データフレーム
    all_tx_events = DataFrame()
    all_terminal_info = DataFrame()

    println("\n各端末の処理を開始します...")
    
    # 5. 端末ループ
    for (idx, terminal) in enumerate(terminals)
        print("\r処理中: 端末 $idx / $(length(terminals)) ... ")

        # クロック初期化
        local terminal_clock::LocalClock
        ref_freq_hz = sim_params.center_frequency_ghz * 1e9
        if sim_params.clock_precision == "high"
            clock_params = create_high_precision_clock_parameters()
        elseif sim_params.clock_precision == "low"
            clock_params = create_low_precision_clock_parameters()
        elseif sim_params.clock_precision == "none"
            clock_params = create_no_drift_clock_parameters()
        else
            clock_params = create_default_clock_parameters()
        end
        terminal_clock = initialize_local_clock(terminal.terminal_id, ref_freq_hz, clock_params)

        # 受信シミュレーション
        terminal_path_loss_params = PathLossParameters(
            terminal.distance_m, sim_params.center_frequency_ghz * 1e9, sim_params.path_loss_exponent,
            sim_params.reference_distance_m, sim_params.reference_path_loss_db
        )
        rx_signal, actual_snr_db, path_loss_db = add_path_loss_and_noise(
            tx_signal_low_rate, terminal_path_loss_params, noise_params, shadowing_params, 
            sim_params.rx_sampling_rate_mhz * 1e6, terminal.shadowing_db
        )
        rx_power = abs2.(rx_signal)

        # ノイズフロア推定
        noise_floor_dbm = estimate_noise_floor_dbm(rx_power, time_axis_rx; pre_signal_end_ms=10.0)
        threshold_dbm = noise_floor_dbm + sim_params.threshold_margin_db
        threshold_w = 10^(threshold_dbm / 10) * 1e-3

        # 非同期起動設定
        startup_time_global_ms = rand() * sim_params.max_startup_delay_ms
        window_end_global_ms = startup_time_global_ms + sim_params.initial_window_duration_ms
        
        if idx == 1
            println("\n[端末1 詳細] SNR: $(round(actual_snr_db, digits=1)) dB, 起動: $(round(startup_time_global_ms, digits=1)) ms")
        end

        # ビーコン検出
        crossings_raw = detect_threshold_crossings(rx_power, time_axis_rx, threshold_w)
        crossings_filtered = filter_by_sample_count(crossings_raw, sim_params.rx_sampling_rate_mhz, 2)
        crossings_initial = filter_by_debounce(crossings_filtered, 1.0)

        valid_crossings = Dict{Symbol, Vector{Float64}}(:peak_times => Float64[], :local_peak_times => Float64[])
        
        for i in 1:length(crossings_initial[:peak_times])
            g_time = crossings_initial[:peak_times][i]
            if startup_time_global_ms <= g_time <= window_end_global_ms
                push!(valid_crossings[:peak_times], g_time)
                g_time_s = g_time / 1000.0
                l_time_s = convert_to_local_time(terminal_clock, g_time_s)
                push!(valid_crossings[:local_peak_times], l_time_s * 1000.0)
            end
        end

        num_valid = length(valid_crossings[:peak_times])

        # 送信シミュレーション
        local tx_events::DataFrame
        if num_valid > 0
            window_end_global_s = window_end_global_ms / 1000.0
            window_end_local_s = convert_to_local_time(terminal_clock, window_end_global_s)
            
            tx_events = simulate_slot_transmission(
                valid_crossings, sim_params, terminal_clock, window_end_local_s * 1000.0
            )
            
            if !isempty(tx_events)
                tx_events[!, :terminal_id] .= terminal.terminal_id
                append!(all_tx_events, tx_events)
            end
        end

        # 端末情報保存
        final_clock_stats = evaluate_clock_quality(terminal_clock, sim_params.total_duration_ms / 1000.0)
        this_info = DataFrame(
            terminal_id = [terminal.terminal_id],
            distance_m = [terminal.distance_m],
            startup_time_ms = [startup_time_global_ms],
            num_detected_beacons = [num_valid],
            clock_ppm = [final_clock_stats["frequency_error_ppm"]],
            # ★★★ 追加項目 ★★★
            rx_power_dbm = [terminal.rx_power_dbm],  # 受信電力 (dBm)
            snr_db = [actual_snr_db],                # SNR (dB)
            noise_floor_dbm = [noise_floor_dbm],     # 推定ノイズフロア (dBm)

        )
        append!(all_terminal_info, this_info)
    end 
    println("\n全端末の処理が完了しました。")

    # ===== 6. 結果解析と保存 =====
    output_dir = "results_multi_terminal"
    mkpath(output_dir)
    execution_timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

    println("\n=== 解析結果と保存 ===")

    if !isempty(all_tx_events)
        # ★★★ 修正点: 戻り値(df)を受け取り、output_dir を渡す ★★★
        analyzed_df = analyze_results(all_tx_events, sim_params.total_duration_ms, sim_params.num_terminals, output_dir)
        
        # 解析済みのデータ(status付き)を保存
        filename = "$(output_dir)/all_tx_events_$(execution_timestamp).csv"
        CSV.write(filename, analyzed_df)
        println("• 保存完了: $filename")
    else
        println("警告: 送信イベントが1つもありませんでした")
    end

    filename_info = "$(output_dir)/all_terminal_info_$(execution_timestamp).csv"
    CSV.write(filename_info, all_terminal_info)
    println("• 保存完了: $filename_info")

    println("="^60)
    println("プログラム完了")
    println("="^60)
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    main_simulation2()
end