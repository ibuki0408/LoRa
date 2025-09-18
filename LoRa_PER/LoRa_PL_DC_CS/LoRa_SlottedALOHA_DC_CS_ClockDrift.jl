using Random, FFTW, Statistics, Printf, DelimitedFiles, LinearAlgebra, StatsBase

# ===== パラメータ設定 =====
# エリアサイズ(km)
const area_size = 0.5

# パスロス係数
const α = 4.0
const β = 9.5
const γ = 4.5

# 搬送波周波数(MHz)
const f_c = 923.2

# シャドウイング標準偏差[dB]
const shadowing_std = 3.48

# 送信電力(dBm)
const Tx_dB = 13.0

# 雑音指数[dB]
const noise_figure = 10

# 帯域幅(Hz)
const band_width = 125e3

# 雑音スペクトラム密度(dBm/Hz)
const noise_power_spectrum_density = -174

# SF値（固定）
const SF = 7

# SNR閾値（SF7用）
const SNR_threshold = -7.5

# キャリアセンス閾値(dBm)
const carrier_sense_threshold = -100.0

# ===== クロックドリフト関連パラメータ =====
# クロックドリフトの範囲（ppm）
const CLOCK_DRIFT_MIN_PPM = -50.0  # 低コスト水晶発振器
const CLOCK_DRIFT_MAX_PPM = 50.0   # 低コスト水晶発振器
const CLOCK_DRIFT_HIGH_ACCURACY_MIN_PPM = -5.0  # 高精度水晶発振器
const CLOCK_DRIFT_HIGH_ACCURACY_MAX_PPM = 5.0   # 高精度水晶発振器

# ===== LoRa 基本関数 =====
@inline function lora_symbol(sf::Int, bw::Float64, m::Int)
    M  = 2^sf
    tc = 1.0 / bw
    ts = M * tc
    x = Vector{ComplexF64}(undef, M)
    @inbounds for k in 1:M
        x[k] = exp(1im * 2π * (bw/(2*ts)) * ((k-1) + m)^2 * tc^2)
    end
    return x
end

const CHIRP_COMPENSATION_CACHE = Dict{Tuple{Int,Float64}, Vector{ComplexF64}}()

function get_chirp_compensation(sf::Int, bw::Float64)
    return get!(CHIRP_COMPENSATION_CACHE, (sf,bw)) do
        M = 2^sf
        tc = 1.0 / bw
        ts = tc * M
        exp.(-1im * 2π * (0:M-1).^2 * bw / (2*ts) * tc^2)
    end
end

# 復調
function demod_lora(sf::Int, bw::Float64, x::AbstractVector{ComplexF64})
    comp_coeff = get_chirp_compensation(sf, bw)
    d = x .* comp_coeff
    Y = fft(d)
    return argmax(abs.(Y)) - 1
end

# AWGN付加
function add_awgn!(y::AbstractVector{ComplexF64}, snr_dB::Float64)
    var = 10^(-snr_dB/10)
    σ = sqrt(var/2)
    @inbounds for i in eachindex(y)
        y[i] += ComplexF64(randn()*σ, randn()*σ)
    end
    return y
end

# ===== 伝搬モデル関数 =====
function distance(x1, y1, x2, y2)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

function path_loss(distance_km::Float64)
    return 10 * α * log10(distance_km) + β + 10 * γ * log10(f_c)
end

# シャドウイング値生成関数
function shadowing_value(rng)
    return randn(rng) * shadowing_std
end

# 受信電力計算関数
function received_power(tx_power_dBm::Float64, distance_km::Float64, shadowing_dB::Float64)
    pl = path_loss(distance_km)
    return tx_power_dBm - pl - shadowing_dB
end

# SNR計算関数
function calculate_snr(distance_km::Float64, shadowing_dB::Float64)
    received_pwr = received_power(Tx_dB, distance_km, shadowing_dB)
    snr = received_pwr - (noise_power_spectrum_density + 10*log10(band_width) + noise_figure)
    return snr
end

# dBm/mW 変換
@inline function dBm_to_mW(p_dBm::Float64)
    return 10.0^((p_dBm) / 10.0)
end

@inline function mW_to_dBm(p_mW::Float64)
    return p_mW > 0 ? 10.0 * log10(p_mW) : -Inf
end

# ===== ポアソン点過程による端末配置 =====
function generate_poisson_positions(num_devices::Int, rng)
    positions = Vector{Tuple{Float64,Float64}}(undef, num_devices)
    
    for i in 1:num_devices
        theta = rand(rng) * 2π
        r = sqrt(rand(rng)) * area_size
        x = r * cos(theta)
        y = r * sin(theta)
        positions[i] = (x, y)
    end
    
    return positions
end

# ===== クロックドリフト管理 =====
mutable struct ClockDrift
    drift_ppm::Float64           # クロックドリフト（ppm）
    accumulated_drift::Float64   # 蓄積された時刻ずれ（秒）
    last_update_time::Float64    # 最後の更新時刻
end

function ClockDrift(drift_ppm::Float64)
    return ClockDrift(drift_ppm, 0.0, 0.0)
end

function update_clock_drift(clock::ClockDrift, current_time::Float64)
    # 前回の更新からの経過時間
    elapsed_time = current_time - clock.last_update_time
    
    # クロックドリフトによる時刻ずれの蓄積
    drift_seconds = clock.drift_ppm * 1e-6 * elapsed_time
    clock.accumulated_drift += drift_seconds
    clock.last_update_time = current_time
    
    return clock.accumulated_drift
end

function get_drifted_time(clock::ClockDrift, reference_time::Float64)
    # 基準時刻にクロックドリフトを適用
    return reference_time + clock.accumulated_drift
end

# クロックドリフト生成関数
function generate_clock_drift(rng, accuracy_level::Symbol=:standard)
    if accuracy_level == :high
        drift_ppm = rand(rng) * (CLOCK_DRIFT_HIGH_ACCURACY_MAX_PPM - CLOCK_DRIFT_HIGH_ACCURACY_MIN_PPM) + CLOCK_DRIFT_HIGH_ACCURACY_MIN_PPM
    else  # :standard
        drift_ppm = rand(rng) * (CLOCK_DRIFT_MAX_PPM - CLOCK_DRIFT_MIN_PPM) + CLOCK_DRIFT_MIN_PPM
    end
    return ClockDrift(drift_ppm)
end

# ===== DC制約管理 =====
mutable struct DCConstraint
    dc_limit::Float64          # DC制約（例: 0.01 = 1%）
    dc_remaining::Float64      # 残りDC
    last_reset_time::Float64   # 最後のリセット時刻
    window_duration::Float64   # DC計算ウィンドウ（秒）
end

function DCConstraint(dc_limit::Float64, window_duration::Float64=3600.0)
    return DCConstraint(dc_limit, dc_limit, 0.0, window_duration)
end

function update_dc_constraint(dc::DCConstraint, current_time::Float64, transmission_duration::Float64)
    # ウィンドウがリセットされるかチェック
    if current_time - dc.last_reset_time >= dc.window_duration
        dc.dc_remaining = dc.dc_limit
        dc.last_reset_time = current_time
    end
    
    # DC消費
    dc_consumption = transmission_duration / dc.window_duration
    dc.dc_remaining = max(0.0, dc.dc_remaining - dc_consumption)
    
    return dc.dc_remaining > 0
end

# ===== キャリアセンス状態管理 =====
mutable struct CarrierSenseState
    is_cs_active::Bool           # キャリアセンス実行中
    cs_detection::Bool           # キャリア検出フラグ
    cs_start_time::Float64       # キャリアセンス開始時刻
    cs_duration::Float64         # キャリアセンス時間
    backoff_until::Float64       # バックオフ終了時刻
    is_in_backoff::Bool          # バックオフ中フラグ
end

function CarrierSenseState(cs_duration::Float64=0.005)  # 5msのキャリアセンス
    return CarrierSenseState(false, false, 0.0, cs_duration, 0.0, false)
end

# ===== キャリアセンス関数 =====
function carrier_sense_check(device_id::Int, 
                           transmitting_devices::Vector{Int},
                           device_positions::Vector{Tuple{Float64,Float64}},
                           device_snrs::Vector{Float64},
                           shadowing_values::Vector{Float64},
                           carrier_sense_threshold::Float64)
    
    # 同じスロットで送信中の端末をチェック
    if length(transmitting_devices) <= 1
        return false  # キャリア未検出
    end
    
    # 受信電力の合計を計算
    total_power_mW = 0.0
    device_pos = device_positions[device_id]
    
    for other_id in transmitting_devices
        if other_id != device_id
            # 距離計算
            other_pos = device_positions[other_id]
            dist_km = distance(device_pos[1], device_pos[2], other_pos[1], other_pos[2])
            
            # 受信電力計算（シャドウイング考慮）
            received_pwr_dBm = received_power(Tx_dB, dist_km, shadowing_values[other_id])
            received_pwr_mW = dBm_to_mW(received_pwr_dBm)
            
            total_power_mW += received_pwr_mW
        end
    end
    
    # キャリアセンス閾値と比較
    if total_power_mW > 0
        total_power_dBm = mW_to_dBm(total_power_mW)
        return total_power_dBm >= carrier_sense_threshold
    else
        return false
    end
end

# ===== バックオフ処理 =====
function calculate_backoff_time(base_backoff::Float64, max_backoff::Float64, rng)
    # 指数バックオフ
    backoff_time = base_backoff * (2.0^rand(rng, 0:10))  # 最大1024倍
    return min(backoff_time, max_backoff)
end

# ===== SNR判定関数 =====
function snr_judge(device_snr::Float64, snr_threshold::Float64)
    return device_snr >= snr_threshold ? 0 : 1  # 0: 良好, 1: 不良
end

# ===== 衝突判定関数 =====
function collision_judge(transmitting_devices::Vector{Int})
    return length(transmitting_devices) > 1
end

# ===== クロックドリフトを考慮したスロット開始時刻計算 =====
function calculate_drifted_slot_start(device_clock::ClockDrift, slot::Int, slot_duration::Float64, fs::Float64)
    # 基準スロット開始時刻
    reference_slot_start = (slot - 1) * slot_duration
    
    # クロックドリフトを適用した時刻
    drifted_time = get_drifted_time(device_clock, reference_slot_start)
    
    # サンプルインデックスに変換
    return Int(round(drifted_time * fs)) + 1
end

# ===== クロックドリフト + DC制約 + キャリアセンス付きSlottedALOHA実装 =====
function slotted_aloha_dc_cs_clock_drift_per(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=10.0,
    slot_duration::Float64=0.1,
    dc_limit::Float64=0.01,  # 1%のDC制約
    dc_window::Float64=3600.0,  # 1時間のDC計算ウィンドウ
    cs_duration::Float64=0.005,  # 5msのキャリアセンス
    backoff_base::Float64=0.1,   # 基本バックオフ時間
    backoff_max::Float64=10.0,   # 最大バックオフ時間
    clock_accuracy::Symbol=:standard,  # :standard または :high
    fixed_snr::Union{Float64, Nothing}=nothing,
    use_shadowing::Bool=true,
    rng=Random.default_rng())

    M  = 2^sf
    Ts = M / bw
    fs = bw

    # ポアソン点過程で端末配置生成
    node_positions = generate_poisson_positions(num_devices, rng)
    
    # シャドウイング値生成
    shadowing_values = use_shadowing ? [shadowing_value(rng) for _ in 1:num_devices] : zeros(num_devices)
    
    # 各端末のSNR計算
    if fixed_snr !== nothing
        device_snrs = fill(fixed_snr, num_devices)
        shadowing_values = zeros(num_devices)
    else
        if use_shadowing
            device_snrs = [calculate_snr(distance(pos[1], pos[2], 0.0, 0.0), shadowing_values[i]) for (i, pos) in enumerate(node_positions)]
        else
            device_snrs = [calculate_snr(distance(pos[1], pos[2], 0.0, 0.0)) for pos in node_positions]
        end
    end
    
    # スロット数計算
    num_slots = Int(ceil(sim_time / slot_duration))
    
    # DC制約管理
    dc_constraints = [DCConstraint(dc_limit, dc_window) for _ in 1:num_devices]
    
    # キャリアセンス状態管理
    cs_states = [CarrierSenseState(cs_duration) for _ in 1:num_devices]
    
    # クロックドリフト管理
    clock_drifts = [generate_clock_drift(rng, clock_accuracy) for _ in 1:num_devices]
    
    # 各スロットでの送信端末記録
    slot_transmissions = Vector{Vector{Int}}()
    slot_payloads = Vector{Vector{Vector{Int}}}()
    slot_packet_starts = Vector{Vector{Int}}()
    
    for slot in 1:num_slots
        push!(slot_transmissions, Int[])
        push!(slot_payloads, Vector{Vector{Int}}[])
        push!(slot_packet_starts, Int[])
    end
    
    # 各端末の送信確率（DC制約を考慮）
    base_transmission_probability = 0.1
    
    # 送信スケジュール生成（クロックドリフト + DC制約 + キャリアセンスを考慮）
    for slot in 1:num_slots
        current_time = slot * slot_duration
        transmitting_devices = Int[]
        
        for d in 1:num_devices
            # クロックドリフト更新
            update_clock_drift(clock_drifts[d], current_time)
            
            # バックオフ中チェック
            if cs_states[d].is_in_backoff
                if current_time >= cs_states[d].backoff_until
                    cs_states[d].is_in_backoff = false
                else
                    continue  # バックオフ中は送信不可
                end
            end
            
            # DC制約チェック
            if dc_constraints[d].dc_remaining > 0
                # 送信確率をDC制約に応じて調整
                adjusted_probability = base_transmission_probability * 
                                     min(1.0, dc_constraints[d].dc_remaining / dc_limit)
                
                if rand(rng) < adjusted_probability
                    # キャリアセンス実行
                    cs_states[d].is_cs_active = true
                    cs_states[d].cs_start_time = current_time
                    
                    # キャリアセンス結果
                    cs_detected = carrier_sense_check(d, transmitting_devices, node_positions, 
                                                   device_snrs, shadowing_values, carrier_sense_threshold)
                    
                    if cs_detected
                        # キャリア検出 → バックオフ
                        cs_states[d].cs_detection = true
                        cs_states[d].is_cs_active = false
                        backoff_time = calculate_backoff_time(backoff_base, backoff_max, rng)
                        cs_states[d].backoff_until = current_time + backoff_time
                        cs_states[d].is_in_backoff = true
                    else
                        # キャリア未検出 → 送信実行
                        cs_states[d].cs_detection = false
                        cs_states[d].is_cs_active = false
                        push!(transmitting_devices, d)
                        
                        # パケット生成
                        payload_len = rand(rng, payload_range[1]:payload_range[2])
                        payload = [rand(rng, 0:M-1) for _ in 1:payload_len]
                        push!(slot_payloads[slot], copy(payload))
                        
                        # クロックドリフトを考慮したパケット開始サンプル計算
                        packet_start = calculate_drifted_slot_start(clock_drifts[d], slot, slot_duration, fs)
                        push!(slot_packet_starts[slot], packet_start)
                        
                        # DC制約更新（送信予定）
                        transmission_duration = payload_len * Ts
                        update_dc_constraint(dc_constraints[d], current_time, transmission_duration)
                    end
                end
            end
        end
        
        slot_transmissions[slot] = transmitting_devices
    end
    
    # 受信バッファ
    total_samples = Int(ceil(sim_time * fs))
    rx_signal = zeros(ComplexF64, total_samples)
    
    # シンボルキャッシュ
    sym_cache = Dict{Int, Vector{ComplexF64}}()
    
    # 送信実行（スロット単位）
    for slot in 1:num_slots
        transmitting_devices = slot_transmissions[slot]
        
        if !isempty(transmitting_devices)
            for (device_idx, d) in enumerate(transmitting_devices)
                payload = slot_payloads[slot][device_idx]
                packet_start = slot_packet_starts[slot][device_idx]
                
                for s in 1:length(payload)
                    m = payload[s]
                    if !haskey(sym_cache, m)
                        sym_cache[m] = lora_symbol(sf, bw, m)
                    end
                    wave = sym_cache[m]
                    
                    start_idx = packet_start + (s-1) * M
                    end_idx = start_idx + M - 1
                    
                    if end_idx <= total_samples && start_idx >= 1
                        for i in 1:M
                            idx = start_idx + i - 1
                            if idx >= 1 && idx <= length(rx_signal) && i <= length(wave)
                                rx_signal[idx] += wave[i]
                            end
                        end
                    end
                end
            end
        end
    end
    
    # 各端末のSNRに基づくAWGN付加
    for slot in 1:num_slots
        transmitting_devices = slot_transmissions[slot]
        
        if !isempty(transmitting_devices)
            for (device_idx, d) in enumerate(transmitting_devices)
                snr = device_snrs[d]
                payload = slot_payloads[slot][device_idx]
                packet_start = slot_packet_starts[slot][device_idx]
                
                for s in 1:length(payload)
                    start_idx = packet_start + (s-1) * M
                    end_idx = start_idx + M - 1
                    
                    if end_idx <= total_samples && start_idx >= 1
                        signal_segment = @view rx_signal[start_idx:end_idx]
                        add_awgn!(signal_segment, snr)
                    end
                end
            end
        end
    end
    
    # パケット分類とPER計算（復調結果ベース）
    total_packets = 0
    success_packets = 0
    collision_packets = 0
    demod_fail_packets = 0  # 復調失敗パケット数
    cs_blocked_packets = 0  # キャリアセンスでブロックされたパケット数
    backoff_packets = 0     # バックオフでブロックされたパケット数
    clock_drift_packets = 0 # クロックドリフトによる送信タイミングずれパケット数
    
    for slot in 1:num_slots
        transmitting_devices = slot_transmissions[slot]
        
        if !isempty(transmitting_devices)
            has_collision = collision_judge(transmitting_devices)
            
            for (device_idx, d) in enumerate(transmitting_devices)
                total_packets += 1
                
                if has_collision  # 衝突
                    collision_packets += 1
                else  # 復調チェック（すべてのSNRで実行）
                    payload = slot_payloads[slot][device_idx]
                    packet_start = slot_packet_starts[slot][device_idx]
                    packet_success = true
                    
                    for s in 1:length(payload)
                        start_idx = packet_start + (s-1) * M
                        end_idx = start_idx + M - 1
                        
                        if end_idx <= total_samples && start_idx >= 1 && end_idx <= length(rx_signal)
                            y = @view rx_signal[start_idx:end_idx]
                            m_hat = demod_lora(sf, bw, y)
                            
                            if m_hat != payload[s]
                                packet_success = false
                                break
                            end
                        else
                            packet_success = false
                            break
                        end
                    end
                    
                    if packet_success
                        success_packets += 1
                    else
                        demod_fail_packets += 1  # 復調失敗
                    end
                end
            end
        end
    end
    
    # キャリアセンスとバックオフでブロックされたパケット数を計算
    for d in 1:num_devices
        for slot in 1:num_slots
            current_time = slot * slot_duration
            if cs_states[d].is_in_backoff && current_time < cs_states[d].backoff_until
                backoff_packets += 1
            elseif cs_states[d].cs_detection
                cs_blocked_packets += 1
            end
        end
    end
    
    # クロックドリフトによる送信タイミングずれを検出
    for d in 1:num_devices
        for slot in 1:num_slots
            current_time = slot * slot_duration
            update_clock_drift(clock_drifts[d], current_time)
            
            # クロックドリフトがスロット時間の一定割合を超えた場合
            if abs(clock_drifts[d].accumulated_drift) > slot_duration * 0.1  # 10%以上ずれ
                clock_drift_packets += 1
            end
        end
    end
    
    # PER計算（復調結果ベース）
    per = total_packets > 0 ? (collision_packets + demod_fail_packets) / total_packets : 0.0
    
    # 統計情報
    stats = Dict(
        :total_packets => total_packets,
        :success_packets => success_packets,
        :collision_packets => collision_packets,
        :demod_fail_packets => demod_fail_packets,
        :cs_blocked_packets => cs_blocked_packets,
        :backoff_packets => backoff_packets,
        :clock_drift_packets => clock_drift_packets,
        :per => per,
        :collision_rate => total_packets > 0 ? collision_packets / total_packets : 0.0,
        :demod_fail_rate => total_packets > 0 ? demod_fail_packets / total_packets : 0.0,
        :cs_blocked_rate => total_packets > 0 ? cs_blocked_packets / total_packets : 0.0,
        :backoff_rate => total_packets > 0 ? backoff_packets / total_packets : 0.0,
        :clock_drift_rate => total_packets > 0 ? clock_drift_packets / total_packets : 0.0
    )
    
    return per, node_positions, device_snrs, stats, clock_drifts
end

# ===== クロックドリフト + DC制約 + キャリアセンス付きSlottedALOHA SNRスイープ機能 =====
function run_slotted_aloha_dc_cs_clock_drift_snr_sweep(sf::Int, bw::Float64, num_devices::Int,
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=10.0,
    slot_duration::Float64=0.1,
    dc_limit::Float64=0.01,
    dc_window::Float64=3600.0,
    cs_duration::Float64=0.005,
    backoff_base::Float64=0.1,
    backoff_max::Float64=10.0,
    clock_accuracy::Symbol=:standard,
    iter::Int=100,
    seed::Int=1234,
    save_path::String="results_LoRa_simple_model/slotted_aloha_dc_cs_clock_drift_snr_sweep_per.csv",
    use_shadowing::Bool=true)

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    per_vals = zeros(Float64, length(snrs))
    
    println("クロックドリフト + DC制約 + キャリアセンス付きSlottedALOHA SNRスイープ実行中...")
    println("SNR範囲: $(snr_min) ~ $(snr_max) dB, ステップ: $(snr_step) dB")
    println("反復回数: $iter")
    println("シミュレーション時間: $sim_time 秒")
    println("スロット時間: $slot_duration 秒")
    println("DC制約: $(dc_limit*100)%")
    println("DC計算ウィンドウ: $dc_window 秒")
    println("キャリアセンス時間: $(cs_duration*1000) ms")
    println("バックオフ時間: $(backoff_base) ~ $(backoff_max) 秒")
    println("クロック精度: $clock_accuracy")
    
    for (i, snr) in enumerate(snrs)
        acc_per = 0.0
        
        for it in 1:iter
            local_rng = MersenneTwister(seed + it + i * 1000)
            
            per, positions, snrs_actual, stats, clock_drifts = slotted_aloha_dc_cs_clock_drift_per(sf, bw, num_devices;
                                               payload_range=payload_range,
                                               sim_time=sim_time,
                                               slot_duration=slot_duration,
                                               dc_limit=dc_limit,
                                               dc_window=dc_window,
                                               cs_duration=cs_duration,
                                               backoff_base=backoff_base,
                                               backoff_max=backoff_max,
                                               clock_accuracy=clock_accuracy,
                                               fixed_snr=snr,
                                               use_shadowing=use_shadowing,
                                               rng=local_rng)
            
            acc_per += per
        end
        
        per_vals[i] = acc_per / iter
        
        @printf("SNR = %5.1f dB, PER = %.6f\n", snr, per_vals[i])
    end

    # CSV保存
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    
    open(save_path, "w") do io
        println(io, "SNR_dB,PER")
        for i in 1:length(snrs)
            println(io, "$(snrs[i]),$(per_vals[i])")
        end
    end
    
    @info "クロックドリフト + DC制約 + キャリアセンス付きSlottedALOHA SNRスイープ結果をCSVに保存しました: $save_path"
    
    return snrs, per_vals
end

# ===== メイン実行部分 =====
println("クロックドリフト + DC制約 + キャリアセンス付きLoRa SlottedALOHA シミュレーション")
println("設定:")
println("  SF: $(SF)")
println("  エリアサイズ: $(area_size) km")
println("  送信電力: $(Tx_dB) dBm")
println("  搬送波周波数: $(f_c) MHz")
println("  パスロス係数: α=$(α), β=$(β), γ=$(γ)")
println("  シャドウイング標準偏差: $(shadowing_std) dB")
println("  SNR閾値: $(SNR_threshold) dB")
println("  キャリアセンス閾値: $(carrier_sense_threshold) dBm")
println("  クロックドリフト範囲: $(CLOCK_DRIFT_MIN_PPM) ~ $(CLOCK_DRIFT_MAX_PPM) ppm")

# パラメータ設定
sf = SF
bw = 125e3
num_devices = 10

# クロックドリフト + DC制約 + キャリアセンス付きSlottedALOHA SNRスイープ実行
println("\n=== クロックドリフト + DC制約 + キャリアセンス付きSlottedALOHA SNRスイープ実行 ===")
snr_min, snr_max, snr_step = -20.0, 0.0, 0.5
iter_sweep = 1000
dc_limit = 0.01  # 1%のDC制約
cs_duration = 0.005  # 5msのキャリアセンス
backoff_base = 0.1   # 基本バックオフ時間
backoff_max = 10.0   # 最大バックオフ時間
clock_accuracy = :standard  # 標準精度のクロック

snrs_clock_drift, per_vals_clock_drift = run_slotted_aloha_dc_cs_clock_drift_snr_sweep(sf, bw, num_devices,
                               snr_min, snr_max, snr_step;
                               payload_range=(16,32),
                               sim_time=10.0,
                               slot_duration=0.1,
                               dc_limit=dc_limit,
                               dc_window=3600.0,
                               cs_duration=cs_duration,
                               backoff_base=backoff_base,
                               backoff_max=backoff_max,
                               clock_accuracy=clock_accuracy,
                               iter=iter_sweep,
                               use_shadowing=true,
                               save_path="LoRa_PER/LoRa_PL_DC_CS/results_SlottedALOHA/slotted_aloha_dc_cs_clock_drift_snr_sweep_per_sf$(sf)_dev$(num_devices)_iter$(iter_sweep).csv")

# 結果の表示
println("\n=== クロックドリフト + DC制約 + キャリアセンス付きSlottedALOHA SNRスイープ結果 ===")
println("SNR範囲: $(snr_min) ~ $(snr_max) dB")
println("DC制約: $(dc_limit*100)%")
println("キャリアセンス時間: $(cs_duration*1000) ms")
println("バックオフ時間: $(backoff_base) ~ $(backoff_max) 秒")
println("クロック精度: $clock_accuracy")
println("最小PER: $(minimum(per_vals_clock_drift))")
println("最大PER: $(maximum(per_vals_clock_drift))")

println("\nシミュレーション完了！")
