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

# ===== SNR判定関数（先輩のコードを参考） =====
function snr_judge(device_snr::Float64, snr_threshold::Float64)
    # SNR閾値と比較して判定
    return device_snr >= snr_threshold ? 0 : 1  # 0: 良好, 1: 不良
end

# ===== 衝突判定関数（先輩のコードを参考） =====
function collision_judge(transmitting_devices::Vector{Int})
    # 複数端末が同じスロットで送信した場合、衝突
    return length(transmitting_devices) > 1
end

# ===== 修正版SlottedALOHA実装 =====
function slotted_aloha_per_fixed(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=10.0,
    slot_duration::Float64=0.1,
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
    
    # 各スロットでの送信端末記録
    slot_transmissions = Vector{Vector{Int}}()  # 各スロットで送信した端末
    slot_payloads = Vector{Vector{Vector{Int}}}()  # 各スロットの送信パケット
    slot_packet_starts = Vector{Vector{Int}}()  # 各スロットのパケット開始サンプル
    
    for slot in 1:num_slots
        push!(slot_transmissions, Int[])
        push!(slot_payloads, Vector{Vector{Int}}[])
        push!(slot_packet_starts, Int[])
    end
    
    # 各端末の送信確率
    transmission_probability = 0.1
    
    # 送信スケジュール生成（スロット単位）
    for slot in 1:num_slots
        transmitting_devices = Int[]
        
        for d in 1:num_devices
            if rand(rng) < transmission_probability
                push!(transmitting_devices, d)
                
                # パケット生成
                payload_len = rand(rng, payload_range[1]:payload_range[2])
                payload = [rand(rng, 0:M-1) for _ in 1:payload_len]
                push!(slot_payloads[slot], copy(payload))
                
                # パケット開始サンプル計算
                slot_start_sample = Int(round((slot - 1) * slot_duration * fs)) + 1
                push!(slot_packet_starts[slot], slot_start_sample)
            end
        end
        
        # 送信端末を記録
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
            # 各送信端末のパケットを送信
            for (device_idx, d) in enumerate(transmitting_devices)
                payload = slot_payloads[slot][device_idx]
                packet_start = slot_packet_starts[slot][device_idx]
                
                # パケット内の各シンボルを送信
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
    
    # パケット分類とPER計算（先輩のコードを参考）
    total_packets = 0
    success_packets = 0
    collision_packets = 0
    snr_fail_packets = 0
    
    for slot in 1:num_slots
        transmitting_devices = slot_transmissions[slot]
        
        if !isempty(transmitting_devices)
            # 衝突判定
            has_collision = collision_judge(transmitting_devices)
            
            for (device_idx, d) in enumerate(transmitting_devices)
                total_packets += 1
                
                # SNR判定
                snr_result = snr_judge(device_snrs[d], SNR_threshold)
                
                if snr_result == 1  # SNR不良
                    snr_fail_packets += 1
                elseif has_collision  # 衝突
                    collision_packets += 1
                else  # 成功の可能性（復調チェック）
                    payload = slot_payloads[slot][device_idx]
                    packet_start = slot_packet_starts[slot][device_idx]
                    packet_success = true
                    
                    # パケット内の各シンボルを復調
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
                        snr_fail_packets += 1  # 復調失敗はSNR不良として扱う
                    end
                end
            end
        end
    end
    
    # PER計算
    per = total_packets > 0 ? (collision_packets + snr_fail_packets) / total_packets : 0.0
    
    # 統計情報
    stats = Dict(
        :total_packets => total_packets,
        :success_packets => success_packets,
        :collision_packets => collision_packets,
        :snr_fail_packets => snr_fail_packets,
        :per => per,
        :collision_rate => total_packets > 0 ? collision_packets / total_packets : 0.0,
        :snr_fail_rate => total_packets > 0 ? snr_fail_packets / total_packets : 0.0
    )
    
    return per, node_positions, device_snrs, stats
end

# ===== 修正版SlottedALOHA SNRスイープ機能 =====
function run_slotted_aloha_snr_sweep_fixed(sf::Int, bw::Float64, num_devices::Int,
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=10.0,
    slot_duration::Float64=0.1,
    iter::Int=100,
    seed::Int=1234,
    save_path::String="results_LoRa_simple_model/slotted_aloha_fixed_snr_sweep_per.csv",
    use_shadowing::Bool=true)

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    per_vals = zeros(Float64, length(snrs))
    collision_rates = zeros(Float64, length(snrs))
    snr_fail_rates = zeros(Float64, length(snrs))
    
    println("修正版SlottedALOHA SNRスイープ実行中...")
    println("SNR範囲: $(snr_min) ~ $(snr_max) dB, ステップ: $(snr_step) dB")
    println("反復回数: $iter")
    println("シミュレーション時間: $sim_time 秒")
    println("スロット時間: $slot_duration 秒")
    
    for (i, snr) in enumerate(snrs)
        acc_per = 0.0
        acc_collision_rate = 0.0
        acc_snr_fail_rate = 0.0
        
        for it in 1:iter
            local_rng = MersenneTwister(seed + it + i * 1000)
            
            per, positions, snrs_actual, stats = slotted_aloha_per_fixed(sf, bw, num_devices;
                                               payload_range=payload_range,
                                               sim_time=sim_time,
                                               slot_duration=slot_duration,
                                               fixed_snr=snr,
                                               use_shadowing=use_shadowing,
                                               rng=local_rng)
            
            acc_per += per
            acc_collision_rate += stats[:collision_rate]
            acc_snr_fail_rate += stats[:snr_fail_rate]
        end
        
        per_vals[i] = acc_per / iter
        collision_rates[i] = acc_collision_rate / iter
        snr_fail_rates[i] = acc_snr_fail_rate / iter
        
        @printf("SNR = %5.1f dB, PER = %.6f, Collision = %.6f, SNR_Fail = %.6f\n", 
                snr, per_vals[i], collision_rates[i], snr_fail_rates[i])
    end

    # CSV保存
    if !isdir(dirname(save_path))
        mkpath(dirname(save_path))
    end
    
    open(save_path, "w") do io
        println(io, "SNR_dB,PER,Collision_Rate,SNR_Fail_Rate")
        for i in 1:length(snrs)
            println(io, "$(snrs[i]),$(per_vals[i]),$(collision_rates[i]),$(snr_fail_rates[i])")
        end
    end
    
    @info "修正版SlottedALOHA SNRスイープ結果をCSVに保存しました: $save_path"
    
    return snrs, per_vals, collision_rates, snr_fail_rates
end

# ===== メイン実行部分 =====
println("修正版LoRa SlottedALOHA シミュレーション")
println("設定:")
println("  SF: $(SF)")
println("  エリアサイズ: $(area_size) km")
println("  送信電力: $(Tx_dB) dBm")
println("  搬送波周波数: $(f_c) MHz")
println("  パスロス係数: α=$(α), β=$(β), γ=$(γ)")
println("  シャドウイング標準偏差: $(shadowing_std) dB")
println("  SNR閾値: $(SNR_threshold) dB")

# パラメータ設定
sf = SF
bw = 125e3
num_devices = 2

# 修正版SlottedALOHA SNRスイープ実行
println("\n=== 修正版SlottedALOHA SNRスイープ実行 ===")
snr_min, snr_max, snr_step = -20.0, 0.0, 0.5
iter_sweep = 1000

snrs_fixed, per_vals_fixed, collision_rates, snr_fail_rates = run_slotted_aloha_snr_sweep_fixed(sf, bw, num_devices,
                               snr_min, snr_max, snr_step;
                               payload_range=(16,32),
                               sim_time=10.0,
                               slot_duration=0.1,
                               iter=iter_sweep,
                               use_shadowing=true,
                               save_path="LoRa_PER/LoRa_PL_DC_CS/results_SlottedALOHA/slotted_aloha_fixed_snr_sweep_per_sf$(sf)_dev$(num_devices)_iter$(iter_sweep).csv")

# 結果の表示
println("\n=== 修正版SlottedALOHA SNRスイープ結果 ===")
println("SNR範囲: $(snr_min) ~ $(snr_max) dB")
println("最小PER: $(minimum(per_vals_fixed))")
println("最大PER: $(maximum(per_vals_fixed))")
println("平均衝突率: $(mean(collision_rates))")
println("平均SNR失敗率: $(mean(snr_fail_rates))")

println("\nシミュレーション完了！")
