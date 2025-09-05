using Random, FFTW, Statistics, Printf, DelimitedFiles

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


# 逆チャープ（復調用キャッシュ）
const CHIRP_COMP = Dict{Tuple{Int,Float64}, Vector{ComplexF64}}()
function get_chirp_comp(sf::Int, bw::Float64)
    get!(CHIRP_COMP, (sf,bw)) do
        M  = 2^sf
        tc = 1.0 / bw
        ts = M * tc
        exp.(-1im * 2π * (bw/(2*ts)) * ((0:M-1).^2) * tc^2)
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
function add_awgn!(y::Vector{ComplexF64}, snr_dB::Float64)
    var = 10^(-snr_dB/10)
    σ = sqrt(var/2)
    @inbounds for i in eachindex(y)
        y[i] += ComplexF64(randn()*σ, randn()*σ)
    end
    return y
end

# ===== CSMA + DC + PER =====
function csma_dc_per(sf::Int, bw::Float64, num_devices::Int;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    snr_dB::Float64=0.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    rng=Random.default_rng())

    M  = 2^sf
    Ts = M / bw     # シンボル長
    fs = bw         # サンプル周波数 ≈ 帯域幅

    # 送信スケジュール
    tx_starts = rand(rng, num_devices) .* sim_time
    payload_lens = rand(rng, payload_range[1]:payload_range[2], num_devices)
    payloads = [rand(rng, 0:M-1, payload_lens[d]) for d in 1:num_devices]

    # 受信バッファ
    total_samples = Int(ceil(sim_time*fs)) + maximum(payload_lens)*M
    rx_signal = zeros(ComplexF64, total_samples)
    channel_busy = falses(total_samples)

    # シンボルキャッシュ
    sym_cache = Dict{Int, Vector{ComplexF64}}()

    # DC制約
    next_time = zeros(Float64, num_devices)

    # 実際に送信
    tx_indices = [Int[] for _ in 1:num_devices]
    for d in sortperm(tx_starts)
        t = max(tx_starts[d], next_time[d])
        L = payload_lens[d]
        first_start = 0
        for s in 1:L
            m = payloads[d][s]
            if !haskey(sym_cache, m)
                sym_cache[m] = lora_symbol(sf, bw, m)
            end
            wave = sym_cache[m]
            n0 = Int(floor(t*fs)) + 1
            n1 = n0 + M - 1

            # CSMA: busyなら待機＋バックオフ
            while n1 <= total_samples && any(channel_busy[n0:n1])
                t += rand(rng)*backoff_max
                n0 = Int(floor(t*fs)) + 1
                n1 = n0 + M - 1
            end

            if n1 <= total_samples
                channel_busy[n0:n1] .= true
                @inbounds for i in 1:M
                    rx_signal[n0+i-1] += wave[i]
                end
                if s==1; first_start=n0; end
                push!(tx_indices[d], n0)
            end
            t += Ts
        end

        # DC制約
        T_on = L * Ts
        tend = (first_start==0) ? t : ((first_start-1)/fs + T_on)
        if dc>0
            T_off = T_on * (1-dc)/dc
            next_time[d] = tend + T_off
        end
    end

    # AWGN
    add_awgn!(rx_signal, snr_dB)

    # 復調＋PER判定
    errors = falses(num_devices)
    for d in 1:num_devices
        for (s,start) in enumerate(tx_indices[d])
            y = @view rx_signal[start:start+M-1]
            m_hat = demod_lora(sf, bw, y)
            if m_hat != payloads[d][s]
                errors[d] = true
                break
            end
        end
    end

    return mean(Float64.(errors))
end

# ===== SNRスイープ (PER vs SNR) =====
function run_per_sweep(sf::Int, bw::Float64,
    num_devices::Int,
    snr_min::Float64, snr_max::Float64, snr_step::Float64;
    payload_range::Tuple{Int,Int}=(16,32),
    sim_time::Float64=1.0,
    backoff_max::Float64=0.01,
    dc::Float64=0.01,
    iter::Int=100,
    seed::Int=1234,
    save_path::String="LoRa_PER_sweep.csv")

    rng = MersenneTwister(seed)
    snrs = collect(snr_min:snr_step:snr_max)
    per_vals = zeros(Float64, length(snrs))

    for (i, snr) in enumerate(snrs)
        acc = 0.0
        for it in 1:iter
            acc += csma_dc_per(sf, bw, num_devices;
                               payload_range=payload_range,
                               sim_time=sim_time,
                               snr_dB=snr,
                               backoff_max=backoff_max,
                               dc=dc,
                               rng=rng)
        end
        per_vals[i] = acc / iter
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
    @info "結果をCSVに保存しました: $save_path"

    return snrs, per_vals
end

sf = 7
bw = 125e3
num_devices = 50
snr_min, snr_max, snr_step = -10.0, 0.0, 1.0
iter = 50   # 反復回数（精度を上げたいときは増やす）

snrs, per_vals = run_per_sweep(sf, bw, num_devices,
                               snr_min, snr_max, snr_step;
                               payload_range=(16,32),
                               sim_time=1.0,
                               backoff_max=0.01,
                               dc=0.01,
                               iter=iter,
                               save_path="LoRa_SER/CSMA_DC/LoRa_PER_sf$(sf)_dev$(num_devices).csv")
