using FFTW
using Statistics
using LinearAlgebra
using Random
using Printf
using DelimitedFiles   # ← CSV 書き出し用
# using Plots         # ← プロット不要なら消す

# ===== LoRa 基本 =====
function lora_symbol(sf::Int, bw::Float64, m::Int)
    M  = 2^sf
    tc = 1.0/bw
    ts = tc*M
    x = Vector{ComplexF64}(undef, M)
    @inbounds for k in 1:M
        x[k] = exp(1im * 2π * (bw/(2*ts)) * ((k-1) + m)^2 * tc^2)
    end
    return x
end

function demod_lora(sf::Int, bw::Float64, x::AbstractVector{ComplexF64})
    M  = 2^sf
    tc = 1.0/bw
    ts = tc*M
    d = Vector{ComplexF64}(undef, M)
    @inbounds for k in 1:M
        d[k] = x[k] * exp(-1im * 2π * (bw/(2*ts)) * (k-1)^2 * tc^2)
    end
    Y = fft(d)
    return argmax(abs.(Y)) - 1
end

@inline function add_awgn!(y::Vector{ComplexF64}, snr_dB::Float64)
    var = 10^(-snr_dB/10)
    noise = sqrt(var/2)
    @inbounds @simd for i in eachindex(y)
        y[i] += ComplexF64(randn()*noise, randn()*noise)
    end
    return y
end

function ser_multi_device_csma_avg(
    snr_dB::Float64, sf::Int, bw::Float64, iter::Int;
    num_devices::Int=4, frame_len::Int=32, rng=Random.default_rng()
)
    M  = 2^sf
    errors = zeros(Int, num_devices)  # 各端末の誤り数
    rx_frame = zeros(ComplexF64, frame_len*M)

    for t in 1:iter
        fill!(rx_frame, 0.0 + 0.0im)
        slot_taken = falses(frame_len)
        m_vec = zeros(Int, num_devices)
        slot_vec = zeros(Int, num_devices)

        # 端末ごとに空きスロットを探して送信
        @inbounds for d in 1:num_devices
            m_vec[d] = rand(rng, 0:M-1)
            available_slots = findall(!, slot_taken)
            if isempty(available_slots)
                s = rand(rng, 0:frame_len-1)
            else
                s = rand(rng, available_slots) - 1
                slot_taken[s+1] = true
            end
            slot_vec[d] = s

            sym  = lora_symbol(sf, bw, m_vec[d])
            i1   = s*M + 1
            i2   = i1 + M - 1
            @views rx_frame[i1:i2] .+= sym
        end

        # AWGN追加
        add_awgn!(rx_frame, snr_dB)

        # 各端末を復調して誤りをカウント
        @inbounds for d in 1:num_devices
            s = slot_vec[d]
            i1 = s*M + 1
            i2 = i1 + M - 1
            @views ywin = rx_frame[i1:i2]
            m_hat = demod_lora(sf, bw, ywin)
            if m_hat != m_vec[d]
                errors[d] += 1
            end
        end
    end

    # 各端末の誤り率の平均を返す
    return mean(errors ./ iter)
end


# ===== SNRスイープ＆CSV出力 =====
function run_sweep_and_save_csv(
    filename::Union{Nothing,String}=nothing;
    save_dir::Union{Nothing,String}=nothing,
    sf::Int=7, bw::Float64=125_000.0,
    snr_min::Float64=-20.0, snr_max::Float64=10.0, snr_step::Float64=2.0,
    iter::Int=2000, device_list::Vector{Int} = [1,2,4,8], frame_len::Int=32
)
    snrs = collect(snr_min:snr_step:snr_max)
    ser_mat = zeros(Float64, length(snrs), length(device_list)+1)
    ser_mat[:,1] .= snrs

    for (j, nd) in enumerate(device_list)
        @info "Simulating num_devices = $nd"
        for (i, snr) in enumerate(snrs)
            ser = ser_multi_device_csma_avg(snr, sf, bw, iter; num_devices=nd, frame_len=frame_len)
            ser_mat[i, j+1] = ser
            @printf "SNR=%5.1f dB  SER=%.6f\n" snr ser
        end
    end

    devices_str = join(device_list, "-")  # 例: [2,4,8,16] → "2-4-8-16"

    # 保存先ファイル名を決定
    if filename === nothing
        base = @sprintf("CS_multi_ser_sf%d_devices%s_frame%d_snr%.1f-%.1f_iter%d.csv",
                        sf, devices_str, frame_len, snr_min, snr_max, iter)
        if save_dir !== nothing
            filename = joinpath(save_dir, base)
        else
            filename = base
        end
    end

    # 保存先フォルダがなければ作成
    mkpath(dirname(filename))

    # 書き込み
    header = ["SNR_dB"; ["devices=$(nd)" for nd in device_list]...]
    open(filename, "w") do io
        println(io, join(header, ","))
        writedlm(io, ser_mat, ',')
    end

    @info "結果をCSVに保存しました: $filename"
    return snrs, ser_mat
end

snrs, ser = run_sweep_and_save_csv(nothing; 
    save_dir="/Users/kimparaibuki/Cursor/LoRa_SER/result",
    sf=7, bw=125_000.0, snr_min=-20.0, snr_max=0.0, snr_step=0.5,
    iter=1000, device_list=[2,4,8,16], frame_len=128)
