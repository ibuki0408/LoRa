using FFTW
using Statistics
using LinearAlgebra
using Plots

# -----------------------------
# 複素正規分布ノイズ生成
# -----------------------------
function nom_dist(mean::Float64, var::Float64)
    x = rand()
    y = rand()
    Re = sqrt(-2.0*var*log(x))*cos(2pi*y) + mean
    Im = sqrt(-2.0*var*log(x))*sin(2pi*y) + mean
    return Re + im*Im
end

function add_noise(mean::Float64,var::Float64)
    return nom_dist(mean,var/2.0)
end

# -----------------------------
# DFT
# -----------------------------
function dft(y::Vector{ComplexF64})
    M = length(y)
    Y = zeros(ComplexF64, M)
    for n in 0:M-1
        for k in 0:M-1
            Y[n+1] += y[k+1]*exp(-im*2pi*n*k/M)
        end
    end
    return Y
end

# -----------------------------
# LoRa 復調（単端末）
# -----------------------------
function demod_lora(sf::Int, bw::Float64, tc::Float64, ts::Float64, x::Vector{ComplexF64})
    M = 2^sf
    # LoRa de-chirp
    d = [x[k] * exp(-im*2pi*(bw/(2*ts))*(k-1)^2*tc^2) for k in 1:M]
    Y = dft(d)
    Y_abs = abs.(Y)/M
    m_star = argmax(Y_abs) - 1
    return m_star
end

# -----------------------------
# SER 計算（単端末）
# -----------------------------
function ser_calc_single(snr::Float64, sf::Int, bw::Float64, iter::Int)
    M = 2^sf
    tc = 1/bw
    ts = tc*M
    var = 1.0/10^(snr/10)
    error = 0

    for i in 1:iter
        m = rand(0:M-1)
        # LoRa チャープ信号生成
        x = [exp(im*2pi*(bw/(2*ts))*((k-1)+m)^2*tc^2) for k in 1:M]
        # ノイズ付加
        for k in 1:M
            x[k] += add_noise(0.0, var)
        end
        m_star = demod_lora(sf, bw, tc, ts, x)
        if m != m_star
            error += 1
        end
    end

    return error/iter
end

# -----------------------------
# メインループ（SNR sweep）
# -----------------------------
function main_lora(sf::Int, bw::Float64, iter::Int)
    minSNR = -20.0
    maxSNR = 10.0
    step = 2.0
    snrs = minSNR:step:maxSNR
    ser_list = Float64[]

    for snr in snrs
        ser = ser_calc_single(snr, sf, bw, iter)
        println("SNR = $snr dB, SER = $ser")
        push!(ser_list, ser)
    end

    plot(snrs, ser_list,
        xlabel="SNR [dB]",
        ylabel="Symbol Error Rate",
        title="Single-device LoRa SER",
        yscale=:log10,
        legend=false,
        grid=true)
end

# -----------------------------
# 実行例
# -----------------------------
main_lora(7, 125_000.0, 100000)
