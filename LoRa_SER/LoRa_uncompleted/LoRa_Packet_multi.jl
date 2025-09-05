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

# ===== パケット単位（複数シンボル連続送信）評価 =====
"""
    per_ser_multi_device_time(
        snr_dB, sf, bw, iter;
        num_devices=4,
        frame_len_syms=4096,
        pkt_syms=256,
        wraparound=false,
        rng=Random.default_rng()
    ) -> (avg_SER, avg_PER)

時間スロット＝シンボル時間。各端末は開始スロットを選び、連続 pkt_syms 個のシンボルを送信。
受信はシンボルごとに復調。端末ごとの SER と PER を計測し平均を返す。

- frame_len_syms : 1フレーム内の総シンボル数（時間スロット数）
- pkt_syms       : パケットの総シンボル数（前置きやヘッダを含めた等価長として扱ってOK）
- wraparound     : パケットが末尾を越えるときに先頭へ回す（true）か、送信スキップ（false）
"""
function per_ser_multi_device_time(
    snr_dB::Float64, sf::Int, bw::Float64, iter::Int;
    num_devices::Int=4, frame_len_syms::Int=4096, pkt_syms::Int=256,
    wraparound::Bool=false, rng=Random.default_rng()
)
    M  = 2^sf                 # 1 LoRaシンボルのサンプル数
    # タイムライン（フレーム全体）: 複素サンプルで frame_len_syms * M
    rx_timeline = zeros(ComplexF64, frame_len_syms * M)

    # 各端末の統計
    total_symbol_tx = zeros(Int, num_devices)
    total_symbol_err = zeros(Int, num_devices)
    total_pkt_tx = zeros(Int, num_devices)
    total_pkt_err = zeros(Int, num_devices)

    for _ in 1:iter
        fill!(rx_timeline, 0.0 + 0.0im)

        # --- 送信予定（開始スロット & シンボル列）を決める ---
        # 各デバイスの開始スロット（シンボルインデックス）
        start_slots = zeros(Int, num_devices)
        # 各デバイスの送るシンボル系列（長さ pkt_syms）
        msg_syms = [zeros(Int, pkt_syms) for _ in 1:num_devices]

        for d in 1:num_devices
            # ランダムな開始スロット
            s0 = rand(rng, 0:frame_len_syms-1)
            # パケットが末尾を越える場合の取り扱い
            if !wraparound && (s0 + pkt_syms > frame_len_syms)
                # 送信スキップ（このイテレーションでは送らない）
                start_slots[d] = -1
                continue
            end
            start_slots[d] = s0
            # ランダムな M-ary シンボル系列
            for n in 1:pkt_syms
                msg_syms[d][n] = rand(rng, 0:M-1)
            end
        end

        # --- タイムラインへ合成（パケット＝連続 pkt_syms シンボル） ---
        for d in 1:num_devices
            s0 = start_slots[d]
            if s0 < 0
                continue
            end
            for n in 0:pkt_syms-1
                m = msg_syms[d][n+1]
                sym = lora_symbol(sf, bw, m)
                # シンボルが配置されるタイムライン上の区間
                if wraparound
                    # ラップアラウンドあり
                    slot = (s0 + n) % frame_len_syms
                    i1 = slot*M + 1
                    i2 = i1 + M - 1
                    @views rx_timeline[i1:i2] .+= sym
                else
                    # ラップなし（範囲外は送らない）
                    slot = s0 + n
                    if slot >= frame_len_syms
                        break
                    end
                    i1 = slot*M + 1
                    i2 = i1 + M - 1
                    @views rx_timeline[i1:i2] .+= sym
                end
            end
        end

        # --- AWGNを付加 ---
        add_awgn!(rx_timeline, snr_dB)

        # --- 受信・復調（シンボルごと）＆誤りカウント ---
        for d in 1:num_devices
            s0 = start_slots[d]
            if s0 < 0
                continue
            end
            pkt_err = false
            sent_syms = 0
            err_syms = 0

            for n in 0:pkt_syms-1
                # 送っていない領域（wraparound=false & オーバーラン）はスキップ
                slot = wraparound ? (s0 + n) % frame_len_syms : (s0 + n)
                if (!wraparound) && (slot >= frame_len_syms)
                    break
                end

                i1 = slot*M + 1
                i2 = i1 + M - 1
                @views ywin = rx_timeline[i1:i2]
                m_hat = demod_lora(sf, bw, ywin)
                m_true = msg_syms[d][n+1]
                sent_syms += 1
                if m_hat != m_true
                    err_syms += 1
                    pkt_err = true
                end
            end

            total_symbol_tx[d] += sent_syms
            total_symbol_err[d] += err_syms
            # 「送信成立したパケット数」= 実際に1シンボル以上送れたときのみカウント
            if sent_syms > 0
                total_pkt_tx[d] += 1
                total_pkt_err[d] += pkt_err ? 1 : 0
            end
        end
    end

    # 端末ごとの平均 SER / PER を取り、さらにデバイス平均を返す
    ser_each = zeros(Float64, num_devices)
    per_each = zeros(Float64, num_devices)
    for d in 1:num_devices
        ser_each[d] = total_symbol_tx[d] == 0 ? NaN :
                      total_symbol_err[d] / total_symbol_tx[d]
        per_each[d] = total_pkt_tx[d] == 0 ? NaN :
                      total_pkt_err[d] / total_pkt_tx[d]
    end
    avg_SER = mean(skipmissing(ser_each))
    avg_PER = mean(skipmissing(per_each))
    return avg_SER, avg_PER
end

# ===== SNRスイープ（SERとPERをCSVに保存） =====
function run_sweep_packet_and_save_csv(
    filename::Union{Nothing,String}=nothing;
    save_dir::Union{Nothing,String}=nothing,
    sf::Int=7, bw::Float64=125_000.0,
    snr_min::Float64=-20.0, snr_max::Float64=10.0, snr_step::Float64=2.0,
    iter::Int=500, device_list::Vector{Int} = [1,2,4,8],
    frame_len_syms::Int=4096, pkt_syms::Int=256, wraparound::Bool=false
)
    snrs = collect(snr_min:snr_step:snr_max)
    # 1列目: SNR、以降: 各デバイス数の SER、さらに PER の列も続ける
    ser_mat = zeros(Float64, length(snrs), 1 + length(device_list))
    per_mat = zeros(Float64, length(snrs), 1 + length(device_list))
    ser_mat[:,1] .= snrs
    per_mat[:,1] .= snrs

    for (j, nd) in enumerate(device_list)
        @info "Simulating nd=$nd, pkt_syms=$pkt_syms, frame_len_syms=$frame_len_syms, wrap=$(wraparound)"
        for (i, snr) in enumerate(snrs)
            avg_SER, avg_PER = per_ser_multi_device_time(
                snr, sf, bw, iter;
                num_devices=nd,
                frame_len_syms=frame_len_syms,
                pkt_syms=pkt_syms,
                wraparound=wraparound
            )
            ser_mat[i, j+1] = avg_SER
            per_mat[i, j+1] = avg_PER
            @printf "SNR=%5.1f dB  SER=%.6e  PER=%.6e\n" snr avg_SER avg_PER
        end
    end

    # ファイル名
    devices_str = join(device_list, "-")
    if filename === nothing
        base = @sprintf("PKT_multi_ser_per_sf%d_devices%s_frame%ds_pkt%ds_wrap%s_snr%.1f-%.1f_iter%d.csv",
                        sf, devices_str, frame_len_syms, pkt_syms,
                        wraparound ? "on" : "off", snr_min, snr_max, iter)
        filename = save_dir !== nothing ? joinpath(save_dir, base) : base
    end
    mkpath(dirname(filename))

    # SERとPERを別CSVに保存
    header_ser = ["SNR_dB"; ["SER_nd=$(nd)" for nd in device_list]...]
    header_per = ["SNR_dB"; ["PER_nd=$(nd)" for nd in device_list]...]
    serfile = replace(filename, ".csv" => "_SER.csv")
    perfile = replace(filename, ".csv" => "_PER.csv")

    open(serfile, "w") do io
        println(io, join(header_ser, ","))
        writedlm(io, ser_mat, ',')
    end
    open(perfile, "w") do io
        println(io, join(header_per, ","))
        writedlm(io, per_mat, ',')
    end

    @info "SER を保存: $serfile"
    @info "PER を保存: $perfile"
    return snrs, ser_mat, per_mat
end

snrs, ser_mat, per_mat = run_sweep_packet_and_save_csv(nothing;
    save_dir="/Users/kimparaibuki/Cursor/LoRa_SER/result",
    sf=7, bw=125_000.0,
    snr_min=-20.0, snr_max=0.0, snr_step=0.5,
    iter=3000, device_list=[2,4,8],
    frame_len_syms=4096,   # 1フレームのシンボル数（時間スロット）
    pkt_syms=64,          # 1パケットのシンボル数（前置き・ヘッダ含む等価長として）
    wraparound=false       # 末尾超えは送信スキップ
)