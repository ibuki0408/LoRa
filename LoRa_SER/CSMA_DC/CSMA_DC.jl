using Random, Distributions, Statistics, DelimitedFiles

# 送信の様子を可視化したコード
# 端末インデックスと送信時刻，送信成功/失敗をCSVに出力

# =============================
# パラメータ設定
# =============================
const SF   = 7                  # spreading factor
const BW   = 125e3              # 帯域幅 (Hz)
const pkt_range = (64, 128)       # パケット長シンボル数の範囲
const sim_time = 1.0           # シミュレーション時間 [s]
const λ = 2.0                   # Poisson到着率 [pkt/s per device]
const num_devices = 1000        # 端末数
const DC = 0.01                 # Duty Cycle 制約 (1%)
const backoff_max = 0.01        # 最大バックオフ時間 [s]
const output_csv = "LoRa_SER/CSMA_DC/csma_dc_SF$(SF)_num$(num_devices)_pkt$(pkt_range[1])-$(pkt_range[2])_results.csv"  # CSV出力先

# シンボル時間
ts = (2.0^SF) / BW

# =============================
# デバイス構造体
# =============================
mutable struct Device
    id::Int
    next_time::Float64    # 次に送信可能な時刻（DC考慮）
end

# =============================
# パケットイベント構造体
# =============================
struct Packet
    dev::Int
    start::Float64
    duration::Float64
    success::Bool
end

# =============================
# シミュレーション本体
# =============================
function simulate_csma_poisson_csv()
    rng = MersenneTwister(0)

    # デバイス初期化
    devices = [Device(d, 0.0) for d in 1:num_devices]
    arrivals = [rand(rng, Exponential(1/λ)) for _ in 1:num_devices]  # 初期到着

    # 結果記録
    all_pkts = Packet[]
    total = 0
    success_count = 0

    t = 0.0
    while t < sim_time
        # 次に送信するデバイスを決定
        next_dev = argmin(arrivals)       # 最小インデックスだけ取得
        # そのデバイスの到着時刻
        t = arrivals[next_dev]            # そのデバイスの到着時刻

        if t > sim_time
            break
        end

        # パケット長（シンボル数 → 時間）
        pkt_syms = rand(rng, pkt_range[1]:pkt_range[2])
        duration = pkt_syms * ts

        # DC制約で送信不可なら次到着に置き換え
        if t < devices[next_dev].next_time
            arrivals[next_dev] = devices[next_dev].next_time + rand(rng, Exponential(1/λ))
            continue
        end

        # チャネル占有チェック
        busy = any(pkt -> (t < pkt.start + pkt.duration && t + duration > pkt.start), all_pkts)

        if busy
            # バックオフ
            t_backoff = t + rand(rng) * backoff_max
            arrivals[next_dev] = t_backoff
            continue
        else
            # 送信開始
            overlap = count(pkt -> (pkt.dev != next_dev && max(pkt.start, t) < min(pkt.start+pkt.duration, t+duration)), all_pkts)
            pkt_success = (overlap == 0)

            push!(all_pkts, Packet(next_dev, t, duration, pkt_success))
            total += 1
            if pkt_success
                success_count += 1
            end

            # DC制約適用
            devices[next_dev].next_time = t + duration / DC

            # 次の到着イベント生成
            arrivals[next_dev] = t + rand(rng, Exponential(1/λ))
        end
    end

    # SER計算
    ser = (total - success_count) / total

    # CSV出力
    header = ["Device", "StartTime", "Duration", "Success"]
    open(output_csv, "w") do io
        println(io, join(header, ","))
        for pkt in all_pkts
            println(io, join([pkt.dev, pkt.start, pkt.duration, pkt.success], ","))
        end
    end

    println("=== Simulation Result ===")
    println("Total packets: $total")
    println("Success: $success_count")
    println("SER: $ser")
    println("Results saved to $output_csv")

    return ser, total, success_count
end

# =============================
# 実行
# =============================
ser, total, success_count = simulate_csma_poisson_csv()
