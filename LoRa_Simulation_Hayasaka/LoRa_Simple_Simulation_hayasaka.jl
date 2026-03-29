using Distributed
using StatsBase
using CSV
using DataFrames
using ProgressMeter
using Plots
using Dates
using Random
using Distributions
using Statistics
using LinearAlgebra
import Base.Filesystem.mkdir

# デバッグモードの切り替え
const debug = false

# ===== パラメータ設定 =====
# シミュレーション時間(分)
const SIM_PERIOD = 24 * 60.0 * 3  # 3日間

# エリアサイズ(km)
const area_size = 0.5

# 周波数幅(Hz)
const BW = 125000.0

# チャネル数
const channel_num = 4

# 光速
const c = 3.0 * 10^8

# パスロス係数
const α = 4.0

# 伝搬損失オフセット
const β = 9.5

# キャリアセンス閾値
const carrier_sense_threshold = -110

# 伝搬周波数係数
const γ = 4.5

# 搬送波周波数
const f_c = 923.2

# シャドウイング標準偏差[dB]
const shadowing_standard_deviation = 3.48

# シャドウイング 素波の数
const SW_N = 500

# シャドウイング相関距離[m]
const dcor = 50.0

const noise_figure = 10 # 雑音指数[dB]
const band_width = 125 * 10^3 # 帯域幅
const noise_power_spectrum_density = -174 # 雑音スペクトラム密度

# 送信電力(dBm)
const Tx_dB = 13

# SNRの閾値[SF, threshold]
const SNR_threshold_list = [[6, -5], [7, -7.5], [8, -10], [9, -12.5], [10, -15], [11, -17.5], [12, -20]]

# SIRの閾値
const SIR_threshold_list = [0 -8 -10 -11 -11 -11 -11; -11 0 -11 -13 -14 -14 -14; -14 -13 0 -14 -16 -17 -17; -17 -17 -16 0 -17 -19 -20; -19 -19 -19 -19 0 -20 -22; -22 -22 -22 -22 -22 0 -23; -24 -24 -24 -25 -25 -25 0]

# パケット送信周期(分)
const packet_period = 2.0

# サブフレーム数
const subframe_number = 4

# スロット
const slot_number = 150

# ペイロードサイズ
const payload_bit = 160

# 符号化率
const CR = 4 / 7

# クロックドリフト平均最小
const CD_mean_min = -1.91 * 10^(-3)

# クロックドリフト平均最大
const CD_mean_max = 0.28 * 10^(-3)

# クロックドリフト分散最小
const CD_variance_min = 9.59 * 10^(-11)

# クロックドリフト分散最大
const CD_variance_max = 3.19 * 10^(-10)

# シミュレーション設定 - 修正版
const SIM_NUM = 100  # デバッグ用に減らす
const SIM_start_size = 50  # 開始端末台数
const SIM_particle_size = 50  # 刻み
const number_devices_max = 100  # 最大端末台数

# 現在時刻
const date = string(Dates.format(now(UTC) + Hour(9), "yyyymmdd_HHMMSS"), "_$SIM_NUM", "_$carrier_sense_threshold")

# ===== データ構造定義 =====
# 端末のパラメーター
Base.@kwdef mutable struct N_parameters
    node_id::Int64
    status::Int64 = 0 # 0:待機，1:CS中，2:Txモード，3:Rxモード
    packet_size::Float64
    sf::Int64
    group_id::Int64
    class::String
    x::Float64
    y::Float64
    P_dB::Int64 = Tx_dB
    channel::Int64
    shadowing::Array{Float64} = []
    offset::Float64 = rand() * packet_period
    last_subframe_number::Int64 = 0
    last_backoff::Int64 = 0
    send_bit::Array{Int64} = zeros(Int(floor(log2(channel_num * slot_number))))
    packet_generation_cycle::Float64
    last_packet_slot::Int64 = 0
    packet_generation_time::Float64
    send_packet_success::Int64 = 0
    packet_collision_count::Int64 = 0
    snr_lost_count::Int64 = 0
    snr_lost_count_inslot::Int64 = 0
    packet_collision_count_inslot::Int64 = 0
    last_packet_success::Int64 = 0
    packet_num::Int64 = 0
    packet_lost_count::Int64 = 0
    cs_detection::Int64 = 0
    tx_collision::Array{Bool} = [false, false]
    cs_detection_time::Array{Float64} = []
    cs_detection_type::Array{Array{Bool}} = []
    cs_detection_node::Array{Array{Int64}} = []
    cs_detection_channel::Array{Int64} = []
    clockdrift_mean::Float64 = rand(Uniform(CD_mean_min, CD_mean_max))
    clockdrift_variance::Float64 = rand(Uniform(CD_variance_min, CD_variance_max))
    CD_last_calc::Float64 = 0.0
    CD_accumulated::Float64 = 0.0
    usable_channel::Array{Int64} = [i for i = 1:channel_num]
    packet_num_cdf::Int64 = 0
    packet_success_cdf::Int64 = 0
end

# シャドウイング用構造体
@kwdef mutable struct s_wave
    f_x::Float64 = 0.0
    f_y::Float64 = 0.0
end

@kwdef mutable struct F_T
    f_T::s_wave = s_wave()
    f_R::s_wave = s_wave()
    theta::Float64 = 0.0
end

# 結果構造体
Base.@kwdef mutable struct Result
    node_xy::Array{Array{Float64}}
    collision::Array{Int64}
    success::Array{Int64}
    rost::Array{Int64}
    packet_num::Array{Int64}
    throughput::Array{Int64}
    pdr_cdf
    pdr_cdf1
    pdr_cdf2
    pdr_cdf3
    pdr_cdf4
    node_num::Array{Int64}
    node_pdf
    cdf_throu
    est_sum
    est_true_sum
    est_interf
    only_insys_interf
    only_exsys_interf
    both_interf
    gw_collision_interf
    mis_detection
    backoff_n
    collision_n
    double_n
    mis_n
end

# ===== 基本関数 =====
function distance(x1, y1, x2, y2)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

function SNR_threshold(SF::Int64)
    if SF < 6 || SF > 12
        println("Error of SF")
        sqrt(-1)
    end
    return filter(x -> x[1] == SF, SNR_threshold_list)[1][2]
end

function SIR_threshold(sf_int::Int64, sf_ref::Int64)
    return SIR_threshold_list[(sf_int-6)*7+(sf_ref-5)]
end

function packet_length_type(packet_type::String, SF::Int64)
    if packet_type == "uplink_toa"
        return ((8 + 4.25 + 8) + ceil((payload_bit) / (CR * SF))) * ((2^SF) / BW) / 60
    elseif packet_type == "downlink_toa"
        return ((8 + 4.25 + 8) + ceil((payload_bit) / (CR * SF))) * ((2^SF) / BW) / 60
    end
    println("packet type name error")
    sqrt(-1)
end

function select_sf(group_id::Int64)
    return 7
end

function number_PLIM_bits(channel_num, slot_number)
    return Int(floor(log2(channel_num * slot_number)))
end

function fim(tx_code::Int64, channel::Array)
    binary_array = [i in channel ? 1 : 0 for i in 1:channel_num]
    num_channel = sum(binary_array)
    S = cumsum(binary_array)
    ch = findfirst(isequal(mod(tx_code, num_channel) + 1), S)
    time_slot = Int(floor(tx_code / num_channel))
    return ch, time_slot
end

function send_packet_generator(Node::N_parameters)
    bit_seq = parse(Int, join(string.(Node.send_bit)), base=2)
    tx_code = mod((bit_seq + Node.node_id + Node.packet_num), length(Node.usable_channel) * slot_number)
    ch, slot = fim(tx_code, Node.usable_channel)
    tx_time = Node.offset + Node.CD_accumulated + (Node.packet_num - 1) * packet_period + 5 / 60000 + slot * (packet_period / (subframe_number * slot_number))
    return ch, slot, tx_time
end

function backoff(Node::N_parameters, sub_num::Int64)
    bit_seq = parse(Int, join(string.(Node.send_bit)), base=2)
    tx_code = mod((bit_seq + Node.node_id + Node.packet_num + sub_num), length(Node.usable_channel) * slot_number)
    ch, slot = fim(tx_code, Node.usable_channel)
    tx_time = Node.offset + Node.CD_accumulated + (Node.packet_num - 1) * packet_period + 5 / 60000 + slot * (packet_period / (subframe_number * slot_number)) + (packet_period / subframe_number * sub_num)
    return ch, slot, tx_time
end

function node_generator(Node_all::Array, num_device_all::Int64, group_id::Int64, packet_type::String)
    for i = 1:num_device_all
        theta::Float64 = rand() * 2 * pi
        r::Float64 = sqrt(rand()) * area_size
        Node_all[i] = N_parameters(node_id=i, packet_size=packet_length_type(packet_type, select_sf(group_id)), channel=rand(1:channel_num), sf=select_sf(group_id), group_id=group_id, class="A", packet_generation_cycle=packet_period, packet_generation_time=0.0, x=r * cos(theta), y=r * sin(theta))
    end
    return Node_all
end

# シャドウイング関数
function shadowing(x, y, u, v, waves)
    c_n::Float64 = sqrt(2.0 / SW_N)
    shadowing_val::Float64 = 0.0
    
    for i in 1:SW_N
        shadowing_val += c_n * cos(2 * π * (waves[i].f_T.f_x * x + waves[i].f_T.f_y * y + waves[i].f_R.f_x * u + waves[i].f_R.f_y * v) + waves[i].theta)
    end
    
    return shadowing_val * shadowing_standard_deviation
end

function generate_waves(waves)
    for i in 1:Int64(SW_N / 2)
        waves[i] = F_T()
        waves[i].f_T = generate_ft()
        waves[i].f_R = generate_ft()
        waves[i].theta = generate_theta()
        
        waves[i+Int64(SW_N / 2)] = F_T()
        waves[i+Int64(SW_N / 2)].f_T = waves[i].f_R
        waves[i+Int64(SW_N / 2)].f_R = waves[i].f_T
        waves[i+Int64(SW_N / 2)].theta = waves[i].theta
    end
    return
end

function generate_ft()
    f = s_wave()
    f_t::Float64 = 0.0
    α = log(2.0) / dcor
    φ = rand() * 2 * π
    
    f_t = abs(α / (2 * π) * (sqrt(1 / ((1 - rand())^2)) - 1))
    
    f.f_x = f_t * cos(φ)
    f.f_y = f_t * sin(φ)
    
    return f
end

function generate_theta()
    return rand() * 2 * π
end

function recieve_power(node_id::Int64, Node_now_coordinate::Array, Node_pre_coordinate::Array, power::Int64, waves::Array, return_dB::Int64, no_shadowing::Int64)
    Node_dis::Float64 = sqrt((Node_now_coordinate[1] - Node_pre_coordinate[1])^2 + (Node_now_coordinate[2] - Node_pre_coordinate[2])^2)
    Passloss = 10 * α * log10(Node_dis) + β + 10 * γ * log10(f_c)
    shadowing_value = shadowing(Node_now_coordinate[1] * 1000, Node_now_coordinate[2] * 1000, Node_pre_coordinate[1] * 1000, Node_pre_coordinate[2] * 1000, waves)
    if no_shadowing == 1
        shadowing_value = 0.0
    end
    if return_dB == 1
        return power - Passloss - shadowing_value
    else
        return 10^((power - Passloss - shadowing_value) / 10)
    end
end

function SNR_judge(node::N_parameters, waves::Array)
    Node_dis::Float64 = sqrt((node.x - 0.0)^2 + (node.y - 0.0)^2)
    Passloss = 10 * α * log10(Node_dis) + β + 10 * γ * log10(f_c)
    shadowing_value = shadowing(node.x * 1000, node.y * 1000, 0.0, 0.0, waves)
    if node.P_dB - Passloss - shadowing_value - (noise_power_spectrum_density + 10 * log10(band_width) + noise_figure) > SNR_threshold(node.sf)
        return 0
    else
        return 1
    end
end

function carrier_sense(node_id::Int64, Node_all::Array, waves::Array)
    now_Node_xy = [Node_all[node_id].x, Node_all[node_id].y]
    now_channel = Node_all[node_id].channel
    Node_xy = []
    power::Float64 = 0.0
    power_db::Float64 = 0.0
    
    for i in Node_all
        if i.channel == now_channel && i.status == 2
            push!(Node_xy, ([i.x, i.y], i.P_dB, i.node_id))
        end
    end
    
    if isempty(Node_xy)
        return 0
    end
    
    for j in eachindex(Node_xy)
        power += recieve_power(Node_xy[j][3], now_Node_xy, [Node_xy[j][1][1], Node_xy[j][1][2]], Node_xy[j][2], waves, 0, 0)
    end
    power_db = 10 * log10(power)
    
    if power_db >= carrier_sense_threshold
        return 1
    else
        return 0
    end
end

function collision_judge(node_id::Int64, Node_all::Array, waves::Array)
    sending = []
    nowchannel = Node_all[node_id].channel
    
    for node in Node_all
        if node.channel == nowchannel && node.status == 2 && node.node_id != node_id
            push!(sending, (node.packet_generation_time, node.node_id))
        end
    end
    
    if isempty(sending)
        return 0
    end
    
    push!(sending, (Node_all[node_id].packet_generation_time, node_id))
    sort!(sending, by=x -> x[1])
    
    for (index, value) in enumerate(sending)
        if value[2] < 10000 && index == 1
            if (sending[index+1][1] - sending[index][1]) >= (20.25 * ((2^Node_all[value[2]].sf) / BW) / 60)
                Node_all[value[2]].tx_collision = [false, true]
            elseif (sending[index+1][1] - sending[index][1]) < (20.25 * ((2^Node_all[value[2]].sf) / BW) / 60)
                Node_all[value[2]].tx_collision = [true, true]
            end
        elseif value[2] < 10000 && index != 1
            Node_all[value[2]].tx_collision = [true, true]
        end
    end
end

function calc_clockdrift(last_time::Float64, now_time::Float64, CD_mean::Float64, CD_variance::Float64)
    elapsed_time = Int(round(now_time - last_time))
    if elapsed_time <= 0
        return 0, 1
    end
    return sum(rand(Normal(CD_mean, CD_variance), elapsed_time)) / 60, 0
end

# ===== メインシミュレーション関数 =====
function main_loop(sim_index::Int64, number_devices::Int64)
    # 結果書き出し用変数
    collision_count = 0
    backoff_count = 0
    mis_count = 0
    double_count = 0
    collision_packet = zeros(Int(ceil(SIM_PERIOD / packet_period)))
    success_packet = zeros(Int(ceil(SIM_PERIOD / packet_period)))
    rost_packet = zeros(Int(ceil(SIM_PERIOD / packet_period)))
    packet_num = zeros(Int(ceil(SIM_PERIOD / packet_period)))
    throughput = zeros(Int(ceil(SIM_PERIOD / packet_period)))
    cdf = Array{Array{Float64}}(undef, 0)
    cdf_1 = Array{Array{Float64}}(undef, 0)
    cdf_2 = Array{Array{Float64}}(undef, 0)
    cdf_3 = Array{Array{Float64}}(undef, 0)
    cdf_4 = Array{Array{Float64}}(undef, 0)
    throughput_cdf = Array{Array{Float64}}(undef, 0)
    throughput_cdf_calc = zeros(number_devices)
    result_node = []
    result_pdf = []
    error = 0
    no_error = 0
    
    # 時刻
    now_time::Float64 = 0.0
    event_time = Array{Tuple{Float64,Int64,Int64}}(undef, 0) # (時刻，コマンド，EN番号)
    gw_time = [Vector{Tuple{Float64,Float64,Int64}}() for _ in 1:channel_num] # (送信開始時刻，送信終了時刻，端末番号)
    
    # 周波数チャネル
    now_channel::Int64 = 0
    
    # ノード確保配列
    Node_all = Array{N_parameters}(undef, number_devices)
    
    # ノード生成
    node_generator(Node_all, number_devices, 0, "uplink_toa")
    
    # シャドウイング値生成
    waves = Array{F_T}(undef, SW_N)
    generate_waves(waves)
    
    # 初期パケット生成
    for node_index = 1:number_devices
        Node_all[node_index].packet_num += 1
        Node_all[node_index].channel, Node_all[node_index].last_packet_slot, packet_generation_time = send_packet_generator(Node_all[node_index])
        Node_all[node_index].last_subframe_number += 1
        Node_all[node_index].packet_size = ((8 + 4.25 + 8) + ceil((160) / Node_all[node_index].sf)) * ((2^Node_all[node_index].sf) / BW) / 60
        
        # クロックドリフト計算
        CD, time_error = calc_clockdrift(Node_all[node_index].CD_last_calc * 60, packet_generation_time * 60, Node_all[node_index].clockdrift_mean, Node_all[node_index].clockdrift_variance)
        if time_error == 0
            no_error += 1
        else
            error += 1
        end
        
        Node_all[node_index].packet_generation_time = packet_generation_time + CD
        Node_all[node_index].CD_accumulated += CD
        Node_all[node_index].CD_last_calc = Node_all[node_index].packet_generation_time
        push!(event_time, (Node_all[node_index].packet_generation_time - 5 / 60000, 1, node_index))
        push!(event_time, (Node_all[node_index].packet_generation_time, 2, node_index))
    end
    
    # SNR規範により各端末のSF決定
    for node_dis = 1:number_devices
        SNR = recieve_power(node_dis, [Node_all[node_dis].x, Node_all[node_dis].y], [0.0, 0.0], Node_all[node_dis].P_dB, waves, 1, 0) - (noise_power_spectrum_density + 10 * log10(band_width) + noise_figure)
        if SNR > -7.5
            Node_all[node_dis].sf = 7
        elseif -7.5 >= SNR > -10.0
            Node_all[node_dis].sf = 8
        elseif -10.0 >= SNR > -12.5
            Node_all[node_dis].sf = 9
        elseif -12.5 >= SNR
            Node_all[node_dis].sf = 10
        end
        Node_all[node_dis].packet_size = packet_length_type("uplink_toa", Node_all[node_dis].sf)
    end
    
    # event_timeをソート
    sort!(event_time, by=x -> x[1])
    
    # 各種カウント変数定義
    counter_gw::Int64 = 0
    counter_en::Int64 = 0
    counter_bk::Int64 = 0
    
    # シミュレーション本体
    while !isempty(event_time)
        # event_timeをソート
        sort!(event_time, by=x -> x[1])
        node_id::Int64 = event_time[1][3]
        now_channel = Node_all[node_id].channel
        now_time = event_time[1][1]
        now_slot = Int(ceil(now_time / packet_period))
        
        if now_slot > Int(ceil(SIM_PERIOD / packet_period))
            break
        end
        
        send_success::Bool = false
        
        # イベント操作
        if event_time[1][2] == 1 ## キャリアセンス開始
            Node_all[node_id].status = 1
            if carrier_sense(node_id, Node_all, waves) == 1
                Node_all[node_id].cs_detection = 1
                push!(Node_all[node_id].cs_detection_time, now_time)
            end
            popfirst!(event_time)
            continue
            
        elseif event_time[1][2] == 2 ## キャリアセンス終了
            Node_all[node_id].status = 0
            if Node_all[node_id].cs_detection == 1 && Node_all[node_id].last_subframe_number == subframe_number
                # パケット破棄処理
                popfirst!(event_time)
                packet_num[now_slot] += 1
                rost_packet[now_slot] += 1
                counter_en += 1
            elseif Node_all[node_id].cs_detection == 1 && Node_all[node_id].last_subframe_number != subframe_number
                # バックオフ処理
                counter_bk += 1
                Node_all[node_id].channel, Node_all[node_id].last_packet_slot, packet_generation_time = backoff(Node_all[node_id], Node_all[node_id].last_subframe_number)
                Node_all[node_id].last_subframe_number += 1
                
                # クロックドリフト計算
                CD, time_error = calc_clockdrift(Node_all[node_id].CD_last_calc * 60, packet_generation_time * 60, Node_all[node_id].clockdrift_mean, Node_all[node_id].clockdrift_variance)
                if time_error == 0
                    no_error += 1
                else
                    error += 1
                end
                
                Node_all[node_id].packet_generation_time = packet_generation_time + CD
                Node_all[node_id].CD_accumulated += CD
                Node_all[node_id].CD_last_calc = Node_all[node_id].packet_generation_time
                push!(event_time, (Node_all[node_id].packet_generation_time - 5 / 60000, 1, node_id))
                push!(event_time, (Node_all[node_id].packet_generation_time, 2, node_id))
                popfirst!(event_time)
                Node_all[node_id].cs_detection = 0
                continue
            else
                # 送信開始
                counter_gw += 1
                Node_all[node_id].status = 2
                push!(event_time, (event_time[1][1] + Node_all[node_id].packet_size, 4, node_id))
                Node_all[node_id].cs_detection = 0
                popfirst!(event_time)
                # 衝突判定
                collision_judge(node_id, Node_all, waves)
                continue
            end
            
        elseif event_time[1][2] == 4
            # 送信終了
            Node_all[node_id].status = 0
            if SNR_judge(Node_all[node_id], waves) == 0 && Node_all[node_id].tx_collision != [true, true]
                if Node_all[node_id].tx_collision == [false, false]
                    send_success = true
                    packet_num[now_slot] += 1
                    success_packet[now_slot] += 1
                    throughput[now_slot] += length(Node_all[node_id].send_bit) + payload_bit
                    push!(gw_time[now_channel], (Node_all[node_id].packet_generation_time, Node_all[node_id].packet_generation_time + Node_all[node_id].packet_size, node_id))
                end
            elseif SNR_judge(Node_all[node_id], waves) == 0 && Node_all[node_id].tx_collision == [true, true]
                packet_num[now_slot] += 1
                collision_packet[now_slot] += 1
            else
                packet_num[now_slot] += 1
            end
            popfirst!(event_time)
        else
            println("予期しないイベント番号です")
            println(event_time[1][2])
            sqrt(-1)
        end
        
        # 次のパケット生成
        Node_all[node_id].send_bit = rand(0:1, length(Node_all[node_id].send_bit))
        Node_all[node_id].packet_num += 1
        Node_all[node_id].channel, Node_all[node_id].last_packet_slot, packet_generation_time = send_packet_generator(Node_all[node_id])
        Node_all[node_id].last_subframe_number = 1
        
        # クロックドリフト計算
        CD, time_error = calc_clockdrift(Node_all[node_id].CD_last_calc * 60, packet_generation_time * 60, Node_all[node_id].clockdrift_mean, Node_all[node_id].clockdrift_variance)
        if time_error == 0
            no_error += 1
        else
            error += 1
        end
        
        Node_all[node_id].packet_generation_time = packet_generation_time + CD
        Node_all[node_id].CD_accumulated += CD
        Node_all[node_id].CD_last_calc = Node_all[node_id].packet_generation_time
        push!(event_time, (Node_all[node_id].packet_generation_time - 5 / 60000, 1, node_id))
        push!(event_time, (Node_all[node_id].packet_generation_time, 2, node_id))
        
        # 初期化
        Node_all[node_id].cs_detection = 0
        Node_all[node_id].tx_collision = [false, false]
        if send_success
            Node_all[node_id].last_packet_success = 1
        else
            Node_all[node_id].last_packet_success = 0
        end
        Node_all[node_id].cs_detection_time = []
        Node_all[node_id].cs_detection_type = []
        Node_all[node_id].cs_detection_node = []
        Node_all[node_id].cs_detection_channel = []
    end
    
    return Result([getfield.(Node_all, :x), getfield.(Node_all, :y)], collision_packet, success_packet, rost_packet, packet_num, throughput, cdf, cdf_1, cdf_2, cdf_3, cdf_4, result_node, result_pdf, throughput_cdf, 0, 0, zeros(Int(ceil(SIM_PERIOD / packet_period))), zeros(Int(ceil(SIM_PERIOD / packet_period))), zeros(Int(ceil(SIM_PERIOD / packet_period))), zeros(Int(ceil(SIM_PERIOD / packet_period))), zeros(Int(ceil(SIM_PERIOD / packet_period))), zeros(Int(ceil(SIM_PERIOD / packet_period))), backoff_count, collision_count, double_count, mis_count)
end

# ===== 結果出力関数 =====
function write_csv(data_type::Int64, result_array::Array, device_num, plot_type::Int64, sim_num::Int64)
    data1 = []
    data2 = []
    result = []
    type::String = "nothing"
    
    if data_type == 1
        push!(data1, sum(getfield.(result_array, :success)) ./ sim_num)
        push!(data2, sum(getfield.(result_array, :packet_num)) ./ sim_num)
        result = [data1[i] ./ data2[i] for i in eachindex(data1)]
        type = "PDR"
    elseif data_type == 2
        push!(data1, sum(getfield.(result_array, :collision)) ./ sim_num)
        push!(data2, sum(getfield.(result_array, :packet_num)) ./ sim_num)
        result = [data1[i] ./ data2[i] for i in eachindex(data1)]
        type = "Collision"
    elseif data_type == 3
        push!(data1, sum(getfield.(result_array, :rost)) ./ sim_num)
        push!(data2, sum(getfield.(result_array, :packet_num)) ./ sim_num)
        result = [data1[i] ./ data2[i] for i in eachindex(data1)]
        type = "Rost"
    elseif data_type == 4
        push!(data1, sum(getfield.(result_array, :throughput)) ./ sim_num)
        result = data1 / (packet_period * 60.0)
        type = "Throughput"
    end
    
    # DataFrameを作成
    df = DataFrame()
    if plot_type == 1
        df[!, :x] = collect(packet_period:packet_period:SIM_PERIOD)
        for (i, col_data) in enumerate(result)
            col_name = "$device_num"
            df[!, col_name] = col_data
        end
        # CSVファイルに書き出し
        if !isdir("result")
            mkdir("result")
        end
        df |> CSV.write("result/$type-$device_num.csv", delim=',', writeheader=false)
        println("結果を保存しました: result/$type-$device_num.csv")
    end
end

# ===== メイン実行部分 =====
println("LoRaシミュレーション開始...")
println("設定:")
println("  シミュレーション回数: $SIM_NUM")
println("  端末台数範囲: $SIM_start_size から $number_devices_max まで ($SIM_particle_size 刻み)")
println("  シミュレーション時間: $SIM_PERIOD 分")

if !isdir("result")
    mkdir("result")
    println("resultフォルダを作成しました")
end

runtime = @elapsed begin
    for number_devices = SIM_start_size:SIM_particle_size:number_devices_max
        println("******* start number_devices=$number_devices *******")
        if debug
            result = [main_loop(1, number_devices)]
            break
        else
            result = [main_loop(i, number_devices) for i in 1:SIM_NUM]
            length_of_results = length(result)
            println("シミュレーション成功回数:$length_of_results")
        end
        
        if !debug
            write_csv(1, result, number_devices, 1, length_of_results) ## 1:PDR
            write_csv(2, result, number_devices, 1, length_of_results) ## 2:衝突率
            write_csv(3, result, number_devices, 1, length_of_results) ## 3:破棄率
            write_csv(4, result, number_devices, 1, length_of_results) ## 4:スループット
        end
        GC.gc()
        println("******* end number_devices=$number_devices *******\n")
    end
end
println("実行時間: $runtime 秒")
println("シミュレーション完了！")
