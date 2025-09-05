#include("param.jl")
#include("Node.jl")

#シャドウイング生成
function gen_multivariate_normal(cov)
  L = cholesky(cov).L
  n = size(cov, 1)
  z = randn(n)
  return L * z
end

# function gen_multivariate_normal(cov)
#   L = cholesky(cov).L
#   n = size(cov, 1)
#   z = randn(n, n)
#   return L * z
# end

function distance(x1, y1, x2, y2)
  return sqrt.((x1 .- x2) .^ 2 .+ (y1 .- y2) .^ 2)
end

function distance_one(x1, y1, x2, y2)
  return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

function gen_varcov_matrix(x, y, dcor, sigma)
  dmat = distance.(x, y, x', y')
  tmp = log(2.0) / dcor
  return sigma .* sigma .* exp.(-dmat .* tmp)
end

#配列平坦化
function chain_from_iterable(iterables)
  return (x for iterable in iterables for x in iterable)
end

#各時間長整合判定
function time_matching(SIM_PERIOD::Float64, packet_period::Float64, grouping_period::Float64)
  if SIM_PERIOD % packet_period != 0.0 || SIM_PERIOD % grouping_period != 0.0
    print("Error of SIM_PERIOD")
    ##これで強制プログラム終了
    sqrt(-1)
  end
  return 0
end

#SNRの閾値
function SNR_threshold(SF::Int64)
  if SF < 6 || SF > 12
    println("Error of SF")
    ##これで強制プログラム終了
    sqrt(-1)
  end
  return filter(x -> x[1] == SF, SNR_threshold_list)[1][2]
end

#SIRの閾値
function SIR_threshold(sf_int::Int64, sf_ref::Int64)
  return SIR_threshold_list[(sf_int-6)*7+(sf_ref-5)]
end

#SIR制約により復調可能か判定(1の場合復調可能，0の場合不可能)
function SIR_judge(sending_packets::Array, Node_all::Array, Interference_all::Array, waves::Array)
  # return 1
  power_int = 0.0
  power = 0.0
  sf = []

  for i in 2:length(sending_packets)
    if sending_packets[i][2] >= 10000
      power_int += recieve_power(sending_packets[i][2], [Interference_all[sending_packets[i][2]%10000].x, Interference_all[sending_packets[i][2]%10000].y], [0.0, 0.0], Interference_all[sending_packets[i][2]%10000].P_dB, waves, 0, 0)
      push!(sf, Interference_all[sending_packets[i][2]%10000].sf)
    else
      power_int += recieve_power(sending_packets[i][2], [Node_all[sending_packets[i][2]].x, Node_all[sending_packets[i][2]].y], [0.0, 0.0], Node_all[sending_packets[i][2]].P_dB, waves, 0, 0)
      push!(sf, Node_all[sending_packets[i][2]].sf)
    end
  end
  power = recieve_power(sending_packets[1][2], [Node_all[sending_packets[1][2]].x, Node_all[sending_packets[1][2]].y], [0.0, 0.0], Node_all[sending_packets[1][2]].P_dB, waves, 1, 0) - 10 * log10(power_int)

  for i in eachindex(sf)
    if power <= SIR_threshold(sf[i], Node_all[sending_packets[1][2]].sf)
      return 0
    end
  end
  return 1
end

function number_PLIM_bits(channel_num, slot_number)
  return Int(floor(log2(channel_num * slot_number)))
end


#SNR制約により復調可能か判定(0の場合復調可能，1の場合不可能)
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

#パケット長判定
function packet_length_type(packet_type::String, SF::Int64)
  if packet_type == "uplink_toa"
    return ((8 + 4.25 + 8) + ceil((payload_bit) / (CR * SF))) * ((2^SF) / BW) / 60
  elseif packet_type == "downlink_toa"
    return ((8 + 4.25 + 8) + ceil((payload_bit) / (CR * SF))) * ((2^SF) / BW) / 60
  end
  println("packet type name error")
  ##これで強制プログラム終了
  sqrt(-1)
end

#送信パケット生成
# function send_packet_generator(packet_cycle::Float64, start_time::Float64, end_time::Float64, class::String, Node::N_parameters)
#   if (end_time - start_time) < 0.0
#     print("error of send packet generation")
#     ##これで強制プログラム終了
#     sqrt(-1)
#   elseif (end_time - start_time) == 0.0
#     return []
#   else
#     number_packet::Int64 = fld((end_time - start_time), packet_cycle)
#     packet_time = Array{Float64}(undef, number_packet)
#     for k = 1:number_packet
#       packet_time[k] = start_time + packet_cycle * (k - 1) + (packet_cycle - Node.packet_size - 5 / 60000) * rand()
#     end
#     return packet_time
#   end
# end

# function fim_old(code::Int64, channel::Array)
#   number_channel = length(channel)
#   ch = channel[mod(code, number_channel)+1]
#   time_slot = mod(Int(floor(code / number_channel)), slot_number)
#   return ch, time_slot
# end

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

#ランダムバックオフ
function backoff(Node::N_parameters, sub_num::Int64)
  # bit_seq = mod(parse(Int, join(string.(Node.send_bit)), base=2) + 1, 2^(length(Node.send_bit)))
  bit_seq = parse(Int, join(string.(Node.send_bit)), base=2)
  tx_code = mod((bit_seq + Node.node_id + Node.packet_num + sub_num), length(Node.usable_channel) * slot_number)
  ch, slot = fim(tx_code, Node.usable_channel)
  tx_time = Node.offset + Node.CD_accumulated + (Node.packet_num - 1) * packet_period + 5 / 60000 + slot * (packet_period / (subframe_number * slot_number)) + (packet_period / subframe_number * sub_num)
  return ch, slot, tx_time
end

#SF選択
function select_sf(group_id::Int64)
  return 7
end

#ノード生成
function node_generator(Node_all::Array, num_device_all::Int64, group_id::Int64, packet_type::String)
  for i = 1:num_device_all
    theta::Float64 = rand() * 2 * pi
    r::Float64 = sqrt(rand()) * area_size
    Node_all[i] = N_parameters(node_id=i, packet_size=packet_length_type(packet_type, select_sf(group_id)), channel=rand(1:channel_num), sf=select_sf(group_id), group_id=group_id, class="A", packet_generation_cycle=packet_period, packet_generation_time=0.0, x=r * cos(theta), y=r * sin(theta))
  end
  return Node_all
end


#衝突判定
function collision_judge(node_id::Int64, Node_all::Array, Interference_all::Array, wave::Array, device_interf::Bool)
  sending = []
  nowchannel = 0
  if device_interf
    nowchannel = Interference_all[node_id%10000].channel
  else
    nowchannel = Node_all[node_id].channel
  end
  for node in Node_all
    if node.channel == nowchannel && node.status == 2 && node.node_id != node_id
      push!(sending, (node.packet_generation_time, node.node_id))
    end
  end
  for node in Interference_all
    if node.channel == nowchannel && node.status == 2 && node.node_id != node_id
      push!(sending, (node.packet_generation_time, node.node_id))
    end
  end
  if isempty(sending)
    return 0
  end
  if device_interf
    push!(sending, (Interference_all[node_id%10000].packet_generation_time, node_id))
  else
    push!(sending, (Node_all[node_id].packet_generation_time, node_id))
  end
  sort!(sending, by=x -> x[1])
  # if Node_all[node_id].packet_generation_time > 2880.0
  #   for i in sending
  #     println(sqrt((Node_all[node_id].x - Node_all[i[2]].x)^2 + (Node_all[node_id].y - Node_all[i[2]].y)^2))
  #     println()
  #   end
  # end
  for (index, value) in enumerate(sending)
    if value[2] < 10000 && index == 1
      if SIR_judge(sending, Node_all, Interference_all, wave) == 0 && (sending[index+1][1] - sending[index][1]) >= (20.25 * ((2^Node_all[value[2]].sf) / BW) / 60)
        Node_all[value[2]].tx_collision = [false, true]
        # println("NG1")
        # else
        #   println("OK")
      elseif SIR_judge(sending, Node_all, Interference_all, wave) == 0 && (sending[index+1][1] - sending[index][1]) < (20.25 * ((2^Node_all[value[2]].sf) / BW) / 60)
        Node_all[value[2]].tx_collision = [true, true]
        #   println("NG2")
        # else
        #   println("OK")
      end
    elseif value[2] < 10000 && index != 1
      Node_all[value[2]].tx_collision = [true, true]
    elseif value[2] >= 10000
      Interference_all[value[2]%10000].tx_collision = 1
    else
      println("error of collision_judge")
    end
  end
end

#受信信号電力
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

#キャリアセンス
function carrier_sense(node_id::Int64, Node_all::Array, Interference_all::Array, shadowing_waves::Array)
  now_Node_xy = []
  now_channel = 0
  if node_id >= 10000
    now_Node_xy = [Interference_all[node_id%10000].x, Interference_all[node_id%10000].y]
    now_channel = Interference_all[node_id%10000].channel
  else
    now_Node_xy = [Node_all[node_id].x, Node_all[node_id].y]
    now_channel = Node_all[node_id].channel
  end
  Node_xy = []
  power::Float64 = 0.0
  power_db::Float64 = 0.0
  for i in Node_all
    if i.channel == now_channel && i.status == 2
      push!(Node_xy, ([i.x, i.y], i.P_dB, i.node_id))
    end
  end
  for i in Interference_all
    if i.channel == now_channel && i.status == 2
      push!(Node_xy, ([i.x, i.y], i.P_dB, i.node_id))
    end
  end
  if isempty(Node_xy)
    return 0
  end
  for j in eachindex(Node_xy)
    power += recieve_power(Node_xy[j][3], now_Node_xy, [Node_xy[j][1][1], Node_xy[j][1][2]], Node_xy[j][2], shadowing_waves, 0, 0)
  end
  power_db = 10 * log10(power)
  if power_db >= carrier_sense_threshold
    if node_id < 10000
      push!(Node_all[node_id].cs_detection_node, getfield.(Node_xy, :3))
      push!(Node_all[node_id].cs_detection_type, [any(x -> x < 10000, getfield.(Node_xy, :3)), any(x -> x >= 10000, getfield.(Node_xy, :3))])
      push!(Node_all[node_id].cs_detection_channel, now_channel)
    end
    return 1
  else
    return 0
  end
end


#パケット長
function packet_size(sender::String, counter::Int64, GW_send::Array, Node_all::Array)
  if sender == "GW"
    return GW_send[counter].packet_size
  elseif sender == "EN"
    return Node_all[counter].packet_size
  end
end

function estimation(gw_time::Array, Node::N_parameters, interference_time::Array, Interference_all::Array, waves, Node_all::Array, gw_collision_time::Array, obs_CD::Array, now_time::Float64)
  result = []

  est_interference::Int64 = 0
  only_insys_interference::Int64 = 0
  only_exsys_interference::Int64 = 0
  both_interference::Int64 = 0
  collision_count::Int64 = 0
  mis_detection::Int64 = 0
  backoff_channel_set = zeros(Int, channel_num)
  double_count::Int64 = 0
  if isempty(Node.cs_detection_time)
    return 0, 0, 0, 0, 0, 0, 0, 0, 0
  end
  #時刻ずれ補正
  mean_time_delay = mean(obs_CD[Node.node_id])

  tmp_plim_seq = mod((Node.last_packet_slot * channel_num + Node.channel - Node.node_id - Node.packet_num - Node.last_subframe_number), (channel_num * slot_number))
  bit_seq = parse(Int, join(string.(Node.send_bit)), base=2)
  if tmp_plim_seq != bit_seq
    println(tmp_plim_seq)
    println(bit_seq)
    println("error of plim seq")
    # sqrt(-1)
  else
    # println("ok")
  end
  tx_code = mod((tmp_plim_seq + Node.node_id + Node.packet_num + Node.last_subframe_number), channel_num * slot_number)
  ch, slot = fim(tx_code, Node.usable_channel)
  ToA = ((8 + 4.25 + 8) + ceil((payload_bit) / (CR * Node.sf))) * ((2^Node.sf) / BW) / 60
  frame_start = now_time - ToA - slot * (packet_period / (subframe_number * slot_number)) - (Node.last_subframe_number - 1) * (packet_period / subframe_number)
  # println("last_subframe_number=", Node.last_subframe_number)
  # println("cs_detection_time=", Node.cs_detection_time)
  for subf_i in 0:(Node.last_subframe_number-2)

    tmp_send_X = (tmp_plim_seq + Node.node_id + Node.packet_num + subf_i) % (channel_num * slot_number)
    tmp_slot = floor(tmp_send_X / channel_num)
    tmp_channel = mod(tmp_send_X, channel_num) + 1
    tmp_index = length(gw_time[tmp_channel])

    #推定時間補正項導出
    slot_diff = Node.last_subframe_number * slot_number + slot - (subf_i * slot_number + tmp_slot)
    # est_time_delay = (slot_diff * (packet_period / (subframe_number * slot_number))) * mean_time_delay
    est_time_delay = 0.0
    #クロックドリフトの信頼区間を計算
    data_size = length(obs_CD[Node.node_id])
    s = std(obs_CD[Node.node_id], corrected=true)
    t = 0.0
    if data_size > 1
      t = quantile(TDist(data_size - 1), 0.98)
    end
    margin_clock_drift = t * s * sqrt(1 + 1 / data_size)
    # EST_margin = (slot_diff * (packet_period / (subframe_number * slot_number))) * margin_clock_drift
    EST_margin = 0.0

    # println("data_size=$data_size")
    # println("slot_diff=$slot_diff")
    # println("mean_time_delay=$mean_time_delay")
    # println("s=$s")
    # println("t=$t")
    # println("est_time_delay=$est_time_delay")
    # println("EST_margin=$EST_margin")
    # println("margin_clock_drift=$margin_clock_drift")
    # println("1slot_diff", (packet_period / (subframe_number * slot_number)))

    # A = (frame_start + (packet_period / subframe_number) * subf_i + (packet_period / (subframe_number * slot_number)) * tmp_slot) - 5 / 60000 - est_time_delay - EST_margin
    # B = (frame_start + (packet_period / subframe_number) * subf_i + (packet_period / (subframe_number * slot_number)) * (tmp_slot)) + ToA - est_time_delay + EST_margin
    # C = Node.cs_detection_time[subf_i+1]
    # println("A=$A")
    # println("B=$B")
    # println("C=$C")
    # println()
    count_flag = 0
    tmp_result = []
    est_time_median = (frame_start + (packet_period / subframe_number) * subf_i + (packet_period / (subframe_number * slot_number)) * (tmp_slot)) + (ToA / 2) - est_time_delay
    diff_est_time_median = []

    while true
      tmp_index -= 1
      if tmp_index == 0
        break
      end
      if gw_time[tmp_channel][tmp_index][2] < (frame_start + (packet_period / subframe_number) * subf_i + (packet_period / (subframe_number * slot_number)) * tmp_slot) - 5 / 60000 - est_time_delay - EST_margin
        break
      end
      if gw_time[tmp_channel][tmp_index][2] < (frame_start + (packet_period / subframe_number) * subf_i + (packet_period / (subframe_number * slot_number)) * (tmp_slot)) + ToA - est_time_delay + EST_margin
        # push!(tmp_result, gw_time[tmp_channel][tmp_index][3])
        push!(tmp_result, (gw_time[tmp_channel][tmp_index][2] - (ToA / 2), gw_time[tmp_channel][tmp_index][3]))
        # collision_count += 1
        # if count_flag == 1
        #   # double_count += 1
        # end
        count_flag = 1

      end
    end
    if length(tmp_result) == 1
      push!(result, tmp_result[1][2])
      collision_count += 1
      if !(tmp_result[1][2] in Node.cs_detection_node[subf_i+1])
        mis_detection += 1
      end
    elseif length(tmp_result) > 1
      double_count += 1
      for i in eachindex(tmp_result)
        push!(diff_est_time_median, abs(tmp_result[i][1] - est_time_median))
      end
      push!(result, tmp_result[argmin(diff_est_time_median)][2])
      collision_count += 1
      if !(tmp_result[argmin(diff_est_time_median)][2] in Node.cs_detection_node[subf_i+1])
        mis_detection += 1
      end
    end

  end


  for i = 1:length(Node.cs_detection_time)
    backoff_channel_set[Node.cs_detection_channel[i]] += 1
    # for ii in eachindex(gw_time)
    #   for jj in 2:length(gw_time[ii])
    #     if gw_time[ii][jj][1] - gw_time[ii][jj-1][2] < 0.0
    #       # print("error:")
    #       # print(gw_time[ii][jj][1])
    #       # print("\t")
    #       # println(gw_time[ii][jj-1][2])
    #       println(gw_time)
    #       sqrt(-1)
    #     end
    #   end
    # end
    # for j = length(gw_time[Node.cs_detection_channel[i]]):-1:1
    #   if gw_time[Node.cs_detection_channel[i]][j][1] <= Node.cs_detection_time[i]
    #     if gw_time[Node.cs_detection_channel[i]][j][2] >= Node.cs_detection_time[i]
    #       push!(result, gw_time[Node.cs_detection_channel[i]][j][3])
    #       if Node.cs_detection_type[i][1] && Node.cs_detection_type[i][2]
    #         mis_detection += 1
    #       end
    #       break
    #     end
    #     # est_interference += 1
    #     if debug && length(Node.cs_detection_node[i]) == 1 && Node.cs_detection_node[i][1] < 10000 && Node_all[Node.cs_detection_node[i][1]].last_packet_success == 1
    #       # if debug && length(Node.cs_detection_node[i]) == 1 && Node.cs_detection_node[i][1] < 10000
    #       # println(gw_time)
    #       # println(Node.cs_detection_time[i])
    #       # println(Node.cs_detection_channel[i])
    #       # println(Node.cs_detection_node[i])
    #       # sqrt(-1)
    #       # print(Node_all[Node.cs_detection_node[i][1]].last_packet_success)
    #       # print("\t")
    #       # println(Node.cs_detection_node[i])
    #     end
    #     # for k = length(gw_collision_time):-1:1
    #     #   if gw_collision_time[k][1] <= Node.cs_detection_time[i]
    #     #     if gw_collision_time[k][2] >= Node.cs_detection_time[i]
    #     #       collision_count += 1
    #     #     end
    #     #     break
    #     #   end
    #     # end
    #     est_interference += 1
    #     if Node.cs_detection_type[i][1] && Node.cs_detection_type[i][2]
    #       both_interference += 1
    #     elseif Node.cs_detection_type[i][1] && !Node.cs_detection_type[i][2]
    #       only_insys_interference += 1
    #     elseif !Node.cs_detection_type[i][1] && Node.cs_detection_type[i][2]
    #       only_exsys_interference += 1
    #     else
    #       println("error of estimation")
    #     end
    #     break
    #   end
    # end
    # est_interference += 1
    # # if debug && length(Node.cs_detection_node[i]) == 1 && Node.cs_detection_node[i][1] < 10000
    # if debug && Node.cs_detection_node[i][1] < 10000
    #   print(Node_all[Node.cs_detection_node[i][1]].packet_collision_count)
    #   print("\t")
    #   println(Node.cs_detection_node[i])
    # end
    # if Node.cs_detection_type[i][1] && Node.cs_detection_type[i][2]
    #   both_interference += 1
    # elseif Node.cs_detection_type[i][1] && !Node.cs_detection_type[i][2]
    #   only_insys_interference += 1
    # elseif !Node.cs_detection_type[i][1] && Node.cs_detection_type[i][2]
    #   only_exsys_interference += 1
    # else
    #   println("error of estimation")
    # end
  end

  if isempty(result)
    return 0, est_interference, only_insys_interference, only_exsys_interference, both_interference, collision_count, mis_detection, backoff_channel_set, double_count
  end


  # if Node.last_subframe_number - 1 - length(result) != 0
  #   println(Node.last_subframe_number - 1 - length(result))
  # end
  return result, est_interference, only_insys_interference, only_exsys_interference, both_interference, collision_count, mis_detection, backoff_channel_set, double_count
end


#スレッド平均化
function sum_results(results_array::Vector{Array{Result}})
  num_arrays = length(results_array)
  num_results = length(results_array[1])
  results = Vector{Result}(undef, num_results)
  # println(results_array)

  for i in 1:num_results
    packet_all = 0
    slot_num = 0
    collision_num = 0
    packet_lost_num = 0
    same_id_count = 0
    same_id_success_count = 0
    difference_id_count = 0
    difference_id_success_count = 0
    snr_lost_num = 0

    for j in 1:num_arrays
      packet_all += results_array[j][i].packet_all
      slot_num += results_array[j][i].slot_num
      collision_num += results_array[j][i].collision_num
      packet_lost_num += results_array[j][i].packet_lost_num
      same_id_count += results_array[j][i].same_id_count
      same_id_success_count += results_array[j][i].same_id_success_count
      difference_id_count += results_array[j][i].difference_id_count
      difference_id_success_count += results_array[j][i].difference_id_success_count
      snr_lost_num += results_array[j][i].snr_lost_num
    end

    results[i] = Result(packet_all, slot_num, collision_num,
      packet_lost_num, same_id_count, same_id_success_count,
      difference_id_count, difference_id_success_count, snr_lost_num)
  end

  return results
end


#結果出力
function result_csv(result_array::Array, axis_name::String, horizontal_axis_unit1::String, horizontal_axis_unit2::String, csv_name::String)
  mkpath(string(csv_name, horizontal_axis_unit1))
  mkpath(string(csv_name, horizontal_axis_unit2))
  Threads.@threads for i = eachindex(result_array)
    result_horizontal_axis = []
    result_vertical_axis = []
    for j = eachindex(result_array[i][1])
      push!(result_horizontal_axis, j)
      push!(result_vertical_axis, result_array[i][1][j])
    end
    csv_name_no = result_array[i][2]
    df_result = DataFrame(horizontal_axis=result_horizontal_axis, vertical_axis=result_vertical_axis)
    rename!(df_result, Symbol.([string(axis_name, "_$horizontal_axis_unit1"), string(axis_name, "_", csv_name_no)]))
    df_result |> CSV.write(string(csv_name, horizontal_axis_unit1, "/$axis_name", "_$csv_name_no.csv"))
  end
  Threads.@threads for k = eachindex(result_array[1][1])
    result_horizontal_axis = []
    result_vertical_axis = []
    for l = eachindex(result_array)
      push!(result_horizontal_axis, result_array[l][2])
      push!(result_vertical_axis, result_array[l][1][k])
    end
    csv_name_no = k
    df_result = DataFrame(horizontal_axis=result_horizontal_axis, vertical_axis=result_vertical_axis)
    rename!(df_result, Symbol.([string(axis_name, "_$horizontal_axis_unit2"), string(axis_name, csv_name_no)]))
    df_result |> CSV.write(string(csv_name, horizontal_axis_unit2, "/$axis_name", "_$csv_name_no.csv"))
  end
  println("結果出力")
end

function write_csv(data_type::Int64, result_array::Array, device_num, plot_type::Int64, sim_num::Int64, clustering_time::Int64, cdf_calc_start_time::Vector{Float64}, cdf_calc_end_time::Vector{Float64})
  data1 = []
  data2 = []
  result = []
  result2 = []
  result3 = []
  type::String = "nothing"
  if data_type == 1
    # println(getfield(result_array[i], :success))
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
  elseif data_type == 5
    push!(data1, getfield.(result_array, :pdr_cdf))
    data2 = hcat(data1...)
    for j in 1:size(hcat(data1...), 1)
      push!(result, vcat(data2[j, :]...))
    end
    result = hcat(result...)
    type = "PDR_CDF"
  elseif data_type == 6
    for j in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :success))./sim_num)[Int(ceil(cdf_calc_start_time[j] / packet_period)):Int(ceil(cdf_calc_end_time[j] / packet_period))]))
      push!(data2, sum((sum(getfield.(result_array, :packet_num))./sim_num)[Int(ceil(cdf_calc_start_time[j] / packet_period)):Int(ceil(cdf_calc_end_time[j] / packet_period))]))
    end
    result = [[data1[i][j] / data2[i][j] for j in eachindex(data1[i])] for i in eachindex(data1)]
    type = "PDR-device"
  elseif data_type == 7
    push!(data1, getfield.(result_array, :pdr_cdf1))
    data2 = hcat(data1...)
    for j in 1:size(hcat(data1...), 1)
      push!(result, vcat(data2[j, :]...))
    end
    result = hcat(result...)
    type = "PDR_CDF1"
  elseif data_type == 8
    push!(data1, getfield.(result_array, :pdr_cdf2))
    data2 = hcat(data1...)
    for j in 1:size(hcat(data1...), 1)
      push!(result, vcat(data2[j, :]...))
    end
    result = hcat(result...)
    type = "PDR_CDF2"
  elseif data_type == 9
    push!(data1, getfield.(result_array, :pdr_cdf3))
    data2 = hcat(data1...)
    for j in 1:size(hcat(data1...), 1)
      push!(result, vcat(data2[j, :]...))
    end
    result = hcat(result...)
    type = "PDR_CDF3"
  elseif data_type == 10
    push!(data1, getfield.(result_array, :pdr_cdf4))
    data2 = hcat(data1...)
    for j in 1:size(hcat(data1...), 1)
      push!(result, vcat(data2[j, :]...))
    end
    result = hcat(result...)
    type = "PDR_CDF4"
  elseif data_type == 11
    for j in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :throughput))./sim_num)[Int(ceil(cdf_calc_start_time[j] / packet_period)):Int(ceil(cdf_calc_end_time[j] / packet_period))]) ./ (60.0 * (cdf_calc_end_time[j] - cdf_calc_start_time[j])))
    end
    result = data1
    type = "Throughput_device"
  elseif data_type == 12
    for j in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :throughput))./sim_num)[Int(ceil(cdf_calc_start_time[j] / packet_period)):Int(ceil(cdf_calc_end_time[j] / packet_period))]) ./ (60.0 * (cdf_calc_end_time[j] - cdf_calc_start_time[j]) * device_num))
    end
    result = data1
    type = "Throughput_device_ave"
  elseif data_type == 13
    push!(data1, getfield.(result_array, :cdf_throu))
    data2 = hcat(data1...)
    for j in 1:size(hcat(data1...), 1)
      push!(result, vcat(data2[j, :]...))
    end
    for k in axes(result, 1)
      for l in eachindex(cdf_calc_start_time)
        result[k][l] /= (60.0 * (cdf_calc_end_time[l] - cdf_calc_start_time[l]))
      end
    end
    result = hcat(result...)
    type = "throu_CDF"
  elseif data_type == 14
    for l in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :only_insys_interf))./sim_num)[Int(ceil(cdf_calc_start_time[l] / packet_period)):Int(ceil(cdf_calc_end_time[l] / packet_period))]))
      push!(data2, sum((sum(getfield.(result_array, :est_interf))./sim_num)[Int(ceil(cdf_calc_start_time[l] / packet_period)):Int(ceil(cdf_calc_end_time[l] / packet_period))]))
    end
    result = [[data1[i][j] / data2[i][j] for j in eachindex(data1[i])] for i in eachindex(data1)]
    type = "only_insys_interf"
  elseif data_type == 15
    for m in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :only_exsys_interf))./sim_num)[Int(ceil(cdf_calc_start_time[m] / packet_period)):Int(ceil(cdf_calc_end_time[m] / packet_period))]))
      push!(data2, sum((sum(getfield.(result_array, :est_interf))./sim_num)[Int(ceil(cdf_calc_start_time[m] / packet_period)):Int(ceil(cdf_calc_end_time[m] / packet_period))]))
    end
    result = [[data1[i][j] / data2[i][j] for j in eachindex(data1[i])] for i in eachindex(data1)]
    type = "only_exsys_interf"
  elseif data_type == 16
    for n in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :both_interf))./sim_num)[Int(ceil(cdf_calc_start_time[n] / packet_period)):Int(ceil(cdf_calc_end_time[n] / packet_period))]))
      push!(data2, sum((sum(getfield.(result_array, :est_interf))./sim_num)[Int(ceil(cdf_calc_start_time[n] / packet_period)):Int(ceil(cdf_calc_end_time[n] / packet_period))]))
    end
    result = [[data1[i][j] / data2[i][j] for j in eachindex(data1[i])] for i in eachindex(data1)]
    type = "both_interf"
  elseif data_type == 17
    for n in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :gw_collision_interf))./sim_num)[Int(ceil(cdf_calc_start_time[n] / packet_period)):Int(ceil(cdf_calc_end_time[n] / packet_period))]))
      push!(data2, sum((sum(getfield.(result_array, :est_interf))./sim_num)[Int(ceil(cdf_calc_start_time[n] / packet_period)):Int(ceil(cdf_calc_end_time[n] / packet_period))]))
    end
    result = [[data1[i][j] / data2[i][j] for j in eachindex(data1[i])] for i in eachindex(data1)]
    type = "gw_collision_interf"
  elseif data_type == 18
    for n in eachindex(cdf_calc_start_time)
      push!(data1, sum((sum(getfield.(result_array, :mis_detection))./sim_num)[Int(ceil(cdf_calc_start_time[n] / packet_period)):Int(ceil(cdf_calc_end_time[n] / packet_period))]))
      push!(data2, sum((sum(getfield.(result_array, :est_interf))./sim_num)[Int(ceil(cdf_calc_start_time[n] / packet_period)):Int(ceil(cdf_calc_end_time[n] / packet_period))]))
    end
    result = [[data1[i][j] / data2[i][j] for j in eachindex(data1[i])] for i in eachindex(data1)]
    type = "mis_detection"
  elseif data_type == 19
    push!(data1, sum(getfield.(result_array, :est_true_sum)) ./ sim_num)
    push!(data2, sum(getfield.(result_array, :est_sum)) ./ sim_num)
    result = [data1[i] ./ data2[i] for i in eachindex(data1)]
    type = "esit_hit"
  elseif data_type == 20
    push!(data1, sum(getfield.(result_array, :backoff_n)))
    result = data1
    type = "backoff"
  elseif data_type == 21
    push!(data1, sum(getfield.(result_array, :collision_n)))
    result = data1
    type = "collision"
  elseif data_type == 22
    push!(data1, sum(getfield.(result_array, :collision_n) ./ sum(getfield.(result_array, :backoff_n))))
    result = data1
    type = "hit_per"
  elseif data_type == 23
    push!(data1, sum(getfield.(result_array, :mis_n)))
    result = data1
    type = "mis"
  elseif data_type == 24
    push!(data1, sum(getfield.(result_array, :mis_n) ./ sum(getfield.(result_array, :collision_n))))
    result = data1
    type = "mis_per"
  else
    return 0
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
    df |> CSV.write("$output_folder/$clustering_time/$type-$device_num.csv", delim=',', writeheader=false)
  elseif plot_type == 2
    for i in axes(result, 1)
      # データをソートしてCDFを計算
      sort_data = vcat(result[i, :]...)
      sorted_data = sort(sort_data)
      cdf_values = cumsum(ones(length(sorted_data))) / length(sorted_data)
      # pushfirst!(sorted_data, 0.0)
      # push!(sorted_data, 1.0)
      # pushfirst!(cdf_values, 0.0)
      # push!(cdf_values, 1.0)

      # df = DataFrame(x=sorted_data, y=cdf_values)
      # CSVファイルに書き出し
      # df |> CSV.write(joinpath(output_folder, string(type, "_$i-$device_num.csv")), delim=',', writeheader=false)

      # Define bin edges from 0 to 1 with a step of 0.01
      bin_edges = 0:0.001:1
      hist = fit(Histogram, sorted_data, bin_edges)

      # Calculate bin centers
      bin_centers = [mean([bin_edges[i], bin_edges[i+1]]) for i in 1:length(bin_edges)-1]

      # Calculate cumulative relative frequencies
      total_count = sum(hist.weights)
      cumulative_frequencies = cumsum(hist.weights) / total_count

      # Create a DataFrame with bin centers and cumulative relative frequencies
      df = DataFrame(BinCenter=bin_centers, CumulativeRelativeFrequency=cumulative_frequencies)

      # CSVファイルに書き出し
      df |> CSV.write(joinpath(output_folder, string(clustering_time), string(type, "_$i-$device_num-classifying.csv")), delim=',', writeheader=false)
    end
  elseif plot_type == 3
    for i in eachindex(result)
      if !isfile(joinpath(output_folder, string(clustering_time), string(type, "_$i.csv")))
        df = DataFrame(x=device_num, y=result[i])
        df |> CSV.write(joinpath(output_folder, string(clustering_time), string(type, "_$i.csv")), delim=',', writeheader=false)
      else
        df = DataFrame(x=device_num, y=result[i])
        df |> CSV.write(joinpath(output_folder, string(clustering_time), string(type, "_$i.csv")), delim=',', writeheader=false, append=true)
      end
    end
  end
end

#結果出力グルーピング
function result_grouprate_csv(result_array::Array, axis_name::String, horizontal_axis_unit1::String, horizontal_axis_unit2::String, csv_name::String)
  mkpath(string(csv_name, horizontal_axis_unit1))
  mkpath(string(csv_name, horizontal_axis_unit2))
  Threads.@threads for i = eachindex(result_array)
    result_horizontal_axis = []
    result_vertical_axis = []
    for j = eachindex(result_array[i][1])
      push!(result_horizontal_axis, j)
      push!(result_vertical_axis, result_array[i][1][j])
    end
    csv_name_no = result_array[i][2]
    df_result = DataFrame(horizontal_axis=result_horizontal_axis, vertical_axis=result_vertical_axis)
    rename!(df_result, Symbol.([string(axis_name, "_$horizontal_axis_unit1"), string(axis_name, "_", csv_name_no)]))
    df_result |> CSV.write(string(csv_name, horizontal_axis_unit1, "/$axis_name", "_$csv_name_no.csv"))
  end
  Threads.@threads for k = eachindex(result_array[end][1])
    result_horizontal_axis = []
    result_vertical_axis = []
    for l = eachindex(result_array)
      push!(result_horizontal_axis, result_array[l][2])
      push!(result_vertical_axis, result_array[l][1][end])
    end
    csv_name_no = k
    df_result = DataFrame(horizontal_axis=result_horizontal_axis, vertical_axis=result_vertical_axis)
    rename!(df_result, Symbol.([string(axis_name, "_$horizontal_axis_unit2"), string(axis_name, csv_name_no)]))
    df_result |> CSV.write(string(csv_name, horizontal_axis_unit2, "/$axis_name", "_end_$csv_name_no.csv"))
  end
end

function result_cdf_csv(result_array::Array, axis_name::String, horizontal_axis_unit::String, csv_name::String)
  mkpath(string(csv_name, horizontal_axis_unit))
  Threads.@threads for i = eachindex(result_array)
    result_horizontal_axis = [n for n in 0:0.001:1]
    result_vertical_axis = []
    for j = eachindex(result_array[i][1])
      push!(result_vertical_axis, result_array[i][1][j])
    end
    csv_name_no = result_array[i][2]
    df_result = DataFrame(horizontal_axis=result_horizontal_axis, vertical_axis=result_vertical_axis)
    rename!(df_result, Symbol.([string(axis_name, "_$horizontal_axis_unit"), string(axis_name, "_", csv_name_no)]))
    df_result |> CSV.write(string(csv_name, horizontal_axis_unit, "/$axis_name", "_$csv_name_no.csv"))
  end
end


function export_grouped_csv(x::Vector, y::Vector, groups::Vector, nowtime::Int64, simindex::Int64, group_num::Int64)
  # グループごとにデータを分割
  grouped_data = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()  # グループごとのデータを格納する辞書
  for i in eachindex(groups)
    group = groups[i]
    x_coord = x[i] * 1000
    y_coord = y[i] * 1000
    if haskey(grouped_data, group)
      push!(grouped_data[group][1], x_coord)
      push!(grouped_data[group][2], y_coord)
    else
      grouped_data[group] = ([x_coord], [y_coord])
    end
  end

  # グループごとにCSVファイルを出力
  for (group, data) in grouped_data
    x_coords, y_coords = data
    # group_csv = CSV.File(joinpath(output_dir, "group_$group.csv"), header=["x", "y"])
    # append!(group_csv, DataFrame(x=x_coords, y=y_coords))
    # CSV.write(joinpath(output_dir, "group-$simindex-$nowtime-$group_num-$group.csv"), DataFrame(x=x_coords, y=y_coords))
    CSV.write("group-$simindex-$nowtime-$group_num-$group.csv", DataFrame(x=x_coords, y=y_coords))
  end
end
function replace_zeros_with_neg_inf(arr)
  for i in eachindex(arr)
    for j in eachindex(arr[1])
      if arr[i][j] == 0.0
        arr[i][j] = -Inf
      end
    end
  end
  return arr
end

function find_zero_rows_and_columns(matrix)
  num_rows, num_columns = size(matrix)
  zero_rows = Int[]
  zero_columns = Int[]

  for i in 1:num_rows
    if all(matrix[i, :] .== 0)
      push!(zero_rows, i)
    end
  end

  for j in 1:num_columns
    if all(matrix[:, j] .== 0)
      push!(zero_columns, j)
    end
  end

  return zero_rows, zero_columns
end

function get_edge_weight(graph, node1, node2)
  if graph.has_edge(node1, node2)
    return graph.get_edge_data(node1, node2)["weight"]
  else
    return 0.0  # エッジが存在しない場合は0を返す
  end
end

# function SC(Node_all, neighborhood_array, number_clusters)
#   similarity_matrix = exp.(-SC_hyper_pram .* neighborhood_array)
#   sc = cluster.spectral_clustering(n_clusters=number_clusters, affinity="precomputed", assign_labels="cluster_qr")
#   C = sc.fit(similarity_matrix)
#   plot = scatter(getfield.(Node_all, :x), getfield.(Node_all, :y), marker_z=C.labels_, color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
#   savefig(plot, string("cluster-1-", number_clusters, ".svg"))
# end

function cluster_qr(vectors)
  k = size(vectors, 1)
  q, r, piv = qr(vectors, Val(true))  # QR分解
  ut, s, v = svd(vectors[:, piv[1:k]])  # 特異値分解
  vectors = abs.(transpose(vectors) * (ut * transpose(conj(v)))) # 行列演算
  result = [idx[2] for idx in argmax(vectors, dims=2)[:, 1]]
  return result, vectors
end

# function cluster_qr(vectors)
#   k = size(vectors, 1)
#   q, r, piv = qr(vectors, Val(true))  # QR分解
#   ut, s, v = svd(vectors[:, piv[1:k]])  # 特異値分解
#   vectors = abs.(transpose(vectors) * (ut * transpose(conj(v)))) # 行列演算
#   result = [idx[2] for idx in argmax(vectors, dims=2)[:, 1]]
#   # 各要素の個数をカウント
#   counts = [count(x -> x == i, result) for i in 1:k]
#   vec_var = [var(vectors[i, :]) for i in eachindex(result)]
#   p = sortperm(vec_var)
#   for i in eachindex(p)
#     if (counts[sortperm(vectors[p[i], :], rev=true)[2]] < length(result) / k) && (counts[argmax(vectors[p[i], :])] > length(result) / k)
#       result[p[i]] = sortperm(vectors[p[i], :], rev=true)[2]
#       counts[sortperm(vectors[p[i], :], rev=true)[2]] += 1
#       counts[argmax(vectors[p[i], :])] -= 1
#     end
#     if maximum(counts) - minimum(counts) <= 1
#       break
#     end
#   end
#   return result, vectors
# end

# function cluster_qr_test(vectors)
#   k = size(vectors, 1)
#   q, r, piv = qr(vectors, Val(true))  # QR分解
#   ut, s, v = svd(vectors[:, piv[1:k]])  # 特異値分解
#   vectors = abs.(transpose(vectors) * (ut * transpose(conj(v)))) # 行列演算
#   result = [idx[2] for idx in argmax(vectors, dims=2)[:, 1]]
#   return result
# end

# function SC(Node_all, neighborhood_array, number_clusters, name::String, sim_num::Int64)
#   result = cluster.spectral_clustering(neighborhood_array, n_clusters=number_clusters, assign_labels="cluster_qr")
#   plot_sc = scatter(getfield.(Node_all, :x), getfield.(Node_all, :y), marker_z=result, color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
#   filename = joinpath(output_folder, string("cluster-", name, "-", number_clusters, "-$sim_num", ".svg"))
#   savefig(plot_sc, filename)
#   return result
# end

function SC_self(Node_all, neighborhood_array, number_clusters, name::String, sim_num::Int64, day::Int64, fig::Bool)
  num_device = length(Node_all)

  #アウテージ対策
  indices = findall(x -> x == 0, vec(sum(neighborhood_array, dims=1)))
  calc_node = setdiff(1:length(Node_all), indices)
  no_outage_array = neighborhood_array[Not(indices), Not(indices)]

  sparse(no_outage_array)

  # D^{-1/2} 行列を計算
  D = spdiagm(vec(sum(no_outage_array, dims=1)))
  dd = Array(sqrt.(diag(D)))
  # println(in(0.0, dd))
  D_sqrt_inv = spdiagm(vec(1 ./ sqrt.(sum(no_outage_array, dims=1))))

  L = D - sparse(no_outage_array)

  L_norm = D_sqrt_inv * L * D_sqrt_inv

  v0 = 2 * rand(size(no_outage_array, 1)) .- 1

  # ラプラシアン行列の最小固有値に対応する最小固有ベクトルを計算
  eigvals, eigvecs = Arpack.eigs(-L_norm, nev=number_clusters, sigma=1.0, v0=v0)
  smallest_eigenvector = eigvecs[:, 1:number_clusters]
  smallest_eigenvector ./= dd

  # k-meansを使用してクラスタリングを実行
  # result = kmeans(transpose(real.(smallest_eigenvector)), number_clusters)

  result, vectors = cluster_qr(transpose(real.(smallest_eigenvector)))
  # df_vec = DataFrame(hcat(hcat(getfield.(Node_all, :x)[Not(indices)], getfield.(Node_all, :y)[Not(indices)]), vectors), :auto)
  # df_vec |> CSV.write("vectors-$sim_num-$num_device.csv", delim=',', writeheader=false)
  # df_vec |> CSV.write(joinpath(output_folder, "vectors-$sim_num-$num_device.csv"), delim=',', writeheader=false)
  if sim_num == 1 && fig
    # if fig
    # plot_sc_self = scatter(getfield.(Node_all, :x)[Not(indices)], getfield.(Node_all, :y)[Not(indices)], marker_z=result, color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
    # plot_sc_self = scatter!(getfield.(Node_all, :x)[Not(calc_node)], getfield.(Node_all, :y)[Not(calc_node)], color="Gray", legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
    # filename = joinpath(output_folder, string("cluster-", name, "-", number_clusters, "-self", "-$sim_num", "-$num_device", "-$day", ".svg"))
    # savefig(plot_sc_self, string("cluster-", name, "-", number_clusters, "-self", "-$sim_num", "-$num_device", "-$day", ".svg"))
  end
  # for i in 1:number_clusters
  #   plot_vec = scatter(getfield.(Node_all, :x)[Not(indices)], getfield.(Node_all, :y)[Not(indices)], marker_z=vectors[:, i], clims=(0.0, 0.005), color=:bluesreds, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
  #   filename = joinpath(output_folder, string("cluster-", name, "-", number_clusters, "-self", "-$sim_num", "-$num_device", "-$day", "-$i", ".svg"))
  #   savefig(plot_vec, filename)
  # end
  return result, calc_node
end

function SC_self_sequentially_2(Node_all, neighborhood_array, same_channel_node, number_clusters, name::String, sim_num::Int64, day::Int64, fig::Bool)
  result_node = zeros(1:length(Node_all))
  result_calc_node = collect(1:length(Node_all))
  now_cluster = getfield.(Node_all, :group_id)
  neighborhood_array_norm = neighborhood_array ./ (same_channel_node ./ minimum(filter(x -> x != 0, same_channel_node)))
  for i in axes(neighborhood_array_norm, 1)
    neighborhood_array_norm[i, i] = 0
  end
  result = []
  for channel in 1:number_clusters
    exclusion_index = findall(x -> x == channel, now_cluster)
    push!(result, exclusion_index)
  end
  println(result)
  cluster_temp1, cluster_temp2, cluster_number = cluster_coupling(result, neighborhood_array_norm, [])
  println(cluster_temp1)
  println()
  println(cluster_temp2)
  println()
  println(cluster_number)
  sqrt(-1)
end

function SC_self_sequentially(Node_all, neighborhood_array, same_channel_node, number_clusters, name::String, sim_num::Int64, day::Int64, fig::Bool)
  result_node = zeros(1:length(Node_all))
  result_calc_node = collect(1:length(Node_all))
  now_cluster = getfield.(Node_all, :group_id)
  neighborhood_array_norm = neighborhood_array ./ (same_channel_node ./ minimum(filter(x -> x != 0, same_channel_node)))
  for i in axes(neighborhood_array_norm, 1)
    neighborhood_array_norm[i, i] = 0
  end
  # println(same_channel_node ./ minimum(filter(x -> x != 0, same_channel_node)))
  # println(neighborhood_array)
  # println(neighborhood_array_norm)
  result = []
  for channel in 1:number_clusters
    exclusion_index = findall(x -> x != channel, now_cluster)
    cluster, calc_node = SC_self(Node_all[Not(exclusion_index)], neighborhood_array_norm[Not(exclusion_index), Not(exclusion_index)], 2, name, sim_num, day, false)
    push!(result, (findall(x -> x == channel, now_cluster)[calc_node[findall(x -> x == 1, cluster)]], findall(x -> x == channel, now_cluster)[calc_node[findall(x -> x == 2, cluster)]]))
  end
  #二つに分けたクラスタ配列の平坦化
  cluster = collect(Iterators.flatten(result))
  node_cluster_index = zeros(length(Node_all))
  for node_cluster in eachindex(cluster)
    for i in eachindex(cluster[node_cluster])
      node_cluster_index[cluster[node_cluster][i]] = node_cluster
    end
  end
  if sim_num == 1 && fig
    # plot_sc_self = scatter(getfield.(Node_all, :x), getfield.(Node_all, :y), marker_z=node_cluster_index, color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
    # filename = joinpath(output_folder, string("cluster-", name, "-", number_clusters, "-double", "-self", "-$sim_num", "-$day", ".svg"))
    # savefig(plot_sc_self, filename)
  end

  if number_clusters == 4
    no_cluster1 = []
    no_cluster2 = []
    cluster_temp1, cluster_temp2, cluster_number = cluster_coupling(cluster, neighborhood_array_norm, [1, 10, 15])
    for i in 1:(size(collect(combinations(1:8, 4))[cluster_number], 1)-1)
      if (collect(combinations(1:8, 4))[cluster_number])[i+1] - (collect(combinations(1:8, 4))[cluster_number])[i] == 1
        if isodd((collect(combinations(1:8, 4))[cluster_number])[i])
          if i == 1 || i == 3
            push!(no_cluster1, 1)
          elseif i == 2
            push!(no_cluster1, 3)
          end
        end
      end
    end
    for ii in 1:(size(setdiff(collect(1:8), collect(combinations(1:8, 4))[cluster_number]), 1)-1)
      if (setdiff(collect(1:8), collect(combinations(1:8, 4))[cluster_number]))[ii+1] - (setdiff(collect(1:8), collect(combinations(1:8, 4))[cluster_number]))[ii] == 1
        if isodd((setdiff(collect(1:8), collect(combinations(1:8, 4))[cluster_number]))[ii])
          if ii == 1 || ii == 3
            push!(no_cluster2, 1)
          elseif ii == 2
            push!(no_cluster2, 3)
          end
        end
      end
    end
    cluster1, cluster2 = cluster_coupling(cluster_temp1, neighborhood_array_norm, [])
    cluster3, cluster4 = cluster_coupling(cluster_temp2, neighborhood_array_norm, [])
    for j in collect(Iterators.flatten(cluster1))
      result_node[j] = 1
    end
    for k in collect(Iterators.flatten(cluster2))
      result_node[k] = 2
    end
    for l in collect(Iterators.flatten(cluster3))
      result_node[l] = 3
    end
    for m in collect(Iterators.flatten(cluster4))
      result_node[m] = 4
    end

    for node_no_clustering in setdiff(setdiff(1:length(Node_all), vcat(collect(Iterators.flatten(result))...)), findall(x -> x == 0, vec(sum(neighborhood_array, dims=1))))
      result_node[node_no_clustering] = argmax([sum(neighborhood_array[node_no_clustering, collect(Iterators.flatten(i))]) for i in [cluster1, cluster2, cluster3, cluster4]])
    end

  else
    println("クラスタ数が予想されたものと異なります")
    println("クラスタ数:$number_clusters")
    sqrt(-1)
  end
  no_calc_node = findall(x -> x == 0, result_node)
  if sim_num == 1
    # plot_sc_self = scatter(getfield.(Node_all, :x)[Not(no_calc_node)], getfield.(Node_all, :y)[Not(no_calc_node)], marker_z=result_node[Not(no_calc_node)], color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
    # # plot_sc_self = scatter!(getfield.(Node_all, :x)[Not(calc_node)], getfield.(Node_all, :y)[Not(calc_node)], color="Gray", legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
    # filename = joinpath(output_folder, string("cluster-", name, "-", number_clusters, "-self", "-$sim_num", "-$day", ".svg"))
    # savefig(plot_sc_self, filename)
  end
  return result_node[Not(no_calc_node)], result_calc_node[Not(no_calc_node)]
end

function cluster_coupling(cluster, neighborhood_array_norm, no_cluster)
  node_comb = collect(combinations(cluster, (length(cluster) ÷ 2)))[collect(1:(binomial(length(cluster), (length(cluster) ÷ 2))÷2))]
  cluster_relation = [(sum(neighborhood_array_norm[collect(Iterators.flatten(i)), collect(Iterators.flatten(i))]) + sum(neighborhood_array_norm[setdiff(collect(Iterators.flatten(cluster)), collect(Iterators.flatten(i))), setdiff(collect(Iterators.flatten(cluster)), collect(Iterators.flatten(i)))])) for i in node_comb]
  if !isempty(no_cluster)
    for i in no_cluster
      cluster_relation[i] = 0
    end
  end
  return collect(combinations(cluster, (length(cluster) ÷ 2)))[argmax(cluster_relation)], setdiff(cluster, collect(combinations(cluster, (length(cluster) ÷ 2)))[argmax(cluster_relation)]), argmax(cluster_relation)
end

function calc_clockdrift(last_time::Float64, now_time::Float64, CD_mean::Float64, CD_variance::Float64)
  elapsed_time = Int(round(now_time - last_time))
  if elapsed_time <= 0
    # println("elapsed_timeは:$elapsed_time")
    return 0, 1
  end
  # println(elapsed_time)
  # return sum(rand(Normal(CD_mean, CD_variance), elapsed_time)), 0
  return sum(rand(Normal(CD_mean, CD_variance), elapsed_time)) / 60, 0
end

function write_png(result_array, num_device)
  result = []
  push!(result, getfield.(result_array, :node_num))
  result_flatten = collect(Iterators.flatten(Iterators.flatten(result)))
  p = histogram(result_flatten, label=false)
  filename = joinpath(output_folder, "cluster-$num_device.svg")
  savefig(p, filename)
end
function write_png2(result_array, num_device)
  result = []
  push!(result, getfield.(result_array, :node_pdf))
  result_flatten = collect(Iterators.flatten(Iterators.flatten(result)))
  p = histogram(result_flatten, label=false)
  filename = joinpath(output_folder, "pdfnode-$num_device.svg")
  savefig(p, filename)
end

function SC_self_nosparse(Node_all, neighborhood_array, number_clusters, name::String, sim_num::Int64)
  #アウテージ対策
  indices = findall(x -> x == 0, vec(sum(neighborhood_array, dims=1)))
  calc_node = setdiff(1:length(Node_all), indices)
  no_outage_array = neighborhood_array[Not(indices), Not(indices)]

  # D^{-1/2} 行列を計算
  D = diagm(vec(sum(no_outage_array, dims=1)))
  dd = sqrt.(diag(D))
  # println(in(0.0, dd))
  D_sqrt_inv = diagm(vec(1 ./ sqrt.(sum(no_outage_array, dims=1))))

  L = D - no_outage_array

  L_norm = D_sqrt_inv * L * D_sqrt_inv

  v0 = 2 * rand(size(no_outage_array, 1)) .- 1

  # ラプラシアン行列の最小固有値に対応する最小固有ベクトルを計算
  eigvals, eigvecs = Arpack.eigs(-L_norm, nev=number_clusters, sigma=1.0, v0=v0, maxiter=100000000)
  smallest_eigenvector = eigvecs[:, 1:number_clusters]
  smallest_eigenvector ./= dd

  # k-meansを使用してクラスタリングを実行
  # result = kmeans(transpose(real.(smallest_eigenvector)), number_clusters)

  result, vectors = cluster_qr(transpose(real.(smallest_eigenvector)))
  # plot_sc_self = scatter(getfield.(Node_all, :x)[Not(indices)], getfield.(Node_all, :y)[Not(indices)], marker_z=result, color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
  # plot_sc_self = scatter!(getfield.(Node_all, :x)[Not(calc_node)], getfield.(Node_all, :y)[Not(calc_node)], color="Gray", legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
  # filename = joinpath(output_folder, string("cluster-", name, "-", number_clusters, "-self", "-$sim_num", ".svg"))
  # savefig(plot_sc_self, filename)
  return result, calc_node
end

function write_shadowingmap(z::Array, x::Float64, y::Float64, Node_all::Array)
  plot_sh = scatter(getfield.(Node_all, :x), getfield.(Node_all, :y), marker_z=z, color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
  savefig(plot_sh, joinpath(output_folder, "shadowingmap-x=$x-y=$y.svg"))
end

function total_csv()
  # DataFrameを初期化
  combined_df = DataFrame()
  for i in clustering_time_size
    df = CSV.read(joinpath(output_folder_top, string(i), "PDR-device_final.csv"), DataFrame, delim=',', header=false, drop=[1])
    rename!(df, Symbol.([string(i)]))
    # 読み込んだDataFrameを結合
    combined_df = hcat(combined_df, df)
  end
  transposed_matrix = hcat(clustering_time_size, transpose(Matrix(combined_df)))
  transposed_df = DataFrame(transposed_matrix, :auto)
  for (index, value) in enumerate(SIM_particle_size:SIM_particle_size:number_devices_max)
    select(transposed_df, [1, index + 1]) |> CSV.write(joinpath(output_folder_top, string("PDR_time_numberdevice=$value.csv")), delim=',', writeheader=false)
  end
end
