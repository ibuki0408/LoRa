using Distributed
# デバッグモードの切り替え
const debug = false

ENV["GKSwstype"] = "100"
@everywhere using StatsBase
@everywhere using BenchmarkTools
@everywhere using CSV
@everywhere using DataFrames
@everywhere using ProgressMeter
@everywhere using Plots
@everywhere using Dates
@everywhere using Random
@everywhere using Distributions
@everywhere using Combinatorics
@everywhere using Statistics
@everywhere using LinearAlgebra
@everywhere import Base.Filesystem.mkdir
@everywhere import Arpack
@everywhere using SparseArrays
@everywhere using StatsPlots
@everywhere using InvertedIndices
@everywhere using FileIO
@everywhere using StaticArrays
@everywhere using DataStructures

@everywhere include("param.jl")
@everywhere include("Node.jl")
@everywhere include("result.jl")
@everywhere include("base_function.jl")
@everywhere include("shadowing.jl")
@everywhere include("ppp.jl")
@everywhere include("interference.jl")

#シミュレーション回数
@everywhere const SIM_NUM = 1000

#シミュレーション粒度
@everywhere const SIM_start_size = 250
@everywhere const SIM_particle_size = 250
@everywhere const number_devices_max = 2000 #端末台数
@everywhere const clustering_start_time = 36
@everywhere const cdf_calc_start_time = [10.0, clustering_start_time * 60.0 + 30.0]
@everywhere const cdf_calc_end_time = [clustering_start_time * 60.0 - 30.0, SIM_PERIOD - 60.0]

#現在時刻
@everywhere const date = string(Dates.format(now(UTC) + Hour(9), "yyyymmdd_HHMMSS"), "_$SIM_NUM", "_$carrier_sense_threshold")

if !isdir("result")
  mkdir("result")
end

# 出力フォルダを指定（現在のディレクトリ内の"result"フォルダに保存）
if !debug
  @everywhere output_folder = joinpath("result", date)
  # 出力フォルダが存在しない場合は作成
  if !isdir(output_folder)
    mkdir(output_folder)
  end
  cp("RCS_test.jl", joinpath(output_folder, "RCS_test.jl"))
  cp("param.jl", joinpath(output_folder, "param.jl"))
  cp("base_function.jl", joinpath(output_folder, "base_function.jl"))
  cp("Node.jl", joinpath(output_folder, "Node.jl"))
  cp("result.jl", joinpath(output_folder, "result.jl"))
end

@everywhere function main_roop(sim_index::Int64, number_devices::Int64)
  #結果書き出し
  collision_count = 0
  backoff_count = 0
  mis_count = 0
  double_count = 0
  collision_packet = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  success_packet = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  rost_packet = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  packet_num = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  throughput = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  estimation_hit = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  estimation_num = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  counter_est_interf = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  counter_only_insys_interf = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  counter_only_exsys_interf = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  counter_both_interf = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  counter_gw_collision_interf = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  counter_mis_detection = zeros(Int(ceil(SIM_PERIOD / packet_period)))
  cdf = Array{Array{Float64}}(undef, 0)
  cdf_1 = Array{Array{Float64}}(undef, 0)
  cdf_2 = Array{Array{Float64}}(undef, 0)
  cdf_3 = Array{Array{Float64}}(undef, 0)
  cdf_4 = Array{Array{Float64}}(undef, 0)
  throughput_cdf = Array{Array{Float64}}(undef, 0)
  throughput_cdf_calc = zeros(number_devices)
  result_node = []
  result_pdf = []
  neighborhood_node_sum = 0
  neighborhood_node_true_sum = 0
  error = 0
  no_error = 0
  clustering_flag = 0
  #時刻
  now_time::Float64 = 0.0
  event_time = Array{Tuple{Float64,Int64,Int64}}(undef, 0) #(時刻，コマンド(1:CS開始，2:CS終了，3:パケット送信開始，4:パケット送信完了，5:受信窓でのCS, 6:クラスタリング, 7:CDF計算開始, 8:CDF計算終了)，EN番号)
  gw_time = [Vector{Tuple{Float64,Float64,Int64}}() for _ in 1:channel_num] #(送信開始時刻，送信終了時刻，端末番号)
  interference_time = Array{Tuple{Float64,Float64,Int64}}(undef, 0) #(送信開始時刻，送信終了時刻，干渉端末番号)
  gw_collision_time = Array{Tuple{Float64,Float64,Int64}}(undef, 0) #(送信開始時刻，送信終了時刻，端末番号)
  GW_obs_CD = [CircularBuffer{Float64}(120) for _ in 1:number_devices]
  GW_tmp_obs_CD = zeros(Float64, number_devices)
  #周波数チャネル
  now_channel::Int64 = 0
  #ノード確保配列
  Node_all = Array{N_parameters}(undef, number_devices)
  #干渉端末生成
  # interference_xy = poisson_point_process(0.0004, @SVector [area_size * 1000, area_size * 1000]) #エリア全体に干渉端末配置
  # interference_xy = [x .+ (area_size * 1000 / 2) for x in poisson_point_process(0.0016, @SVector [area_size * 1000 / 2, area_size * 1000 / 2])]  #エリア第一象限に干渉端末配置
  interference_xy = []
  Interfere_all = Array{interference}(undef, length(interference_xy))
  #ノード生成
  node_generator(Node_all, number_devices, 0, "uplink_toa")
  interference_generator(Interfere_all, interference_xy, "uplink_toa")
  # if sim_index == 1
  #   plot_sc_self = scatter(getfield.(Node_all, :x), getfield.(Node_all, :y), color=:blue, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
  #   plot_sc_self = scatter!(getfield.(Interfere_all, :x), getfield.(Interfere_all, :y), color="Gray", legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
  #   savefig(plot_sc_self, joinpath(output_folder, string("interf_node_posi-", number_devices, ".svg")))
  # end
  # sqrt(-1)
  #シャドウイング値生成
  waves = Array{F_T}(undef, SW_N)
  generate_waves(waves)

  #初期パケット生成
  for node_index = 1:number_devices
    Node_all[node_index].packet_num += 1
    Node_all[node_index].channel, Node_all[node_index].last_packet_slot, packet_generation_time = send_packet_generator(Node_all[node_index])
    Node_all[node_index].last_subframe_number += 1
    Node_all[node_index].packet_size = ((8 + 4.25 + 8) + ceil((160) / Node_all[node_index].sf)) * ((2^Node_all[node_index].sf) / BW) / 60
    #クロックドリフト計算
    CD, time_error = calc_clockdrift(Node_all[node_index].CD_last_calc * 60, packet_generation_time * 60, Node_all[node_index].clockdrift_mean, Node_all[node_index].clockdrift_variance)
    if time_error == 0
      no_error += 1
    else
      error += 1
    end
    GW_tmp_obs_CD[node_index] += CD
    Node_all[node_index].packet_generation_time = packet_generation_time + CD
    Node_all[node_index].CD_accumulated += CD
    Node_all[node_index].CD_last_calc = Node_all[node_index].packet_generation_time
    push!(event_time, (Node_all[node_index].packet_generation_time - 5 / 60000, 1, node_index))
    push!(event_time, (Node_all[node_index].packet_generation_time, 2, node_index))
  end
  for inf_index in eachindex(Interfere_all)
    push!(event_time, (Interfere_all[inf_index].offset - 5 / 60000, 1, Interfere_all[inf_index].node_id))
    push!(event_time, (Interfere_all[inf_index].offset, 2, Interfere_all[inf_index].node_id))
    Interfere_all[inf_index].packet_generation_time = Interfere_all[inf_index].offset
  end
  #SNR規範により各端末のSF決定
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
  #クラスタリング時間およびCDF計算時間の出力
  push!(event_time, (clustering_start_time * 60.0, 6, 0))
  for cdf_start_time in cdf_calc_start_time
    push!(event_time, (cdf_start_time, 7, 0))
  end
  for cdf_end_time in cdf_calc_end_time
    push!(event_time, (cdf_end_time, 8, 0))
  end
  #event_timeをソート
  sort!(event_time, by=x -> x[1])

  #各種カウント変数定義
  counter_gw::Int64 = 0
  counter_en::Int64 = 0
  counter_bk::Int64 = 0
  counter_backoff_channel_set = [zeros(Int, channel_num) for _ in 1:number_devices]
  #フラグ
  calc_cdf::Bool = false
  device_interf::Bool = false

  #隣接行列
  neighborhood_node = zeros(Int64, number_devices, number_devices)
  neighborhood_node_cs = zeros(Int64, number_devices, number_devices)
  for i in axes(neighborhood_node_cs, 1)
    for j in axes(neighborhood_node_cs, 1)
      if i == j
        neighborhood_node_cs[i, j] = 0
      else
        if recieve_power(i, [Node_all[i].x, Node_all[i].y], [Node_all[j].x, Node_all[j].y], Tx_dB, waves, 1, 0) >= carrier_sense_threshold
          neighborhood_node_cs[i, j] = 1
        else
          neighborhood_node_cs[i, j] = 0
        end
      end
    end
  end

  #シミュレイーション本体
  while !isempty(event_time)
    #event_timeをソート
    sort!(event_time, by=x -> x[1])
    node_id::Int64 = event_time[1][3]
    if node_id >= 10000
      device_interf = true
    else
      device_interf = false
    end
    if device_interf
      now_channel = Interfere_all[node_id%10000].channel
    elseif !device_interf && node_id == 0
      now_channel = 0
    else
      now_channel = Node_all[node_id].channel
    end
    now_time = event_time[1][1]
    # if sim_index == 1
    #   println(now_time)
    # end
    now_slot = Int(ceil(now_time / packet_period))
    if now_slot > Int(ceil(SIM_PERIOD / packet_period))
      break
    end
    send_success::Bool = false

    #イベント操作
    if event_time[1][2] == 1 ##キャリアセンス開始
      if device_interf
        Interfere_all[node_id%10000].status = 1
        if carrier_sense(node_id, Node_all, Interfere_all, waves) == 1
          Interfere_all[node_id%10000].cs_detection = 1
          push!(Interfere_all[node_id%10000].cs_detection_time, now_time)
        end
      else
        Node_all[node_id].status = 1
        if carrier_sense(node_id, Node_all, Interfere_all, waves) == 1
          Node_all[node_id].cs_detection = 1
          push!(Node_all[node_id].cs_detection_time, now_time)
        end
      end
      popfirst!(event_time)
      continue
    elseif event_time[1][2] == 2 ##キャリアセンス終了
      if device_interf
        Interfere_all[node_id%10000].status = 0
        if Interfere_all[node_id%10000].cs_detection == 0
          Interfere_all[node_id%10000].status = 2
          push!(event_time, (event_time[1][1] + Interfere_all[node_id%10000].packet_size, 4, node_id))
          for cs_update = 1:number_devices
            if Node_all[cs_update].status == 1 && Node_all[cs_update].cs_detection == 0
              if carrier_sense(cs_update, Node_all, Interfere_all, waves) == 1
                Node_all[cs_update].cs_detection = 1
                push!(Node_all[cs_update].cs_detection_time, now_time)
              end
            end
          end
          for cs_update in eachindex(Interfere_all)
            if Interfere_all[cs_update].status == 1 && Interfere_all[cs_update].cs_detection == 0
              if carrier_sense(cs_update, Node_all, Interfere_all, waves) == 1
                Interfere_all[cs_update].cs_detection = 1
                push!(Interfere_all[cs_update].cs_detection_time, now_time)
              end
            end
          end
          #衝突判定
          collision_judge(node_id, Node_all, Interfere_all, waves, device_interf)
          push!(interference_time, (event_time[1][1], event_time[1][1] + Interfere_all[node_id%10000].packet_size, node_id))
        else
          push!(event_time, (Interfere_all[node_id%10000].packet_generation_time + packet_period - 5 / 60000, 1, Interfere_all[node_id%10000].node_id))
          push!(event_time, (Interfere_all[node_id%10000].packet_generation_time + packet_period, 2, Interfere_all[node_id%10000].node_id))
          Interfere_all[node_id%10000].packet_generation_time += packet_period
        end
        #初期化
        Interfere_all[node_id%10000].cs_detection = 0
        popfirst!(event_time)
        continue
      else
        Node_all[node_id].status = 0
        if Node_all[node_id].cs_detection == 1 && Node_all[node_id].last_subframe_number == subframe_number
          #パケット破棄処理
          popfirst!(event_time)
          packet_num[now_slot] += 1
          rost_packet[now_slot] += 1
          counter_en += 1
        elseif Node_all[node_id].cs_detection == 1 && Node_all[node_id].last_subframe_number != subframe_number
          #バックオフ処理
          counter_bk += 1
          Node_all[node_id].channel, Node_all[node_id].last_packet_slot, packet_generation_time = backoff(Node_all[node_id], Node_all[node_id].last_subframe_number)
          Node_all[node_id].last_subframe_number += 1
          #クロックドリフト計算
          CD, time_error = calc_clockdrift(Node_all[node_id].CD_last_calc * 60, packet_generation_time * 60, Node_all[node_id].clockdrift_mean, Node_all[node_id].clockdrift_variance)
          if time_error == 0
            no_error += 1
          else
            error += 1
            # println("クロックドリフト計算エラーbackoff")
          end
          GW_tmp_obs_CD[node_id] += CD
          Node_all[node_id].packet_generation_time = packet_generation_time + CD
          Node_all[node_id].CD_accumulated += CD
          Node_all[node_id].CD_last_calc = Node_all[node_id].packet_generation_time
          push!(event_time, (Node_all[node_id].packet_generation_time - 5 / 60000, 1, node_id))
          push!(event_time, (Node_all[node_id].packet_generation_time, 2, node_id))
          popfirst!(event_time)
          #初期化
          Node_all[node_id].cs_detection = 0
          continue
        else
          #送信開始
          counter_gw += 1
          Node_all[node_id].status = 2
          push!(event_time, (event_time[1][1] + Node_all[node_id].packet_size, 4, node_id))
          # push!(event_time, (event_time[1][1] + Node_all[node_id].packet_size + 1 / 60, 5, node_id))
          for cs_update = 1:number_devices
            if Node_all[cs_update].status == 1 && Node_all[cs_update].cs_detection == 0
              if carrier_sense(cs_update, Node_all, Interfere_all, waves) == 1
                Node_all[cs_update].cs_detection = 1
                push!(Node_all[cs_update].cs_detection_time, now_time)
              end
            end
          end
          for cs_update in eachindex(Interfere_all)
            if Interfere_all[cs_update].status == 1 && Interfere_all[cs_update].cs_detection == 0
              if carrier_sense(cs_update, Node_all, Interfere_all, waves) == 1
                Interfere_all[cs_update].cs_detection = 1
                push!(Interfere_all[cs_update].cs_detection_time, now_time)
              end
            end
          end
          #初期化
          Node_all[node_id].cs_detection = 0
          popfirst!(event_time)
          #衝突判定
          collision_judge(node_id, Node_all, Interfere_all, waves, device_interf)
          continue
        end
      end
    elseif event_time[1][2] == 4
      if device_interf
        Interfere_all[node_id%10000].status = 0
        popfirst!(event_time)
        push!(event_time, (Interfere_all[node_id%10000].packet_generation_time + packet_period - 5 / 60000, 1, Interfere_all[node_id%10000].node_id))
        push!(event_time, (Interfere_all[node_id%10000].packet_generation_time + packet_period, 2, Interfere_all[node_id%10000].node_id))
        Interfere_all[node_id%10000].packet_generation_time += packet_period
        Interfere_all[node_id%10000].cs_detection = 0
        Interfere_all[node_id%10000].tx_collision = 0
        Interfere_all[node_id%10000].cs_detection_time = []
        continue
      else
        #送信終了
        Node_all[node_id].status = 0
        if SNR_judge(Node_all[node_id], waves) == 0 && Node_all[node_id].tx_collision != [true, true]
          if Node_all[node_id].tx_collision == [false, false]
            send_success = true
            packet_num[now_slot] += 1
            success_packet[now_slot] += 1
            if calc_cdf
              Node_all[node_id].packet_success_cdf += 1
              throughput_cdf_calc[node_id] += length(Node_all[node_id].send_bit) + payload_bit
            end
            throughput[now_slot] += length(Node_all[node_id].send_bit) + payload_bit
            push!(gw_time[now_channel], (Node_all[node_id].packet_generation_time, Node_all[node_id].packet_generation_time + Node_all[node_id].packet_size, node_id)) #五行下とどちらかをコメントアウト（プリアンブル衝突判断）
            push!(GW_obs_CD[node_id], GW_tmp_obs_CD[node_id])
            GW_tmp_obs_CD[node_id] = 0
          end
          # push!(gw_time[now_channel], (Node_all[node_id].packet_generation_time, Node_all[node_id].packet_generation_time + Node_all[node_id].packet_size, node_id))
          # filter!(x -> x[2] >= now_time - 4.0, gw_time)
          # sort!(gw_time, by=x -> x[2])
          if clustering_flag == 0
            set_result, est_interf, only_insys_interf, only_exsys_interf, both_interf, collision_interf, mis_detection, backoff_channel_set, double_hit = estimation(gw_time, Node_all[node_id], interference_time, Interfere_all, waves, Node_all, gw_collision_time, GW_obs_CD, now_time)
            counter_est_interf[now_slot] += est_interf
            counter_only_insys_interf[now_slot] += only_insys_interf
            counter_only_exsys_interf[now_slot] += only_exsys_interf
            counter_both_interf[now_slot] += both_interf
            counter_gw_collision_interf[now_slot] += collision_interf
            counter_mis_detection[now_slot] += mis_detection
            counter_backoff_channel_set[node_id] .+= backoff_channel_set

            backoff_count += length(Node_all[node_id].cs_detection_time)
            collision_count += collision_interf
            double_count += double_hit
            mis_count += mis_detection
            if set_result != 0
              for set in set_result
                # push!(estimation_node, [Node_all[node_id].packet_num, node_id, set])
                neighborhood_node[node_id, set] += 1
                neighborhood_node[set, node_id] += 1
                # push!(estimation_node, [node_id, set])
              end
            end
          end
          # elseif SNR_judge(Node_all[node_id]) == 1
          #   push!(sended_packets, (Node_all[node_id].packet_generation_time, Node_all[node_id].packet_num, 3, node_id))
        elseif SNR_judge(Node_all[node_id], waves) == 0 && Node_all[node_id].tx_collision == [true, true]
          # push!(sended_packets, (Node_all[node_id].packet_generation_time, Node_all[node_id].packet_num, 2, node_id))
          push!(gw_collision_time, (Node_all[node_id].packet_generation_time, now_time, node_id))
          packet_num[now_slot] += 1
          collision_packet[now_slot] += 1
        else
          packet_num[now_slot] += 1
          # println("SNR_Drop")
        end
        popfirst!(event_time)
      end
    elseif event_time[1][2] == 5
      #将来ダウンリンク追加予定
      popfirst!(event_time)
      continue
    elseif event_time[1][2] == 6

      # neighborhood_hist = histogram(vec(neighborhood_node), bins=nbins = maximum(vec(neighborhood_node)) - minimum(vec(neighborhood_node)))
      # savefig(neighborhood_hist, joinpath(output_folder, string("neighborhood_hist", "-$sim_index", ".svg")))

      # println("1")
      # # 閾値を行列の平均値として設定
      # # hist = fit(Histogram, Float64.(vec(neighborhood_node) .+ 1); nbins=maximum(Float64.(vec(neighborhood_node) .+ 1)) - minimum(Float64.(vec(neighborhood_node) .+ 1)))
      # println("2")
      # # threshold_neiborhood = find_threshold(hist.weights, (hist.edges[1])[1:end-1], Otsu()) - 1
      # threshold_neiborhood = mean(neighborhood_node)
      # println("3")
      # # 閾値に基づいて行列の要素を置き換える
      # # neighborhood_node_norm = map(x -> x >= threshold_neiborhood ? 1 : 0, neighborhood_node)
      # neighborhood_node_norm = map(x -> x >= 1 ? 1 : 0, neighborhood_node)
      # println("4")
      # # 行列 A と B の要素が異なるかどうかを比較
      # diff_neighborhood = neighborhood_node_norm .!= neighborhood_node_cs
      # println("5")
      # # 異なる要素の数を数える
      # num_diff_elements = sum(diff_neighborhood[5, :])
      # println("6")
      # # 行列の総要素数を計算
      # total_elements = sum(neighborhood_node_cs[5, :])
      # println("7")
      # # 一致する要素の割合を計算
      # hit_probability = 1 - (num_diff_elements / total_elements)
      # println("8")
      # # 異なる要素の割合を計算
      # println("隠れ端末推定率: $hit_probability")

      # if sim_index == 1
      #   plot_sc_self = scatter(getfield.(Node_all, :x), getfield.(Node_all, :y), marker_z=neighborhood_node_norm[5, :], color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
      #   plot_sc_self = scatter!([(getfield.(Node_all, :x))[5]], [(getfield.(Node_all, :y))[5]], color="Green", legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
      #   filename = joinpath(output_folder, string("est_hidden", "-$sim_index", ".svg"))
      #   savefig(plot_sc_self, filename)
      #   plot_sc_self1 = scatter(getfield.(Node_all, :x), getfield.(Node_all, :y), marker_z=neighborhood_node_cs[5, :], color=:jet, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
      #   plot_sc_self1 = scatter!([(getfield.(Node_all, :x))[5]], [(getfield.(Node_all, :y))[5]], color="Green", legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")
      #   filename = joinpath(output_folder, string("hidden", "-$sim_index", ".svg"))
      #   savefig(plot_sc_self1, filename)
      # end

      # CSVファイルに行列を書き出す
      # if !debug && sim_index == 1
      #   df = DataFrame(neighborhood_node, :auto)
      #   df |> CSV.write(joinpath(output_folder, "neighborhood_node$sim_index-$number_devices.csv"), delim=',', writeheader=false)
      # end
      neighborhood_node_true = neighborhood_node .* neighborhood_node_cs
      neighborhood_node_true_sum = sum(neighborhood_node_true)
      neighborhood_node_sum = sum(neighborhood_node)
      # println("CS可否推定一致率: $(neighborhood_node_true_sum / neighborhood_node_sum * 100) %")
      #クラスタリング
      # SC_self(Node_all, neighborhood_node, cluster_num, "$now_time", sim_index)
      clustering_flag = 1
      cluster, calc_node = SC_self(Node_all, neighborhood_node, cluster_num, "$now_time", sim_index, Int(now_time), true)
      cluster_backoff_channel = [zeros(Int, channel_num) for _ in 1:cluster_num]
      cluster_channel_allocation = zeros(Int, cluster_num)
      for i in cluster
        cluster_backoff_channel[cluster[i]][argmax(counter_backoff_channel_set[calc_node[i]])] += 1
      end
      matrix_cluster_backoff_channel = hcat(cluster_backoff_channel...)'
      for _ in 1:cluster_num
        if !(argmax(matrix_cluster_backoff_channel)[2] in cluster_channel_allocation)
          cluster_channel_allocation[argmax(matrix_cluster_backoff_channel)[1]] = argmax(matrix_cluster_backoff_channel)[2]
          matrix_cluster_backoff_channel[:, argmax(matrix_cluster_backoff_channel)[2]] .= 0
        end
      end
      existing_values = setdiff(Set(cluster_channel_allocation), 0)
      missing_values = setdiff(1:cluster_num, existing_values)
      shuffle!(missing_values)
      for i in 1:cluster_num
        if cluster_channel_allocation[i] == 0
          cluster_channel_allocation[i] = pop!(missing_values)
        end
      end
      result_node = collect(values(countmap(cluster)))
      # export_grouped_csv(getfield.(Node_all, :x), getfield.(Node_all, :y), cluster, Int(now_time / 60), sim_index, cluster_num)
      # cluster, calc_node = SC_self(Node_all, neighborhood_node, 8, "$now_time", sim_index, Int(now_time))
      # result_node = collect(values(countmap(cluster)))
      # export_grouped_csv(getfield.(Node_all, :x), getfield.(Node_all, :y), cluster, output_folder, Int(now_time / 60), sim_index, 8)
      for (index, value) in enumerate(calc_node)
        Node_all[value].usable_channel = [cluster[index]]
        Node_all[value].group_id = cluster[index]
        # Node_all[value].usable_channel = [cluster_channel_allocation[cluster[index]]]
        # Node_all[value].group_id = cluster[index]
        Node_all[value].send_bit = zeros(number_PLIM_bits(length(Node_all[value].usable_channel), slot_number))
      end
      popfirst!(event_time)
      continue
    elseif event_time[1][2] == 7
      calc_cdf = true
      for i = 1:number_devices
        Node_all[i].packet_num_cdf = 0
        Node_all[i].packet_success_cdf = 0
        throughput_cdf_calc[i] = 0
      end
      popfirst!(event_time)
      continue
    elseif event_time[1][2] == 8
      calc_cdf = false
      push!(cdf, getfield.(Node_all, :packet_success_cdf) ./ getfield.(Node_all, :packet_num_cdf))
      push!(throughput_cdf, copy(throughput_cdf_calc))
      push!(cdf_1, [getfield(Node, :packet_success_cdf) / getfield(Node, :packet_num_cdf) for Node in Node_all if getfield(Node, :usable_channel) == [1]])
      push!(cdf_2, [getfield(Node, :packet_success_cdf) / getfield(Node, :packet_num_cdf) for Node in Node_all if getfield(Node, :usable_channel) == [2]])
      push!(cdf_3, [getfield(Node, :packet_success_cdf) / getfield(Node, :packet_num_cdf) for Node in Node_all if getfield(Node, :usable_channel) == [3]])
      push!(cdf_4, [getfield(Node, :packet_success_cdf) / getfield(Node, :packet_num_cdf) for Node in Node_all if getfield(Node, :usable_channel) == [4]])
      result_pdf = [sum([getfield(Node, :packet_success_cdf) for Node in Node_all if getfield(Node, :usable_channel) == [i]]) / sum([getfield(Node, :packet_num_cdf) for Node in Node_all if getfield(Node, :usable_channel) == [i]]) for i in 1:channel_num]
      popfirst!(event_time)
      continue
    else
      println("予期しないイベント番号です")
      println(event_time[1][2])
      sqrt(-1)
    end
    #次のパケット生成
    Node_all[node_id].send_bit = rand(0:1, length(Node_all[node_id].send_bit))
    Node_all[node_id].packet_num += 1
    Node_all[node_id].channel, Node_all[node_id].last_packet_slot, packet_generation_time = send_packet_generator(Node_all[node_id])
    Node_all[node_id].last_subframe_number = 1
    if calc_cdf
      Node_all[node_id].packet_num_cdf += 1
    end
    #クロックドリフト計算
    CD, time_error = calc_clockdrift(Node_all[node_id].CD_last_calc * 60, packet_generation_time * 60, Node_all[node_id].clockdrift_mean, Node_all[node_id].clockdrift_variance)
    if time_error == 0
      no_error += 1
    else
      error += 1
      # println("クロックドリフト計算エラーnextpacket")
    end
    Node_all[node_id].packet_generation_time = packet_generation_time + CD
    Node_all[node_id].CD_accumulated += CD
    Node_all[node_id].CD_last_calc = Node_all[node_id].packet_generation_time
    # if Node_all[node_id].packet_num % Int(floor(60.0 / packet_period)) == 0
    #   Node_all[node_id].last_backoff = 0
    # end
    push!(event_time, (Node_all[node_id].packet_generation_time - 5 / 60000, 1, node_id))
    push!(event_time, (Node_all[node_id].packet_generation_time, 2, node_id))
    #初期化
    # Node_all[node_id].send_bit = zeros(number_PLIM_bits(length(Node_all[i].prohibited_channel), slot_number))
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

  # CSVファイルに行列を書き出す
  # df = DataFrame(neighborhood_node, :auto)
  # df |> CSV.write(joinpath(output_folder, "neighborhood_node.csv"), delim=',', writeheader=false)

  # df = DataFrame([collect(packet_period:packet_period:SIM_PERIOD), success_packet ./ packet_num], :auto)
  # df |> CSV.write(joinpath(output_folder, "PDR$sim_index.csv"), delim=',', writeheader=false)
  # println("クロックドリフト計算エラー率:", error / (error + no_error))
  return Result([getfield.(Node_all, :x), getfield.(Node_all, :y)], collision_packet, success_packet, rost_packet, packet_num, throughput, cdf, cdf_1, cdf_2, cdf_3, cdf_4, result_node, result_pdf, throughput_cdf, neighborhood_node_sum, neighborhood_node_true_sum, counter_est_interf, counter_only_insys_interf, counter_only_exsys_interf, counter_both_interf, counter_gw_collision_interf, counter_mis_detection, backoff_count, collision_count, double_count, mis_count)
end

runtime = @elapsed begin
  for number_devices = SIM_start_size:SIM_particle_size:number_devices_max
    println("******* start number_devices=$number_devices *******")
    if debug
      result = [main_roop(1, 1000)]
      # result = progress_pmap(x -> main_roop(x, number_devices), 1:SIM_NUM)
      break
    else
      # result = progress_pmap(x -> main_roop(x, number_devices), 1:SIM_NUM, retry_delays=ExponentialBackOff(n=5), on_error=nothing)
      result = progress_pmap(x -> main_roop(x, number_devices), 1:SIM_NUM)
      # nothingを除去
      filtered_results = filter(!isnothing, result)
      length_of_results = length(filtered_results)
      println("シミュレーション成功回数:$length_of_results")
    end
    if !debug
      write_csv(1, filtered_results, number_devices, 1, length_of_results) ## 1:PDR,1:Timeseries
      write_csv(2, filtered_results, number_devices, 1, length_of_results) ## 2:衝突率,1:Timeseries
      write_csv(3, filtered_results, number_devices, 1, length_of_results) ## 3:破棄率,1:Timeseries
      write_csv(4, filtered_results, number_devices, 1, length_of_results) ## 4:スループット,1:Timeseries
      write_csv(5, filtered_results, number_devices, 2, length_of_results) ## 5:端末ごとのPDRのCDF,2:CDF
      write_csv(6, filtered_results, number_devices, 3, length_of_results) ## 6:PDR, 3:端末台数対
      write_csv(7, filtered_results, number_devices, 2, length_of_results) ## 5:端末ごとのPDRのCDF,2:CDF
      write_csv(8, filtered_results, number_devices, 2, length_of_results) ## 5:端末ごとのPDRのCDF,2:CDF
      write_csv(9, filtered_results, number_devices, 2, length_of_results) ## 5:端末ごとのPDRのCDF,2:CDF
      write_csv(10, filtered_results, number_devices, 2, length_of_results) ## 5:端末ごとのPDRのCDF,2:CDF
      write_csv(11, filtered_results, number_devices, 3, length_of_results) ## 4:スループット,3:端末台数対
      write_csv(12, filtered_results, number_devices, 3, length_of_results) ## 4:端末平均スループット,3:端末台数対
      write_csv(13, filtered_results, number_devices, 2, length_of_results) ## 4:端末平均スループット,2:CDF
      write_csv(14, filtered_results, number_devices, 3, length_of_results) ## 干渉推定の結果出力
      write_csv(15, filtered_results, number_devices, 3, length_of_results) ## 干渉推定の結果出力
      write_csv(16, filtered_results, number_devices, 3, length_of_results) ## 干渉推定の結果出力
      write_csv(17, filtered_results, number_devices, 3, length_of_results) ## 干渉推定の結果出力
      write_csv(18, filtered_results, number_devices, 3, length_of_results) ## 干渉推定の結果出力
      write_csv(19, filtered_results, number_devices, 3, length_of_results) ## 無線環境情報推定一致率
      write_csv(20, filtered_results, number_devices, 3, length_of_results) ## 無線環境情報推定精度
      write_csv(21, filtered_results, number_devices, 3, length_of_results) ## 無線環境情報推定精度
      write_csv(22, filtered_results, number_devices, 3, length_of_results) ## 無線環境情報推定精度
      write_csv(23, filtered_results, number_devices, 3, length_of_results) ## 無線環境情報推定精度
      write_csv(24, filtered_results, number_devices, 3, length_of_results) ## 無線環境情報推定精度
      # write_png(result, number_devices)
      # write_png2(result, number_devices)
    end
    GC.gc()
    println("******* end number_devices=$number_devices *******\n")
  end
end
println(runtime)
