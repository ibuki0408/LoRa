#端末のパラメーター
Base.@kwdef mutable struct N_parameters
  node_id::Int64
  status::Int64 = 0 #0:待機，1:CS中，2:Txモード，3:Rxモード
  packet_size::Float64
  sf::Int64
  group_id::Int64
  class::String
  x::Float64
  y::Float64
  P_dB::Int64 = Tx_dB
  channel::Int64
  # shadowing::Float64 = rand(Normal(0, shadowing_standard_deviation), 1)[1]
  # shadowing::Float64 = 0.0
  shadowing::Array{Float64} = []
  #オフセット
  offset::Float64 = rand() * packet_period
  #最終フレーム番号
  last_subframe_number::Int64 = 0
  #最終パケットバックオフ回数
  last_backoff::Int64 = 0
  #送信データ
  send_bit::Array{Int64} = zeros(Int(floor(log2(channel_num * slot_number))))
  # パケット生起周期
  packet_generation_cycle::Float64
  # 送信決定時刻
  # transmit_time::Float64
  #最終パケット送信スロット
  last_packet_slot::Int64 = 0
  # パケット送信時刻
  packet_generation_time::Float64
  #パケット送信成功回数
  send_packet_success::Int64 = 0
  #パケット衝突回数
  packet_collision_count::Int64 = 0
  #SNR規範によりパケット復調失敗
  snr_lost_count::Int64 = 0
  snr_lost_count_inslot::Int64 = 0
  #スロット内パケット衝突回数
  packet_collision_count_inslot::Int64 = 0
  #直前パケット成功フラグ
  last_packet_success::Int64 = 0
  #パケットカウンタ
  packet_num::Int64 = 0
  #パケット破棄回数
  packet_lost_count::Int64 = 0
  #ACKパケットサイズ
  #ACK_packet_size::Float64 = packet_length_type("ackEN")
  #CS中の検知フラグ
  cs_detection::Int64 = 0
  #送信中の衝突フラグ
  tx_collision::Array{Bool} = [false, false]
  #CSにより検知した時間配列
  cs_detection_time::Array{Float64} = []
  #CSにより検知した時間配列に対する，システム内干渉かシステム外干渉かどちらもか(1:システム内，2:システム外)
  cs_detection_type::Array{Array{Bool}} = []
  #CSにより検知した時間配列に対する，真の端末番号
  cs_detection_node::Array{Array{Int64}} = []
  #CSにより検知した時間配列に対する，使用チャネル番号
  cs_detection_channel::Array{Int64} = []
  #クロックドリフト平均
  clockdrift_mean::Float64 = rand(Uniform(CD_mean_min, CD_mean_max))
  #クロックドリフト分散
  clockdrift_variance::Float64 = rand(Uniform(CD_variance_min, CD_variance_max))
  #最終クロックドリフト計算時刻
  CD_last_calc::Float64 = 0.0
  #累積クロックドリフト
  CD_accumulated::Float64 = 0.0
  #使用可能周波数
  usable_channel::Array{Int64} = [i for i = 1:channel_num]

  #CDF計算用
  #パケットカウンタ
  packet_num_cdf::Int64 = 0
  #送信成功回数
  packet_success_cdf::Int64 = 0
end
