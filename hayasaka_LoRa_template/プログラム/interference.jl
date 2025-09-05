#端末のパラメーター
Base.@kwdef mutable struct interference
  node_id::Int64
  status::Int64 = 0 #0:待機，1:CS中，2:Txモード，3:Rxモード
  packet_size::Float64
  sf::Int64
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
  # パケット生起周期
  packet_generation_cycle::Float64
  # 送信決定時刻
  # transmit_time::Float64
  # パケット送信時刻
  packet_generation_time::Float64
  #クロックドリフト平均
  clockdrift_mean::Float64 = rand(Uniform(CD_mean_min, CD_mean_max))
  #クロックドリフト分散
  clockdrift_variance::Float64 = rand(Uniform(CD_variance_min, CD_variance_max))
  #最終クロックドリフト計算時刻
  CD_last_calc::Float64 = 0.0
  #使用可能周波数
  usable_channel::Array{Int64} = [i for i = 1:channel_num]
  #CS中の検知フラグ
  cs_detection::Int64 = 0
  #送信中の衝突フラグ
  tx_collision::Int64 = 0
  #CSにより検知した時間配列
  cs_detection_time::Array{Float64} = []
end


function interference_generator(Interference_all, interference_xy, packet_type)
  for i in eachindex(Interference_all)
    Interference_all[i] = interference(node_id=i + 10000, packet_size=packet_length_type(packet_type, 7), channel=3, sf=7, class="A", packet_generation_cycle=packet_period, packet_generation_time=0.0, x=interference_xy[i][1] / 1000, y=interference_xy[i][2] / 1000)
  end
  return Interference_all
end
