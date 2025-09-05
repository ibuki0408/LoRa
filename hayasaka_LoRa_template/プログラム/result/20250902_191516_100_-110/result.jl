Base.@kwdef mutable struct Result
  #端末位置
  node_xy::Array{Array{Float64}}

  #総パケット数
  # total_packets::Int64

  #全パケットデータ
  # packet_all

  #衝突パケット数
  collision::Array{Int64}
  # collision

  #正常送信パケット数
  success::Array{Int64}
  # success

  #破棄パケット数
  rost::Array{Int64}
  # rost

  #各範囲におけるパケット数
  packet_num::Array{Int64}
  # packet_num

  #スループット
  throughput::Array{Int64}

  #PDR_CDF
  pdr_cdf

  #PDR_CDF1
  pdr_cdf1
  #PDR_CDF2
  pdr_cdf2
  #PDR_CDF3
  pdr_cdf3
  #PDR_CDF4
  pdr_cdf4

  #推定ヒット数
  # est_hit

  #推定全数
  # est_all

  #各クラスタ端末台数
  node_num::Array{Int64}

  #各クラスタPDF
  node_pdf

  #スループットCDF
  cdf_throu

  #推定一致率計測
  est_sum
  est_true_sum

  #干渉
  est_interf
  only_insys_interf
  only_exsys_interf
  both_interf
  gw_collision_interf
  mis_detection

  #推定精度
  backoff_n
  collision_n
  double_n
  mis_n
end
