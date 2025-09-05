# using Distributions, StaticArrays, Random, Plots

"""
    poisson_point_process(λ, area)

2次元ポアソン点過程を生成する関数。
λ: 単位面積あたりの平均点の個数
area: 点過程を生成する区画の面積 (例: @SVector [10.0, 10.0])
"""
function poisson_point_process(λ, area)
  # 区画の面積を計算
  area_size = prod(area)

  # ポアソン分布からランダムな点の個数を生成
  num_points = rand(Poisson(λ * area_size))

  # 区画内の一様ランダムな位置に点を配置
  points = [@SVector [area[1] * 2 * rand() - area[1], area[2] * 2 * rand() - area[2]] for _ in 1:num_points]

  return points
end

# # 使用例
# λ = 0.0004   # 単位面積あたりの平均点の個数
# area = @SVector [500.0, 500.0]  # 10m x 10mの区画
# points = poisson_point_process(λ, area)

# println("生成された点の座標:")
# # for p in points
# #   println(p)
# # end
# println(length(points))

# plot_sc_self = scatter([loc[1] for loc in points], [loc[2] for loc in points], color=:blue, legend=false, aspect_ratio=1.0, xlabel="\$x\$", ylabel="\$y\$")

# savefig(plot_sc_self, "node_locations.svg")
