### STRUCT ###

@kwdef mutable struct s_wave
  f_x::Float64 = 0.0
  f_y::Float64 = 0.0
end

@kwdef mutable struct F_T
  f_T::s_wave = s_wave()
  f_R::s_wave = s_wave()
  theta::Float64 = 0.0
end

### FUNCTION ###
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
