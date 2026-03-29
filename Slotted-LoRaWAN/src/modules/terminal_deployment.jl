module TerminalDeployment

using Random, Statistics, Distributions, Printf

export TerminalInfo, TerminalDeploymentParameters, deploy_terminals, analyze_terminal_deployment

# ===== 端末情報 =====
struct TerminalInfo
    terminal_id::Int        # 端末ID
    x_m::Float64            # X座標（m）
    y_m::Float64            # Y座標（m）
    distance_m::Float64     # 基地局からの距離（m）
    path_loss_db::Float64   # パスロス（dB）
    shadowing_db::Float64   # シャドウイング（dB）
    total_loss_db::Float64  # 総損失（パスロス+シャドウイング）（dB）
    rx_power_dbm::Float64   # 受信電力（dBm）
    clock_drift_ppm::Float64  # クロックドリフト（ppm）
    clock_drift_factor::Float64  # ドリフト係数（1.0 + ppm/1e6）
end

# ===== 端末配置パラメータ =====
struct TerminalDeploymentParameters
    deployment_mode::String   # 配置モード: "poisson" または "fixed"
    lambda::Float64           # ポアソン点過程の密度（点/m²）
    num_terminals::Int        # 固定端末数（fixedモード時）
    area_size_m::Float64      # エリアサイズ（m）
    min_distance_m::Float64   # 最小距離（m）
    max_distance_m::Float64  # 最大距離（m）
    frequency_hz::Float64    # 周波数（Hz）
    path_loss_exponent::Float64  # パスロス指数
    reference_distance_m::Float64  # 参照距離（m）
    reference_path_loss_db::Float64  # 参照距離でのパスロス（dB）
    bs_x_m::Float64          # 基地局X座標（m）- 距離計算の基準点
    bs_y_m::Float64          # 基地局Y座標（m）- 距離計算の基準点
end

# ===== パスロス計算（簡易版） =====
function calculate_path_loss_simple(distance_m::Float64, path_loss_exponent::Float64=2.0, 
                                  reference_distance_m::Float64=1.0, reference_path_loss_db::Float64=0.0)
    distance_ratio = distance_m / reference_distance_m
    path_loss_db = reference_path_loss_db + 10 * path_loss_exponent * log10(distance_ratio)
    return path_loss_db
end

# ===== シャドウイング計算（簡易版） =====
function calculate_shadowing_simple(std_db::Float64, enabled::Bool=true)
    if !enabled
        return 0.0
    end
    shadowing_db = randn() * std_db
    return shadowing_db
end

# ===== ポアソン点過程による端末配置 =====
function deploy_terminals_poisson(deployment_params::TerminalDeploymentParameters, 
                                 shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    lambda = deployment_params.lambda
    area_size = deployment_params.area_size_m
    num_points = rand(Poisson(lambda * area_size^2))
    
    for i in 1:num_points
        r_max = deployment_params.max_distance_m
        r_min = deployment_params.min_distance_m
        r = sqrt(r_min^2 + (r_max^2 - r_min^2) * rand())
        theta = 2 * π * rand()
        x = r * cos(theta)
        y = r * sin(theta)
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)
        
        if distance >= deployment_params.min_distance_m && distance <= deployment_params.max_distance_m
            path_loss_db = calculate_path_loss_simple(distance, deployment_params.path_loss_exponent, deployment_params.reference_distance_m, deployment_params.reference_path_loss_db)
            shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)
            total_loss_db = path_loss_db + shadowing_db
            rx_power_dbm = tx_power_dbm - total_loss_db
            drift_ppm = rand(-4.0:0.1:4.0)
            drift_factor = 1.0 + drift_ppm / 1e6
            push!(terminals, TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor))
        end
    end
    return terminals
end

function deploy_terminals_fixed(deployment_params::TerminalDeploymentParameters, 
                                shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    fixed_positions = [(250.0, 250.0), (-50.0, 0.0), (0.0, 50.0), (0.0, -50.0), (35.4, 35.4), (-35.4, 35.4), (35.4, -35.4), (-35.4, -35.4)]
    
    for i in 1:min(deployment_params.num_terminals, length(fixed_positions))
        x, y = fixed_positions[i]
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)
        path_loss_db = calculate_path_loss_simple(distance, deployment_params.path_loss_exponent, deployment_params.reference_distance_m, deployment_params.reference_path_loss_db)
        shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)
        total_loss_db = path_loss_db + shadowing_db
        rx_power_dbm = tx_power_dbm - total_loss_db
        drift_ppm = rand(-4.0:0.1:4.0)
        drift_factor = 1.0 + drift_ppm / 1e6
        push!(terminals, TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor))
    end
    return terminals
end

function deploy_terminals_uniform(deployment_params::TerminalDeploymentParameters, 
                                 shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    num_terminals = deployment_params.num_terminals
    max_distance = min(deployment_params.max_distance_m, deployment_params.area_size_m / 2)
    
    for i in 1:num_terminals
        angle = 2π * (i - 1) / num_terminals
        radius = max_distance * 0.8
        x = radius * cos(angle)
        y = radius * sin(angle)
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)
        path_loss_db = calculate_path_loss_simple(distance, deployment_params.path_loss_exponent, deployment_params.reference_distance_m, deployment_params.reference_path_loss_db)
        shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)
        total_loss_db = path_loss_db + shadowing_db
        rx_power_dbm = tx_power_dbm - total_loss_db
        drift_ppm = rand(-4.0:0.1:4.0)
        drift_factor = 1.0 + drift_ppm / 1e6
        push!(terminals, TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor))
    end
    return terminals
end

function deploy_terminals_random_fixed(deployment_params::TerminalDeploymentParameters,
                                     shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    num_terminals = deployment_params.num_terminals
    r_min = deployment_params.min_distance_m
    r_max = deployment_params.max_distance_m

    for i in 1:num_terminals
        r = sqrt(r_min^2 + (r_max^2 - r_min^2) * rand())
        theta = 2 * π * rand()
        x = r * cos(theta)
        y = r * sin(theta)
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)
        path_loss_db = calculate_path_loss_simple(distance, deployment_params.path_loss_exponent, deployment_params.reference_distance_m, deployment_params.reference_path_loss_db)
        shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)
        total_loss_db = path_loss_db + shadowing_db
        rx_power_dbm = tx_power_dbm - total_loss_db
        drift_ppm = rand(-4.0:0.1:4.0)
        drift_factor = 1.0 + drift_ppm / 1e6
        push!(terminals, TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor))
    end
    return terminals
end

function deploy_terminals(deployment_params::TerminalDeploymentParameters, 
                         shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    if deployment_params.deployment_mode == "poisson"
        return deploy_terminals_poisson(deployment_params, shadowing_std_db, shadowing_enabled, tx_power_dbm)
    elseif deployment_params.deployment_mode == "fixed"
        return deploy_terminals_fixed(deployment_params, shadowing_std_db, shadowing_enabled, tx_power_dbm)
    elseif deployment_params.deployment_mode == "uniform"
        return deploy_terminals_uniform(deployment_params, shadowing_std_db, shadowing_enabled, tx_power_dbm)
    elseif deployment_params.deployment_mode == "random_fixed"
        return deploy_terminals_random_fixed(deployment_params, shadowing_std_db, shadowing_enabled, tx_power_dbm)
    else
        error("Unknown deployment mode: $(deployment_params.deployment_mode)")
    end
end

function analyze_terminal_deployment(terminals::Vector{TerminalInfo})
    if isempty(terminals)
        println("端末が配置されていません。")
        return
    end
    distances = [t.distance_m for t in terminals]
    path_losses = [t.path_loss_db for t in terminals]
    shadowings = [t.shadowing_db for t in terminals]
    rx_powers = [t.rx_power_dbm for t in terminals]
    
    @printf("Terminal Deployment Stats:\n")
    @printf("  Count: %d\n", length(terminals))
    @printf("  Distance: %.1f - %.1f m (Avg: %.1f m)\n", minimum(distances), maximum(distances), mean(distances))
    @printf("  Path Loss: %.1f - %.1f dB (Avg: %.1f dB)\n", minimum(path_losses), maximum(path_losses), mean(path_losses))
    @printf("  Shadowing: %.1f - %.1f dB (Avg: %.1f dB)\n", minimum(shadowings), maximum(shadowings), mean(shadowings))
    @printf("  Rx Power: %.1f - %.1f dBm (Avg: %.1f dBm)\n", minimum(rx_powers), maximum(rx_powers), mean(rx_powers))
end

end # module
