# ===== 端末配置関数 =====

using Random, Statistics, Distributions

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
    
    # ポアソン点過程モード（ランダム配置）
    lambda = deployment_params.lambda
    area_size = deployment_params.area_size_m
    num_points = rand(Poisson(lambda * area_size^2))
    println("ポアソン点過程モード: 密度$(lambda)点/m², 生成点数$(num_points)")
    
    for i in 1:num_points
        # 円内のランダムな位置を生成（面積均等分布）
        r_max = deployment_params.max_distance_m
        r_min = deployment_params.min_distance_m
        
        # 面積均等分布のための距離生成
        r = sqrt(r_min^2 + (r_max^2 - r_min^2) * rand())
        theta = 2 * π * rand()  # 0から2πのランダムな角度
        
        # 直交座標に変換
        x = r * cos(theta)
        y = r * sin(theta)
        # BSからの距離を計算
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)
        
        # 距離制限をチェック
        if distance >= deployment_params.min_distance_m && distance <= deployment_params.max_distance_m
            # パスロス計算
            path_loss_db = calculate_path_loss_simple(
                distance, deployment_params.path_loss_exponent,
                deployment_params.reference_distance_m, deployment_params.reference_path_loss_db
            )
            
            # シャドウイング計算
            shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)
            
            # 総損失
            total_loss_db = path_loss_db + shadowing_db
            
            # 受信電力計算
            rx_power_dbm = tx_power_dbm - total_loss_db
            
            # クロックドリフト初期化
            drift_ppm = rand(-4.0:0.1:4.0)
            drift_factor = 1.0 + drift_ppm / 1e6
            
            terminal = TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor)
            push!(terminals, terminal)
            
            println("• 端末$(i): 位置($(round(x, digits=1)), $(round(y, digits=1))) m, 距離$(round(distance, digits=1)) m")
            println("  - パスロス: $(round(path_loss_db, digits=2)) dB")
            println("  - シャドウイング: $(round(shadowing_db, digits=2)) dB")
            println("  - 総損失: $(round(total_loss_db, digits=2)) dB")
            println("  - 受信電力: $(round(rx_power_dbm, digits=1)) dBm")
        end
    end
    
    return terminals
end

# ===== 固定位置による端末配置 =====
function deploy_terminals_fixed(deployment_params::TerminalDeploymentParameters, 
                               shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    
    # 固定端末数モード
    println("固定端末数モード: 端末数$(deployment_params.num_terminals)（固定位置配置）")
    
    # 固定位置を定義
    fixed_positions = [
        (250.0, 250.0),    # 東方向 50m
        (-50.0, 0.0),   # 西方向 50m
        (0.0, 50.0),    # 北方向 50m
        (0.0, -50.0),   # 南方向 50m
        (35.4, 35.4),   # 北東方向 50m
        (-35.4, 35.4),  # 北西方向 50m
        (35.4, -35.4),  # 南東方向 50m
        (-35.4, -35.4)  # 南西方向 50m
    ]
    
    for i in 1:min(deployment_params.num_terminals, length(fixed_positions))
        x, y = fixed_positions[i]
        # BSからの距離を計算
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)
        
        # パスロス計算
        path_loss_db = calculate_path_loss_simple(
            distance, deployment_params.path_loss_exponent,
            deployment_params.reference_distance_m, deployment_params.reference_path_loss_db
        )
        
        # シャドウイング計算
        shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)
        
        # 総損失
        total_loss_db = path_loss_db + shadowing_db
        
        # 受信電力計算
        rx_power_dbm = tx_power_dbm - total_loss_db
        
        # クロックドリフト初期化
        drift_ppm = rand(-4.0:0.1:4.0)
        drift_factor = 1.0 + drift_ppm / 1e6
        
        terminal = TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor)
        push!(terminals, terminal)
        
        println("• 端末$(i): 位置($(x), $(y)) m, 距離$(round(distance, digits=1)) m")
        println("  - パスロス: $(round(path_loss_db, digits=2)) dB")
        println("  - シャドウイング: $(round(shadowing_db, digits=2)) dB")
        println("  - 総損失: $(round(total_loss_db, digits=2)) dB")
        println("  - 受信電力: $(round(rx_power_dbm, digits=1)) dBm")
    end
    
    return terminals
end

# ===== 固定数均等配置 =====
function deploy_terminals_uniform(deployment_params::TerminalDeploymentParameters, 
                                 shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    terminals = TerminalInfo[]
    
    # 固定数均等配置モード
    num_terminals = deployment_params.num_terminals
    max_distance = min(deployment_params.max_distance_m, deployment_params.area_size_m / 2)
    
    println("固定数均等配置モード: 端末数$(num_terminals)（面積均等配置）")
    
    # 円周上に均等配置
    for i in 1:num_terminals
        angle = 2π * (i - 1) / num_terminals  # 均等な角度
        radius = max_distance * 0.8  # 最大距離の80%の位置
        x = radius * cos(angle)
        y = radius * sin(angle)
        # BSからの距離を計算
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)
        
        # パスロス計算
        path_loss_db = calculate_path_loss_simple(
            distance, deployment_params.path_loss_exponent,
            deployment_params.reference_distance_m, deployment_params.reference_path_loss_db
        )
        
        # シャドウイング計算
        shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)
        
        # 総損失
        total_loss_db = path_loss_db + shadowing_db
        
        # 受信電力計算
        rx_power_dbm = tx_power_dbm - total_loss_db
        
        # クロックドリフト初期化
        drift_ppm = rand(-4.0:0.1:4.0)
        drift_factor = 1.0 + drift_ppm / 1e6
        
        terminal = TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor)
        push!(terminals, terminal)
        
        println("• 端末$(i): 位置($(round(x, digits=1)), $(round(y, digits=1))) m, 距離$(round(distance, digits=1)) m")
        println("  - パスロス: $(round(path_loss_db, digits=2)) dB")
        println("  - シャドウイング: $(round(shadowing_db, digits=2)) dB")
        println("  - 総損失: $(round(total_loss_db, digits=2)) dB")
        println("  - 受信電力: $(round(rx_power_dbm, digits=1)) dBm")
    end
    
    return terminals
end

# ===== 固定数ランダム配置（面積均等・ポアソン風、個数固定） =====
function deploy_terminals_random_fixed(deployment_params::TerminalDeploymentParameters,
                                     shadowing_std_db::Float64, shadowing_enabled::Bool, tx_power_dbm::Float64)
    terminals = TerminalInfo[]

    num_terminals = deployment_params.num_terminals
    r_min = deployment_params.min_distance_m
    r_max = deployment_params.max_distance_m

    println("固定数ランダム配置モード: 端末数$(num_terminals)（面積均等ランダム）")

    for i in 1:num_terminals
        # 面積均等（円環一様）に半径をサンプルし、角度は一様
        r = sqrt(r_min^2 + (r_max^2 - r_min^2) * rand())
        theta = 2 * π * rand()

        x = r * cos(theta)
        y = r * sin(theta)
        # BSからの距離を計算
        distance = sqrt((x - deployment_params.bs_x_m)^2 + (y - deployment_params.bs_y_m)^2)

        # パスロス計算
        path_loss_db = calculate_path_loss_simple(
            distance, deployment_params.path_loss_exponent,
            deployment_params.reference_distance_m, deployment_params.reference_path_loss_db
        )

        # シャドウイング計算
        shadowing_db = calculate_shadowing_simple(shadowing_std_db, shadowing_enabled)

        # 総損失
        total_loss_db = path_loss_db + shadowing_db

        # 受信電力計算
        rx_power_dbm = tx_power_dbm - total_loss_db

        # クロックドリフト初期化
        drift_ppm = rand(-4.0:0.1:4.0)
        drift_factor = 1.0 + drift_ppm / 1e6
        
        terminal = TerminalInfo(i, x, y, distance, path_loss_db, shadowing_db, total_loss_db, rx_power_dbm, drift_ppm, drift_factor)
        push!(terminals, terminal)

            # println("• 端末$(i): 位置($(round(x, digits=1)), $(round(y, digits=1))) m, 距離$(round(distance, digits=1)) m")
            # println("  - パスロス: $(round(path_loss_db, digits=2)) dB")
            # println("  - シャドウイング: $(round(shadowing_db, digits=2)) dB")
            # println("  - 総損失: $(round(total_loss_db, digits=2)) dB")
            # println("  - 受信電力: $(round(rx_power_dbm, digits=1)) dBm")
    end

    return terminals
end

# ===== 統合端末配置関数 =====
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

# ===== 端末配置統計分析 =====
function analyze_terminal_deployment(terminals::Vector{TerminalInfo})
    if isempty(terminals)
        println("端末が配置されていません。")
        return
    end
    
    distances = [t.distance_m for t in terminals]
    path_losses = [t.path_loss_db for t in terminals]
    shadowings = [t.shadowing_db for t in terminals]
    rx_powers = [t.rx_power_dbm for t in terminals]
    
    println("端末配置統計分析:")
    println("• 端末数: $(length(terminals))")
    println("• 距離範囲: $(round(minimum(distances), digits=1)) - $(round(maximum(distances), digits=1)) m")
    println("• 平均距離: $(round(mean(distances), digits=1)) m")
    println("• パスロス範囲: $(round(minimum(path_losses), digits=1)) - $(round(maximum(path_losses), digits=1)) dB")
    println("• 平均パスロス: $(round(mean(path_losses), digits=1)) dB")
    println("• シャドウイング範囲: $(round(minimum(shadowings), digits=1)) - $(round(maximum(shadowings), digits=1)) dB")
    println("• 平均シャドウイング: $(round(mean(shadowings), digits=1)) dB")
    println("• 受信電力範囲: $(round(minimum(rx_powers), digits=1)) - $(round(maximum(rx_powers), digits=1)) dBm")
    println("• 平均受信電力: $(round(mean(rx_powers), digits=1)) dBm")
end

# ===== 端末配置のテスト関数 =====
function test_terminal_deployment()
    Random.seed!(1234)  # ランダムシードを固定
    println("=== 端末配置テスト ===")
    
    # パラメータ設定
    deployment_params = TerminalDeploymentParameters(
        "poisson",          # ポアソン点過程モード
        0.001,              # 密度（点/m²）
        0,                  # 固定端末数（使用しない）
        100.0,              # エリアサイズ（m）
        10.0,               # 最小距離（m）
        500.0,              # 最大距離（m）
        4.7e9,              # 周波数（Hz）
        3.0,                # パスロス指数
        1.0,                # 参照距離（m）
        0.0,                # 参照距離でのパスロス（dB）
        0.0,                # bs_x_m（原点）
        0.0                 # bs_y_m（原点）
    )
    
    # 端末配置
    terminals = deploy_terminals(deployment_params, 8.0, true, 43.0)
    
    # 統計分析
    analyze_terminal_deployment(terminals)
    
    return terminals
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    test_terminal_deployment()
end
