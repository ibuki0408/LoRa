# ===== パスロス計算関数 =====

using Statistics

# ===== パスロスパラメータ =====
struct PathLossParameters
    distance_m::Float64       # 送信機と受信機の距離（m）
    frequency_hz::Float64     # 周波数（Hz）
    path_loss_exponent::Float64  # パスロス指数（通常2.0-4.0）
    reference_distance_m::Float64  # 参照距離（m）
    reference_path_loss_db::Float64  # 参照距離でのパスロス（dB）
end

# ===== パスロス計算 =====
function calculate_path_loss(path_loss_params::PathLossParameters)
    distance_ratio = path_loss_params.distance_m / path_loss_params.reference_distance_m
    path_loss_db = path_loss_params.reference_path_loss_db + 
                   10 * path_loss_params.path_loss_exponent * log10(distance_ratio)
    return path_loss_db
end

# ===== パスロス計算（簡易版） =====
function calculate_path_loss_simple(distance_m::Float64, path_loss_exponent::Float64=2.0, 
                                  reference_distance_m::Float64=1.0, reference_path_loss_db::Float64=0.0)
    distance_ratio = distance_m / reference_distance_m
    path_loss_db = reference_path_loss_db + 10 * path_loss_exponent * log10(distance_ratio)
    return path_loss_db
end

# ===== パスロス計算のテスト関数 =====
function test_path_loss()
    println("=== パスロス計算テスト ===")
    
    # パラメータ設定
    path_loss_params = PathLossParameters(
        50.0,               # 距離（m）
        4.7e9,             # 周波数（Hz）
        3.0,               # パスロス指数
        1.0,               # 参照距離（m）
        0.0                # 参照距離でのパスロス（dB）
    )
    
    # パスロス計算
    path_loss_db = calculate_path_loss(path_loss_params)
    
    println("パスロス計算結果:")
    println("• 距離: $(path_loss_params.distance_m) m")
    println("• パスロス指数: $(path_loss_params.path_loss_exponent)")
    println("• パスロス: $(round(path_loss_db, digits=2)) dB")
    
    # 距離によるパスロスの変化
    distances = [10.0, 50.0, 100.0, 200.0, 500.0]
    println("\n距離によるパスロスの変化:")
    for d in distances
        pl = calculate_path_loss_simple(d, 3.0)
        println("• $(d) m → $(round(pl, digits=2)) dB")
    end
    
    return path_loss_db
end

# ===== 実行 =====
if abspath(PROGRAM_FILE) == @__FILE__
    test_path_loss()
end
