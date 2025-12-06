# ===== LoRa Time on Air (ToA) 計算モジュール =====

using Printf

"""
LoRa パラメータ構造体

# フィールド
- `spreading_factor::Int`: 拡散率 (SF7-SF12)
- `bandwidth_hz::Float64`: 帯域幅 (Hz)
- `coding_rate::Int`: 符号化率 (1-4, 実効レートは 4/(4+CR))
- `preamble_length::Int`: プリアンブル長 (シンボル数)
- `payload_size_bytes::Int`: ペイロードサイズ (バイト)
- `crc_enabled::Bool`: CRC 有効化
- `header_enabled::Bool`: 明示的ヘッダー有効化
- `low_data_rate_opt::Bool`: Low Data Rate Optimization 有効化
"""
struct LoRaParameters
    spreading_factor::Int
    bandwidth_hz::Float64
    coding_rate::Int
    preamble_length::Int
    payload_size_bytes::Int
    crc_enabled::Bool
    header_enabled::Bool
    low_data_rate_opt::Bool
end

"""
    create_lora_params(sf::Int, payload_bytes::Int; kwargs...)

LoRa パラメータを生成する便利関数

# 引数
- `sf::Int`: 拡散率 (7-12)
- `payload_bytes::Int`: ペイロードサイズ (バイト)

# キーワード引数
- `bandwidth_hz::Float64=125e3`: 帯域幅 (Hz)
- `coding_rate::Int=1`: 符号化率 (1-4)
- `preamble_length::Int=8`: プリアンブル長
- `crc_enabled::Bool=true`: CRC 有効化
- `header_enabled::Bool=true`: 明示的ヘッダー有効化
- `low_data_rate_opt::Bool=false`: Low Data Rate Optimization

# 例
```julia
# SF10, 20バイトペイロード, デフォルト設定
params = create_lora_params(10, 20)

# SF7, 10バイトペイロード, BW250kHz
params = create_lora_params(7, 10; bandwidth_hz=250e3)
```
"""
function create_lora_params(sf::Int, payload_bytes::Int;
                           bandwidth_hz::Float64=125e3,
                           coding_rate::Int=1,
                           preamble_length::Int=8,
                           crc_enabled::Bool=true,
                           header_enabled::Bool=true,
                           low_data_rate_opt::Bool=false)
    return LoRaParameters(
        sf,
        bandwidth_hz,
        coding_rate,
        preamble_length,
        payload_bytes,
        crc_enabled,
        header_enabled,
        low_data_rate_opt
    )
end

"""
    calculate_lora_airtime(lora_params::LoRaParameters) -> Float64

LoRa パケットの送信時間 (Time on Air) を計算

# 引数
- `lora_params::LoRaParameters`: LoRa パラメータ

# 戻り値
- `Float64`: 送信時間 (ミリ秒)

# 計算式
```
T_symbol = 2^SF / BW
T_preamble = (n_preamble + 4.25) × T_symbol
n_payload = 8 + max(⌈(8PL - 4SF + 28 + 16CRC - 20H) / (4(SF - 2DE))⌉ × (CR + 4), 0)
T_payload = n_payload × T_symbol
T_packet = T_preamble + T_payload
```

# 例
```julia
params = create_lora_params(10, 20)
airtime_ms = calculate_lora_airtime(params)
println("送信時間: \$(airtime_ms) ms")
```
"""
function calculate_lora_airtime(lora_params::LoRaParameters)
    SF = lora_params.spreading_factor
    BW = lora_params.bandwidth_hz
    CR = lora_params.coding_rate
    PL = lora_params.payload_size_bytes
    CRC = lora_params.crc_enabled ? 1 : 0
    H = lora_params.header_enabled ? 0 : 1  # implicit header = 1
    DE = lora_params.low_data_rate_opt ? 1 : 0
    n_preamble = lora_params.preamble_length
    
    # シンボル時間 (秒)
    T_symbol = (2^SF) / BW
    
    # ペイロードシンボル数の計算
    numerator = 8 * PL - 4 * SF + 28 + 16 * CRC - 20 * H
    denominator = 4 * (SF - 2 * DE)
    
    # 分母が0の場合の処理（通常は発生しない）
    if denominator <= 0
        denominator = 1
    end
    
    n_payload_symbols = 8 + max(ceil(numerator / denominator) * (CR + 4), 0)
    
    # 総送信時間 (ミリ秒)
    T_preamble = (n_preamble + 4.25) * T_symbol
    T_payload = n_payload_symbols * T_symbol
    T_packet_ms = (T_preamble + T_payload) * 1000
    
    return T_packet_ms
end

"""
    print_lora_airtime_table(; sf_range=7:12, payload_bytes=20, bandwidth_hz=125e3)

各 SF における送信時間の一覧表を表示

# キーワード引数
- `sf_range`: SF の範囲 (デフォルト: 7:12)
- `payload_bytes::Int`: ペイロードサイズ (デフォルト: 20)
- `bandwidth_hz::Float64`: 帯域幅 (デフォルト: 125e3)

# 例
```julia
print_lora_airtime_table()
print_lora_airtime_table(payload_bytes=10, bandwidth_hz=250e3)
```
"""
function print_lora_airtime_table(; sf_range=7:12, payload_bytes=20, bandwidth_hz=125e3)
    println("="^60)
    println("LoRa Time on Air (ToA) 一覧表")
    println("ペイロード: $(payload_bytes) bytes, 帯域幅: $(bandwidth_hz/1e3) kHz")
    println("="^60)
    println("SF  | 送信時間 (ms) | シンボル時間 (ms)")
    println("-"^60)
    
    for sf in sf_range
        params = create_lora_params(sf, payload_bytes; bandwidth_hz=bandwidth_hz)
        airtime_ms = calculate_lora_airtime(params)
        symbol_time_ms = (2^sf / bandwidth_hz) * 1000
        
        @printf("%-3d | %13.2f | %17.4f\n", sf, airtime_ms, symbol_time_ms)
    end
    
    println("="^60)
end
