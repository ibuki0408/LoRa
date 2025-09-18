# LoRa PER (Packet Error Rate) Analysis

## 概要
LoRa通信システムのパケットエラー率（PER）を分析するJuliaシミュレーション群です。実際のLoRa復調処理を使用して、様々な通信プロトコルと環境要因を考慮したSNRとPERの関係を詳細に分析します。

## ディレクトリ構成

```
LoRa_PER/
├── LoRa_PL_DC_CS/                    # メインシミュレーションディレクトリ
│   ├── LoRa_Simple_Model.jl         # 基本LoRaシミュレーション
│   ├── LoRa_SlottedALOHA.jl         # SlottedALOHA実装
│   ├── LoRa_SlottedALOHA_DC.jl     # DC制約付きSlottedALOHA
│   ├── LoRa_SlottedALOHA_DC_CS.jl   # DC制約+キャリアセンス付きSlottedALOHA
│   ├── LoRa_SlottedALOHA_DC_CS_ClockDrift.jl  # クロックドリフト考慮版
│   ├── LoRa_SlottedALOHA_Sync.jl    # 時間同期機能付き版
│   ├── Continuous_SNR_Sweep.jl      # 連続SNRスイープ分析
│   ├── test_clock_drift.jl          # クロックドリフトテスト
│   ├── results_LoRa_simple_model/   # 基本シミュレーション結果
│   └── results_SlottedALOHA/        # SlottedALOHA結果
```

## ファイル詳細解説

### 🔧 **基本シミュレーション**

#### **`LoRa_Simple_Model.jl`** - 基本LoRaシミュレーション
- **機能**: LoRa通信の物理層シミュレーション
- **特徴**: 
  - 実際のLoRa変調/復調処理
  - CSMA/DC + Duty Cycle制約
  - シャドウイング効果の考慮
  - SNRスイープ分析
- **用途**: 基本的なLoRa通信性能の評価

#### **`Continuous_SNR_Sweep.jl`** - 連続SNRスイープ分析
- **機能**: 連続的なSNR範囲でのPER分析
- **特徴**:
  - 細かいSNRステップでの分析
  - 統計的評価（複数回実行）
  - CSV形式での結果保存
- **用途**: 詳細なSNR-PER特性の分析

### 📡 **SlottedALOHA実装**

#### **`LoRa_SlottedALOHA.jl`** - 基本SlottedALOHA
- **機能**: 基本的なSlottedALOHAプロトコル
- **特徴**:
  - スロット単位での送信制御
  - 衝突検出とPER計算
  - 復調結果ベースのエラー判定
- **用途**: SlottedALOHAの基本性能評価

#### **`LoRa_SlottedALOHA_DC.jl`** - DC制約付きSlottedALOHA
- **機能**: Duty Cycle制約を考慮したSlottedALOHA
- **特徴**:
  - DC制約管理（1%制約など）
  - DC計算ウィンドウ（1時間）
  - DC制約による送信ブロック
- **用途**: 実用的なDC制約下での性能評価

#### **`LoRa_SlottedALOHA_DC_CS.jl`** - DC制約+キャリアセンス付きSlottedALOHA
- **機能**: DC制約とキャリアセンスを組み合わせたSlottedALOHA
- **特徴**:
  - キャリアセンス機能
  - バックオフ処理
  - 衝突回避メカニズム
- **用途**: 実用的なLoRa通信システムの性能評価

### ⏰ **高度な機能実装**

#### **`LoRa_SlottedALOHA_DC_CS_ClockDrift.jl`** - クロックドリフト考慮版
- **機能**: クロックドリフトを考慮したSlottedALOHA
- **特徴**:
  - 各端末の個別クロックドリフト（±50 ppm）
  - 時間経過に伴う時刻ずれの蓄積
  - スロット開始時刻のずれ考慮
- **用途**: 現実的なクロック精度での性能評価

#### **`LoRa_SlottedALOHA_Sync.jl`** - 時間同期機能付き版
- **機能**: 外部信号による時間同期を考慮したSlottedALOHA
- **特徴**:
  - GPS同期（1μs精度）
  - ネットワーク同期（1ms精度）
  - ビーコン同期（100μs精度）
  - ガードインターバル設定
- **用途**: 実用的な時間同期システムの性能評価

### 🧪 **テスト・分析**

#### **`test_clock_drift.jl`** - クロックドリフトテスト
- **機能**: クロックドリフト実装のテスト
- **特徴**:
  - 短時間シミュレーション
  - クロックドリフト情報の表示
  - 性能影響の確認
- **用途**: クロックドリフト実装の検証

## 主要機能

### 1. **LoRa物理層シミュレーション**
```julia
# LoRa変調/復調
function lora_symbol(sf::Int, bw::Float64, m::Int)
function demod_lora(sf::Int, bw::Float64, x::AbstractVector{ComplexF64})

# AWGN雑音
function add_awgn!(y::AbstractVector{ComplexF64}, snr_dB::Float64)
```

### 2. **伝搬モデル**
```julia
# パスロス計算
function path_loss(distance_km::Float64)
    return 10 * α * log10(distance_km) + β + 10 * γ * log10(f_c)
end

# シャドウイング
function shadowing_value(rng)
    return randn(rng) * shadowing_std
end
```

### 3. **通信プロトコル**

#### **SlottedALOHA**
```julia
# スロット単位での送信制御
for slot in 1:num_slots
    current_time = slot * slot_duration
    # 送信判定とパケット生成
end
```

#### **DC制約管理**
```julia
mutable struct DCConstraint
    dc_limit::Float64          # DC制約（例: 0.01 = 1%）
    dc_remaining::Float64      # 残りDC
    last_reset_time::Float64   # 最後のリセット時刻
    window_duration::Float64   # DC計算ウィンドウ（秒）
end
```

#### **キャリアセンス**
```julia
function carrier_sense_check(device_id::Int, 
                           transmitting_devices::Vector{Int},
                           device_positions::Vector{Tuple{Float64,Float64}},
                           carrier_sense_threshold::Float64)
```

### 4. **クロックドリフト管理**
```julia
mutable struct ClockDrift
    drift_ppm::Float64           # クロックドリフト（ppm）
    accumulated_drift::Float64   # 蓄積された時刻ずれ（秒）
    last_update_time::Float64    # 最後の更新時刻
end
```

### 5. **時間同期管理**
```julia
@enum SyncType GPS_SYNC NETWORK_SYNC BEACON_SYNC

mutable struct TimeSync
    sync_type::SyncType           # 同期方式
    sync_accuracy::Float64        # 同期精度（秒）
    sync_interval::Float64        # 同期間隔（秒）
    guard_interval::Float64       # ガードインターバル（秒）
end
```

## パラメータ設定

### 基本パラメータ
```julia
const area_size = 0.5          # エリアサイズ [km]
const SF = 7                   # スプレディングファクタ
const Tx_dB = 13.0            # 送信電力 [dBm]
const f_c = 923.2             # 搬送波周波数 [MHz]
const band_width = 125e3       # 帯域幅 [Hz]
```

### 伝搬パラメータ
```julia
const α = 4.0                 # パスロス係数
const β = 9.5                 # 伝搬損失オフセット
const γ = 4.5                 # 伝搬周波数係数
const shadowing_std = 3.48    # シャドウイング標準偏差 [dB]
```

### シミュレーション設定
```julia
num_devices = 2               # 端末数
payload_range = (16,32)       # ペイロード長範囲 [シンボル]
sim_time = 10.0              # シミュレーション時間 [秒]
slot_duration = 0.1          # スロット時間 [秒]
dc_limit = 0.01              # DC制約（1%）
```

## 実行方法

### 1. 基本シミュレーション
```bash
julia LoRa_Simple_Model.jl
```

### 2. SlottedALOHAシミュレーション
```bash
julia LoRa_SlottedALOHA.jl
julia LoRa_SlottedALOHA_DC.jl
julia LoRa_SlottedALOHA_DC_CS.jl
```

### 3. クロックドリフト考慮シミュレーション
```bash
julia LoRa_SlottedALOHA_DC_CS_ClockDrift.jl
```

### 4. 時間同期シミュレーション
```bash
julia LoRa_SlottedALOHA_Sync.jl
```

### 5. テスト実行
```bash
julia test_clock_drift.jl
```

## 出力結果

### 1. 可視化ファイル
- **端末配置プロット**: `node_positions_*.png`
- **エラー状態プロット**: `node_positions_errors_*.png`
- **SNR分布プロット**: `snr_distribution_*.png`

### 2. データファイル
- **SNRスイープ結果**: `*_snr_sweep_per_*.csv`
- **統計結果**: コンソール出力

### 3. 結果ディレクトリ
- **`results_LoRa_simple_model/`**: 基本シミュレーション結果
- **`results_SlottedALOHA/`**: SlottedALOHA結果

## 性能比較

### 実装の進化
1. **基本SlottedALOHA** → 基本的なスロット制御
2. **DC制約付き** → 実用的なDC制約考慮
3. **キャリアセンス付き** → 衝突回避機能追加
4. **クロックドリフト考慮** → 現実的な時刻ずれ考慮
5. **時間同期機能** → 外部信号による同期補正

### 主要な改善点
- **PER性能**: クロックドリフト考慮により現実的な性能評価
- **同期精度**: 時間同期により時刻ずれを補正
- **実用性**: 実際のLoRa通信環境に近い条件での評価

## 技術的特徴

### 1. **現実的なモデリング**
- 実際のLoRa変調/復調処理
- 現実的な伝搬モデル
- 実用的な通信プロトコル

### 2. **包括的な分析**
- 物理層からMAC層まで
- 環境要因の考慮
- 統計的評価

### 3. **拡張性**
- モジュール化された設計
- 新機能の追加が容易
- パラメータの柔軟な設定

## 今後の発展

### 1. **高度な同期機能**
- GPS同期の詳細実装
- ネットワーク同期プロトコル
- ビーコン同期システム

### 2. **性能最適化**
- 動的パラメータ調整
- 適応的ガードインターバル
- 最適同期間隔の決定

### 3. **実用化**
- 実際のLoRaWANシステムとの連携
- 実測データとの比較
- 商用システムへの応用

---

このシミュレーション群により、LoRa通信システムの包括的な性能評価が可能です。基本機能から高度な時間同期まで、段階的に機能を拡張しながら、現実的な通信環境での性能を詳細に分析できます。
