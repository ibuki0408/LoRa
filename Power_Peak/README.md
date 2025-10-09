# Power_Peak - 同期信号受信シミュレーション

このディレクトリには、無線通信における同期信号の受信とピーク検出をシミュレーションするJuliaプログラムが含まれています。

## 📁 ファイル構成

### メインシミュレーションファイル

#### 1. `sync_signal_simulation.jl`
- **目的**: 基本的な同期信号受信シミュレーション
- **特徴**: 
  - 単一サンプリングレート
  - パスロス、シャドウイング、ノイズを考慮
  - 動的ピーク検出
- **出力**: `results_sync_simulation/` ディレクトリ

#### 2. `sync_fixed_terminal.jl`
- **目的**: 固定端末位置での同期信号受信シミュレーション
- **特徴**:
  - 端末位置を固定座標で定義
  - 同期信号と端末の帯域幅を分離設定
  - 再現性の高いシミュレーション
- **出力**: `results_sync_simulation_fixed/` ディレクトリ

#### 3. `sync_resampling.jl`
- **目的**: 送信側と受信側で異なるサンプリングレートを使用するシミュレーション
- **特徴**:
  - 送信側: 高サンプリングレート（理想信号生成）
  - 受信側: 低サンプリングレート（効率的受信）
  - `DSP.jl`の`resample`関数を使用
- **出力**: `results_sync_simulation_resampling/` ディレクトリ

### 参考ファイル

#### 4. `continuous_QPSK_simulation.jl`
- **目的**: 連続QPSK信号生成の参考実装
- **特徴**: 基本的なOFDM信号生成手法

#### 5. `LoRa_Simple_Simulation_hayasaka.jl`
- **目的**: LoRa通信の参考実装
- **特徴**: LoRa特有のパスロスモデル

## 🔧 主要機能

### 信号生成
- **OFDM信号**: QPSK変調による同期信号
- **周期的送信**: 設定可能な間隔での信号送信
- **帯域幅制御**: 送信側と受信側で異なる帯域幅設定

### 伝搬モデル
- **パスロス**: 自由空間パスロスモデル
- **シャドウイング**: 対数正規分布によるランダム減衰
- **ノイズ**: AWGN（加法性白色ガウス雑音）

### 端末配置
- **固定配置**: 事前定義された座標での端末配置
- **ランダム配置**: ポアソン点過程によるランダム配置
- **円形配置**: 基地局中心の円形エリア内配置

### ピーク検出
- **動的閾値**: 平均電力に基づく適応的閾値
- **時間窓**: 期待信号タイミング周辺での検索
- **検出率**: 統計的な検出性能評価

## 📊 シミュレーションパラメータ

### 信号パラメータ
```julia
signal_duration_us::Float64      # 信号持続時間（μs）
center_frequency_ghz::Float64    # 中心周波数（GHz）
signal_bandwidth_mhz::Float64    # 同期信号帯域幅（MHz）
terminal_bandwidth_mhz::Float64  # 端末受信帯域幅（MHz）
tx_sampling_rate_mhz::Float64    # 送信側サンプリングレート（MHz）
rx_sampling_rate_mhz::Float64    # 受信側サンプリングレート（MHz）
tx_power_dbm::Float64            # 送信電力（dBm）
```

### 受信環境パラメータ
```julia
snr_db::Float64                  # 信号対雑音比（dB）
shadowing_enabled::Bool          # シャドウイング有効/無効
shadowing_std_db::Float64        # シャドウイング標準偏差（dB）
```

### 端末配置パラメータ
```julia
deployment_mode::String          # 配置モード: "fixed" または "poisson"
num_terminals::Int               # 端末数
area_size_m::Float64            # エリアサイズ（m）
path_loss_exponent::Float64      # パスロス指数
```

## 🚀 実行方法

### 基本的な実行
```bash
# 基本シミュレーション
julia sync_signal_simulation.jl

# 固定端末シミュレーション
julia sync_fixed_terminal.jl

# リサンプリング対応シミュレーション
julia sync_resampling.jl
```

### パラメータ変更
各ファイル内の`create_simulation_parameters()`関数でパラメータを変更できます：

```julia
function create_simulation_parameters()
    return SimulationParameters(
        # 信号パラメータ
        142.8,              # 信号持続時間（μs）
        4.7,                # 中心周波数（GHz）
        1.0,                # 同期信号帯域幅（MHz）
        0.01,               # 端末受信帯域幅（MHz）
        2.0,                # 送信側サンプリングレート（MHz）
        0.02,               # 受信側サンプリングレート（MHz）
        20.0,               # 送信電力（dBm）
        # ... その他のパラメータ
    )
end
```

## 📈 出力結果

### CSVファイル
- **`received_power_data_*.csv`**: 受信電力データ（時間、電力）
- **`sync_simulation_parameters_*.csv`**: シミュレーションパラメータ

### コンソール出力
- シミュレーションパラメータ
- 端末配置情報
- 受信電力ピーク検出結果
- 検出率とピーク電力範囲

## 🔬 技術的特徴

### リサンプリング技術
`sync_resampling.jl`では、`DSP.jl`の`resample`関数を使用して：
- 送信側: 高サンプリングレート（16.0 MHz）で理想信号生成
- 受信側: 低サンプリングレート（0.02 MHz）で効率的受信
- エイリアシング防止: 適切なフィルタリングによる品質保持

### 帯域幅分離
- **同期信号帯域幅**: 広帯域信号（1.0 MHz）
- **端末受信帯域幅**: 狭帯域受信（0.01 MHz）
- **ノイズ低減**: 狭帯域受信による低ノイズ電力

### 動的ピーク検出
- **適応的閾値**: 平均電力の3倍を基準
- **時間窓検索**: 期待信号タイミング±10%の範囲
- **統計評価**: 検出率とピーク電力分布

## 📚 依存パッケージ

```julia
using Random, Statistics, Printf, FFTW, LinearAlgebra, DSP, 
      Distributions, CSV, DataFrames, Plots, Dates
```

## 🎯 使用例

### 1. 基本的な同期信号検出
```bash
julia sync_signal_simulation.jl
```

### 2. 固定端末での検出性能評価
```bash
julia sync_fixed_terminal.jl
```

### 3. 異なるサンプリングレートでの現実的シミュレーション
```bash
julia sync_resampling.jl
```

## 📝 注意事項

1. **再現性**: `Random.seed!(1234)`で固定シードを使用
2. **メモリ使用量**: 高サンプリングレートでは大量のメモリが必要
3. **実行時間**: リサンプリング処理により実行時間が増加
4. **出力ディレクトリ**: 各シミュレーションで異なるディレクトリに保存

## 🔧 カスタマイズ

### パラメータ調整
- **SNR**: `snr_db`パラメータで調整
- **端末数**: `num_terminals`で変更
- **検出閾値**: `power_threshold_dbm`で調整
- **帯域幅**: `signal_bandwidth_mhz`と`terminal_bandwidth_mhz`で分離設定

### 新しい機能の追加
- マルチパス伝搬の実装
- 異なる変調方式の追加
- より複雑な端末配置パターン
- リアルタイム処理の最適化

## 📊 結果の解釈

### 検出率
- **100%**: 全ての同期信号を検出
- **0%**: 信号が検出されない（閾値が高すぎる可能性）
- **部分検出**: 信号品質とノイズレベルのバランス

### 受信電力
- **高い値**: 良好な受信環境
- **低い値**: パスロスやシャドウイングの影響
- **変動**: ノイズとシャドウイングによる統計的変動

このシミュレーションシステムは、無線通信における同期信号の検出性能を評価し、様々な通信環境での動作を理解するための包括的なツールです。
