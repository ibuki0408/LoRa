# LoRa_Pathloss_CSMA_DC フォルダ

## 概要
このフォルダには、LoRa通信のパスロス、CSMA/CA、デューティサイクル制約を考慮したシミュレーションモデルが含まれています。複雑な要素を排除し、以下の基本機能に焦点を当てたシンプルなモデルです：

- **固定SF値**: 全端末が同じSF（Spreading Factor）を使用
- **CSMA/CA**: キャリアセンス多重アクセス/衝突回避
- **DC制約**: デューティサイクル制約
- **ポアソン点過程**: 端末のランダム配置
- **パスロス**: 距離に基づく伝搬損失
- **SNRスイープ**: 各SNR値でのSER/PER計算

## 主な特徴

### 1. シンプルな設計
- SF値は固定（デフォルト: SF7）
- シャドウイングなし
- 複雑なSNR計算なし

### 2. 端末配置
- ポアソン点過程に基づくランダム配置
- 円形エリア内（デフォルト: 半径0.5km）
- ゲートウェイは原点（0,0）に配置

### 3. 通信プロトコル
- **CSMA/CA**: チャネルがビジーの場合はバックオフ
- **DC制約**: 送信時間に基づく待機時間
- **PER計算**: パケットエラー率の計算

### 4. 伝搬モデル
- パスロスモデル: `PL = 10α*log10(d) + β + 10γ*log10(fc)`
- パラメータ:
  - α = 4.0 (パスロス係数)
  - β = 9.5 (伝搬損失オフセット)
  - γ = 4.5 (伝搬周波数係数)
  - fc = 923.2 MHz (搬送波周波数)

## パラメータ設定

```julia
# 基本パラメータ
const area_size = 0.5        # エリアサイズ(km)
const SF = 7                 # 固定SF値
const Tx_dB = 13.0           # 送信電力(dBm)
const f_c = 923.2            # 搬送波周波数(MHz)

# パスロス係数
const α = 4.0
const β = 9.5
const γ = 4.5

# 通信パラメータ
const band_width = 125e3     # 帯域幅(Hz)
const SNR_threshold = -7.5   # SNR閾値(dB)
```

## 使用方法

### 基本的な実行
```julia
julia LoRa_Simple_Model.jl
```

### カスタムパラメータでの実行
```julia
# シミュレーション実行
per, positions, snrs, errors = simple_csma_dc_per(
    sf=7,                    # SF値
    bw=125e3,               # 帯域幅
    num_devices=50,         # 端末数
    payload_range=(16,32),  # ペイロード長範囲
    sim_time=1.0,           # シミュレーション時間
    backoff_max=0.01,       # 最大バックオフ時間
    dc=0.01                 # デューティサイクル
)
```

## 出力結果

### 1. シミュレーション結果
- 総端末数
- エラー数
- PER（パケットエラー率）
- 平均距離
- 平均SNR
- SNR範囲

### 2. SNRスイープ結果（CSV）
- **ファイル**: `results_simple_model/snr_sweep_ser_sf7_dev50.csv`
- **内容**: 
  - SNR_dB: SNR値（-20.0 ~ 0.0 dB、0.5 dBステップ）
  - SER: シンボルエラー率
  - PER: パケットエラー率
- **用途**: 数値解析、グラフ作成

### 3. 可視化
- 端末配置プロット（エラー状態で色分け）
- ゲートウェイ位置の表示
- PNG形式で保存

### 4. 統計分析
- 複数回実行による統計
- 平均PER、標準偏差
- 最小・最大PER

## フォルダ構成

```
LoRa_Pathloss_CSMA_DC/
├── LoRa_Simple_Model.jl                    # メインシミュレーションファイル（シンプル版）
├── LoRa_CSMA_DC_PathLoss.jl               # 複雑なシミュレーションファイル（元版）
├── README_Simple_Model.md                 # このファイル
└── results_simple_model/                  # 結果保存ディレクトリ
    ├── simple_model_plot.png              # 端末配置プロット
    ├── snr_sweep_ser_sf7_dev50.csv        # SNRスイープ結果（CSV）
    └── ser_vs_snr_sf7_dev50.png           # SER/PER vs SNRプロット（削除済み）
```

## ファイル詳細

### 1. LoRa_Simple_Model.jl
- **目的**: シンプルなLoRaシミュレーションモデル
- **特徴**: 
  - 固定SF値（SF7）
  - ポアソン点過程による端末配置
  - CSMA/CA + DC制約
  - パスロス考慮
  - SNRスイープ機能（SER/PER計算）
- **出力**: CSV形式の数値結果

### 2. LoRa_CSMA_DC_PathLoss.jl
- **目的**: 複雑なLoRaシミュレーションモデル（元版）
- **特徴**:
  - 動的SF選択
  - シャドウイング考慮
  - 詳細な分析機能
  - 可視化機能
- **用途**: より詳細な研究・分析用

### 3. results_simple_model/
- **目的**: シンプルモデルの結果保存
- **内容**:
  - 端末配置プロット
  - SNRスイープ結果（CSV）
  - 各種統計データ

## 元の複雑モデルとの違い

| 要素 | 複雑モデル | シンプルモデル |
|------|------------|----------------|
| SF選択 | 動的（距離・SNRベース） | 固定 |
| シャドウイング | あり | なし |
| 端末配置 | ランダム | ポアソン点過程 |
| SNR計算 | 複雑 | シンプル |
| 分析機能 | 詳細 | 基本 |

## 使用方法の詳細

### シンプルモデルの実行
```bash
cd LoRa_Pathloss_CSMA_DC
julia LoRa_Simple_Model.jl
```

### 複雑モデルの実行
```bash
julia LoRa_CSMA_DC_PathLoss.jl
```

### 結果の確認
```bash
# CSV結果の確認
cat results_simple_model/snr_sweep_ser_sf7_dev50.csv

# プロットの確認
open results_simple_model/simple_model_plot.png
```

## パラメータ調整

### SNRスイープ範囲の変更
```julia
# LoRa_Simple_Model.jl の実行部分で変更
snr_min, snr_max, snr_step = -20.0, 0.0, 0.5  # 現在の設定
# snr_min, snr_max, snr_step = -15.0, 5.0, 1.0  # 例：変更後
```

### 端末数の変更
```julia
num_devices = 50  # 現在の設定
# num_devices = 100  # 例：変更後
```

### 反復回数の変更
```julia
iter_sweep = 50  # 現在の設定
# iter_sweep = 100  # 例：変更後（より正確な結果）
```

## 今後の拡張可能性

1. **SF適応**: 距離に基づくSF選択の追加
2. **シャドウイング**: ログ正規シャドウイングの追加
3. **干渉モデル**: より詳細な干渉計算
4. **移動性**: 端末の移動モデル
5. **エネルギー消費**: バッテリー消費モデル
6. **複数チャネル**: 複数周波数チャネルの考慮

## 注意事項

- このモデルは教育・研究目的で作成されています
- 実際のLoRaデバイスの動作とは異なる場合があります
- パラメータは調整可能ですが、現実的な値の範囲内で設定してください
- 結果の保存先パスは、必要に応じて変更可能です
