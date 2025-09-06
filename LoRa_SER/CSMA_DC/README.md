# LoRa CSMA/DC シミュレーション

## 概要
LoRa通信におけるCSMA/CA（Carrier Sense Multiple Access with Collision Avoidance）とDC（Duty Cycle）制約を考慮した物理層シミュレーションです。

## ファイル構成

### メインシミュレーションファイル
- **`CSMA_DC_v2.jl`** - 最新版のCSMA/DCシミュレーション
  - LoRa物理層の忠実な再現（チャープ変調、FFT復調）
  - CSMA/CA + DC制約の実装
  - SNRスイープによるPER（Packet Error Rate）評価

- **`CSMA_DC.jl`** - 旧版のCSMA/DCシミュレーション
  - 基本的なCSMA/DC機能

- **`LoRa_CSMA_DC_PER.jl`** - PER評価専用
  - パケットエラー率の詳細分析
  - 物理層変調・復調の実装

- **`LoRa_CSMA_DC_PathLoss.jl`** - パスロス考慮版
  - 端末配置とパスロスモデルの実装
  - 距離依存のSF選択

### 結果ファイル（results/）
- **`CSMA_DC_PER_sf7_iter100_dev1000.csv`** - SF7でのPER結果
- **`csma_dc_SF7_num1000_pkt64-128_results.csv`** - 1000台端末での結果
- **`csma_dc_SF7_num1000_pkt64-128_iter100_summary.csv`** - 統計サマリー

## 主な機能

### 1. LoRa物理層
- **チャープ変調**: LoRaの特徴的なチャープ変調を実装
- **FFT復調**: 実際のLoRa受信機と同様のFFT復調
- **SF対応**: Spreading Factor 7-12に対応

### 2. CSMA/CA
- **キャリアセンス**: 送信前のチャネル監視
- **バックオフ**: 衝突回避のためのランダムバックオフ
- **衝突検知**: パケット衝突の検出と処理

### 3. DC制約
- **Duty Cycle制限**: EU規制に準拠した送信時間制限
- **オフ時間計算**: 送信時間に基づくオフ時間の自動計算

### 4. 伝搬モデル
- **パスロス**: 距離依存の減衰
- **シャドウイング**: 空間相関のあるシャドウイング効果
- **端末配置**: 円形エリア内のランダム配置

## 実行方法

### 基本実行
```bash
julia CSMA_DC_v2.jl
```

### パスロス考慮版
```bash
julia LoRa_CSMA_DC_PathLoss.jl
```

## 出力結果
- **PER vs SNR**: パケットエラー率とSNRの関係
- **端末数依存**: 異なる端末数での性能比較
- **SF依存**: 異なるSpreading Factorでの性能

## パラメータ設定
- **帯域幅**: 125kHz
- **送信電力**: 13dBm
- **搬送波周波数**: 923.2MHz
- **DC制約**: 1%（EU規制準拠）
