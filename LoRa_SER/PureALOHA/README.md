# LoRa Pure ALOHA シミュレーション

## 概要
LoRa通信におけるPure ALOHAプロトコルのシミュレーションです。CSMA/CAと比較して、より単純なアクセス制御方式の性能を評価します。

## ファイル構成

### メインシミュレーションファイル
- **`PureALOHA_v2.jl`** - 最新版のPure ALOHAシミュレーション
  - キャリアセンス機能付き
  - 高精度な物理層モデル

- **`PureALOHA_CS.jl`** - キャリアセンス付きPure ALOHA
  - CSMA/CAとの性能比較用

- **`PureALOHA_fast.jl`** - 高速実行版
  - 計算最適化済み

- **`PureALOHA.jl`** - 基本版Pure ALOHA
  - シンプルな実装

### 結果ファイル
#### results-CS/
- **`CS_PureALOHA_devices2-8_sf7_iter10000.csv`** - キャリアセンス付き結果
- **`CSMA_SER_sf7_iter10000_dev2-8.csv`** - CSMA/CAとの比較結果

#### results-PureALOHA/
- **`PureALOHA_iter10000_sf7_dev1-2-4-8.csv`** - 基本Pure ALOHA結果
- **`purealoha_iter100_sf7_dev2-8.csv`** - 短時間シミュレーション結果

## 主な機能

### 1. Pure ALOHA
- **ランダムアクセス**: 送信タイミングのランダム化
- **衝突処理**: パケット衝突時の再送制御
- **スループット計算**: 理論値との比較

### 2. キャリアセンス機能
- **送信前監視**: チャネル状態の確認
- **衝突回避**: 干渉の軽減

### 3. 性能評価
- **SER（Symbol Error Rate）**: シンボルエラー率
- **PER（Packet Error Rate）**: パケットエラー率
- **スループット**: データ転送効率

## 実行方法

### 基本実行
```bash
julia PureALOHA.jl
```

### キャリアセンス付き
```bash
julia PureALOHA_CS_v2.jl
```

## 出力結果
- **端末数依存**: 2-8台での性能比較
- **SF依存**: SF7での詳細分析
- **時間依存**: シミュレーション時間による性能変化

## 理論値との比較
- **Pure ALOHA**: 最大スループット 18.4%
- **Slotted ALOHA**: 最大スループット 36.8%
- **CSMA/CA**: より高いスループット（実測値）
