# LoRa Simple Simulation

## 概要
RCS_test.jlを参考にして、クラスタリングと干渉推定機能を除いたシンプルなLoRa通信シミュレーションです。基本的なLoRa通信機能に焦点を当てています。

## ファイル構成

### メインシミュレーションファイル
- **`LoRa_Simple_Simulation.jl`** - メインシミュレーションファイル
  - 基本的なLoRa通信シミュレーション
  - キャリアセンス機能
  - パスロスとシャドウイングモデル
  - クラスタリング機能なし（簡略化版）

### 結果ファイル（result/）
- **`PDR-50.csv`** - 50台端末でのパケット配信率
- **`PDR-100.csv`** - 100台端末でのパケット配信率
- **`Collision-50.csv`** - 50台端末での衝突率
- **`Collision-100.csv`** - 100台端末での衝突率
- **`Rost-50.csv`** - 50台端末での破棄率
- **`Rost-100.csv`** - 100台端末での破棄率
- **`Throughput-50.csv`** - 50台端末でのスループット
- **`Throughput-100.csv`** - 100台端末でのスループット

## 主な機能

### 1. LoRa通信シミュレーション
- **キャリアセンス**: CSMA/CA機能
- **バックオフ**: 衝突回避制御
- **マルチチャネル**: 4チャネル対応
- **SNR/SIR判定**: 信号品質判定

### 2. 伝搬モデル
- **パスロス**: 距離依存減衰
- **シャドウイング**: 空間相関シャドウイング
- **クロックドリフト**: 端末間の時刻ずれ

### 3. 性能評価
- **PDR**: パケット配信率
- **衝突率**: パケット衝突率
- **破棄率**: パケット破棄率
- **スループット**: データ転送効率

## 実行方法

### 基本実行
```bash
julia LoRa_Simple_Simulation.jl
```

### 並列実行
```bash
julia -p auto LoRa_Simple_Simulation.jl
```

## 出力結果

### データ形式
```csv
時間(分), 指標値
2.0, 0.9975664165483674
4.0, 0.992866696388765
6.0, 0.9945454545454545
...
4320.0, 0.9931672025723474
```

### 指標の説明
- **PDR**: 成功パケット数 ÷ 総パケット数
- **衝突率**: 衝突パケット数 ÷ 総パケット数
- **破棄率**: 破棄パケット数 ÷ 総パケット数
- **スループット**: 総データ量 ÷ (パケット周期 × 60秒)

## パラメータ設定

### シミュレーション設定
- **シミュレーション時間**: 3日間（4320分）
- **端末台数**: 50台、100台
- **シミュレーション回数**: 1000回
- **パケット送信周期**: 2分間隔

### LoRa設定
- **帯域幅**: 125kHz
- **送信電力**: 13dBm
- **搬送波周波数**: 923.2MHz
- **SF範囲**: 7-10
- **チャネル数**: 4

### 伝搬モデル
- **パスロス係数**: α=4.0, β=9.5, γ=4.5
- **シャドウイング**: 標準偏差3.48dB
- **相関距離**: 50m

## 特徴
- **シンプル**: クラスタリング機能を除去
- **基本機能**: LoRa通信の基本機能に集中
- **理解しやすい**: コード構造が簡潔
- **学習用途**: LoRa通信の理解に適している

## RCS_test.jlとの違い
- **クラスタリング機能**: 削除
- **干渉推定機能**: 削除
- **隠れ端末推定**: 削除
- **基本通信機能**: 保持
- **伝搬モデル**: 保持

## **Gitプッシュの手順**

### **1. 現在の状態確認**
```bash
git status
```

### **2. 変更されたファイルをステージング**
```bash
# すべての変更を追加
git add .

# または特定のファイルのみ
git add README.md
git add LoRa_SER/CSMA_DC/README.md
git add hayasaka_LoRa_template/README.md
git add LoRa_Simulation/README.md
git add LoRa_SER/PureALOHA/README.md
git add LoRa_SER/LoRa_uncompleted/README.md
```

### **3. コミット**
```bash
git commit -m "Add README files for each LoRa simulation folder

- Added comprehensive documentation for LoRa_SER/CSMA_DC
- Added documentation for hayasaka_LoRa_template
- Added documentation for LoRa_Simulation
- Added documentation for PureALOHA simulations
- Added documentation for LoRa_uncompleted simulations
- Added documentation for Distribution folder"
```

### **4. プッシュ**
```bash
git push origin main
```

### **5. エラーが発生した場合**

#### **リモートブランチが進んでいる場合**
```bash
git pull origin main
git push origin main
```

#### **初回プッシュの場合**
```bash
git push -u origin main
```

### **6. 確認**
```bash
git log --oneline -5
git status
```

## **推奨手順**

1. **まず現在の状態を確認**：
   ```bash
   git status
   ```

2. **変更をステージング**：
   ```bash
   git add .
   ```

3. **コミット**：
   ```bash
   git commit -m "Add comprehensive README documentation for all LoRa simulation folders"
   ```

4. **プッシュ**：
   ```bash
   git push origin main
   ```

これで、作成したREADMEファイルとその他の変更がGitHubにプッシュされます。
