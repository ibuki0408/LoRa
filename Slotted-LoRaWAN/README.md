# Power_Peak: Integrated LoRa & 5G Synchronization Simulation

此のプロジェクトは、LoRa通信と5GのOut-of-Band同期（SSB信号を利用した同期）を統合したシミュレーション環境です。

## 主な特徴

- **Out-of-Band同期**: 5G帯（3.7GHz）の同期信号とLoRa帯（920MHz）のデータ通信を分離
- **高精度同期検出**: 線形補間、デバウンス処理、サンプル数フィルタによる正確な同期時刻検出
- **同期精度評価**: 理想ビーコン時刻との誤差計算、統計分析、相関分析
- **スロットベースMAC**: Slotted ALOHAベースのMAC層実装
- **改良型バックオフ**: スロット境界に合わせたバックオフ + チャネル再選択
- **マルチチャネル対応**: AS923準拠の複数チャネル（1-16ch）サポート
- **ポアソン分布起動**: 現実的なIoTデバイス起動シナリオ（指数分布）
- **LoRa ToA自動計算**: SF、ペイロード長から正確なパケット送信時間を計算

## ディレクトリ構成

```
Power_Peak/
├── src/                      # シミュレーションのコアロジック
│   ├── Prop.jl                 # 【提案手法】統合シミュレーションロジック (Sync)
│   ├── Prop_Async.jl           # 【比較手法】非同期シミュレーションロジック (Async Optimized)
│   ├── Pure_ALOHA.jl           # 【比較手法】Pure ALOHAロジック
│   └── modules/                # 機能別モジュール群 (信号生成, パスロス, ToA計算など)
├── scripts/                  # 実行用スクリプト
│   ├── run_prop_per.jl         # 提案手法 (Sync) の実行
│   ├── run_prop_async_per.jl   # 比較手法 (Async) の実行
│   └── run_pure_aloha_per.jl   # 比較手法 (Pure ALOHA) の実行
├── analysis/                 # 分析用スクリプト
│   └── collision_distance_analysis.jl  # 衝突距離・隠れ端末問題の分析
├── results/                  # シミュレーション結果出力先 (CSV, PNG)
│   ├── parallel_evaluation/    # 並列実行結果
│   └── analysis/               # 分析結果
├── figures/                  # グラフ・プロット出力先
└── test/                     # テストコード
```

## 実行方法

### 1. 提案手法（Sync）の実行
統合シミュレーション（同期 + Slotted ALOHA + LBT）を実行します。
```bash
julia scripts/run_prop_per.jl
```

### 2. 比較手法（Async）の実行
非同期シミュレーション（Unslotted ALOHA + LBT）を実行します（旧 `Prop_Async_Optimized.jl`）。
```bash
julia scripts/run_prop_async_per.jl
```

### 3. Pure ALOHAの実行
Pure ALOHA（キャリアセンスなし、非同期）を実行します。
```bash
julia scripts/run_pure_aloha_per.jl
```

## 分析ツールの実行

### 衝突距離分析
衝突した端末ペアの距離分布を分析し、隠れ端末問題の影響を評価します。
```bash
julia analysis/collision_distance_analysis.jl
```

## 主要ファイルの説明

### `src/Prop.jl` (提案手法)
- **Phase 1: 高精度同期**: 5G SSB信号を利用した同期検出
- **Phase 2: スケジューリング**: 同期時刻に基づくスロット割り当て
- **Phase 3: MAC層実行**: スロット同期 + LBT + バックオフ
- **Phase 4: 結果集計**: PER, スループット, 同期精度などの評価

### `src/Prop_Async.jl` (比較手法)
- **非同期動作**: 同期フェーズをスキップ（または強制非同期モード）
- **Unslotted ALOHA**: スロット境界を意識せず送信
- **LBT**: キャリアセンスは有効

### `src/modules/`
- `packet_generation.jl`: ポアソン分布パケット生成
- `collision_detection.jl`: SINR/Capture効果考慮の衝突判定
- `lora_airtime.jl`: LoRa物理層パラメータからToA精算
- `local_clock.jl`: クロックドリフトシミュレーション

## 出力ファイル

### `results/parallel_evaluation/`
各手法ごとにサブディレクトリが作成され、CSV結果が保存されます。
- `mean_per_*`: PER（パケット誤り率）の平均値
- `summary_*`: シミュレーション設定と結果のサマリー

### `results/analysis/`
分析スクリプトの出力結果が保存されます。
