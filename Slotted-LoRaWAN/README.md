# SlottedLoRaWAN: LoRaWAN Slotted ALOHA Simulation

このプロジェクトは、LoRaWANにおけるSlotted ALOHAの性能を評価するためのシミュレーション環境です。5GのSSB信号を利用したOut-of-Band同期メカニズムを統合し、高精度なスロット同期を実現する機能を備えています。

## 主な特徴

- **スロット同期 (Slotted ALOHA)**: 5G帯（3.7GHz）の同期信号（SSB）を利用した高精度なLoRa通信（920MHz）同期
- **高精度同期検出**: 線形補間、デバウンス処理、サンプル数フィルタによる正確な同期時刻検出
- **同期精度評価**: 理想ビーコン時刻との誤差計算、統計分析、相関分析
- **LBT (Listen Before Talk)**: キャリアセンスに基づく送信制御
- **改良型バックオフ**: スロット境界に合わせたバックオフ + チャネル再選択
- **マルチチャネル対応**: AS923準拠の複数チャネル（1-16ch）サポート
- **ポアソン分布起動**: 現実的なIoTデバイス起動シナリオ（指数分布）
- **LoRa ToA自動計算**: SF、ペイロード長から正確なパケット送信時間を計算

## ディレクトリ構成

```
Slotted-LoRaWAN/
├── src/                      # シミュレーションのコアロジック
│   ├── SlottedLoRaWAN.jl     # メインモジュール・エントリポイント
│   ├── Prop.jl               # 【提案手法】統合シミュレーション（同期あり）
│   ├── Prop_Async.jl         # 【比較手法】非同期シミュレーション（LBTあり）
│   ├── Pure_ALOHA.jl         # 【比較手法】Pure ALOHAロジック
│   └── modules/              # 機能別モジュール群 (信号生成, パスロス, ToA計算など)
├── scripts/                  # 実行用スクリプト
│   ├── run_prop_per.jl       # 提案手法 (Sync) の実行
│   ├── run_prop_async_per.jl # 比較手法 (Async) の実行
│   └── run_pure_aloha_per.jl # 比較手法 (Pure ALOHA) の実行
├── theory/                   # 理論分析スクリプト
│   ├── theory.jl             # スロット同期方式の理論解析
│   ├── theory_ppp.jl         # PPP（Poisson Point Process）に基づく解析
│   └── pure_aloha.jl         # Pure ALOHAの理論解析
├── results/                  # シミュレーション結果出力先 (CSV, PNG)
│   ├── parallel_evaluation/  # 並列実行結果
│   └── analysis/             # 分析結果
├── figure/                   # グラフ・プロット出力先
├── compare.jl                # 理論値とシミュレーション値の比較スクリプト
└── Project.toml              # Juliaプロジェクト設定
```

## 実行方法

### 1. 提案手法（Sync）の実行
統合シミュレーション（同期 + Slotted ALOHA + LBT）を実行します。
```bash
julia scripts/run_prop_per.jl
```

### 2. 比較手法（Async）の実行
非同期シミュレーション（Unslotted ALOHA + LBT）を実行します。
```bash
julia scripts/run_prop_async_per.jl
```

### 3. Pure ALOHAの実行
Pure ALOHA（キャリアセンスなし、非同期）を実行します。
```bash
julia scripts/run_pure_aloha_per.jl
```

### 4. 理論値との比較
シミュレーション結果と理論解析結果を比較します。
```bash
julia compare.jl
```

## 主要ファイルの説明

### `src/SlottedLoRaWAN.jl`
プロジェクトのメインモジュールです。必要なサブモジュールを統合し、公開関数や型を定義しています。

### `src/Prop.jl` (提案手法)
- **Phase 1: 高精度同期**: 5G SSB信号を利用した同期検出
- **Phase 2: スケジューリング**: 同期時刻に基づくスロット割り当て
- **Phase 3: MAC層実行**: スロット同期 + LBT + バックオフ
- **Phase 4: 結果集計**: PER, スループット, 同期精度などの評価

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

### `figure/`
シミュレーションや理論解析によって生成されたグラフが保存されます。
