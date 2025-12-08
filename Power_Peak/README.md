# Power_Peak: Integrated LoRa & 5G Synchronization Simulation

このプロジェクトは、LoRa通信と5GのOut-of-Band同期（SSB信号を利用した同期）を統合したシミュレーション環境です。

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
├── integrated_lora_sim.jl    # 【メイン】統合シミュレーション実行スクリプト
├── csma_ca_sim.jl            # CSMA/CA比較用シミュレーション
├── modules/                  # 機能別モジュール群
│   ├── signal_generation.jl    # 信号生成 (5G SSB)
│   ├── path_loss.jl            # パスロス計算 (Log-Distanceモデル)
│   ├── shadowing.jl            # シャドウイング計算 (対数正規分布)
│   ├── noise_generation.jl     # ノイズ生成 (AWGN)
│   ├── terminal_deployment.jl  # 端末配置 (ランダム配置)
│   ├── lora_airtime.jl         # LoRa Time on Air (ToA) 計算
│   ├── collision_detection.jl  # SINR ベース衝突判定
│   ├── packet_generation.jl    # ポアソン分布パケット生成
│   └── local_clock.jl          # クロックドリフト管理
├── result_integrated/        # シミュレーション結果出力先 (CSV, PNG)
└── result_csma_ca/          # CSMA/CA結果出力先
```

## 主要ファイルの説明

### 1. `integrated_lora_sim.jl`
統合シミュレーションのメインスクリプト。以下の処理を実行：

1. **パラメータ設定**
   - 5G同期信号: 3.7GHz, 43dBm
   - LoRaデータ通信: 920MHz, 13dBm
   - マルチチャネル: 8チャネル（AS923準拠）
   - ToA自動計算: SF10, 10バイトペイロード → 約289ms

2. **Phase 1: 高精度同期**
   - ポアソン分布起動（平均15秒、最大30秒）
   - 5G SSB信号受信（35秒観察）
   - 線形補間による正確な同期時刻検出
   - デバウンス処理とサンプル数フィルタ
   - **同期精度分析**: 誤差統計、距離・SNR相関

3. **Phase 2: スケジューリング**
   - 同期時刻に基づくスロット割り当て
   - ポアソン分布によるパケット生成
   - クロックドリフト考慮

4. **Phase 3: MAC層実行**
   - キャリアセンス（LBT）
   - **スロットベースバックオフ**: 1-5スロット後に再試行
   - **チャネル再選択**: バックオフ時に新しいチャネルを選択
   - Duty Cycle制約（1%）

5. **Phase 4: 結果分析**
   - SINR ベース衝突判定
   - PER（Packet Error Rate）計算
   - 統計出力とグラフ生成

### 2. `csma_ca_sim.jl`
CSMA/CA比較用シミュレーション。Integrated LoRa Simとの性能比較に使用。

### 3. `modules/` (モジュール群)
- **`packet_generation.jl`**: ポアソン分布に従ったパケット生成
- **`collision_detection.jl`**: SINR ベースの衝突判定
- **`lora_airtime.jl`**: LoRa物理層パラメータからToA計算
- **`local_clock.jl`**: クロックドリフト（±20ppm）のシミュレーション

## 実行方法

### 統合シミュレーション
```bash
julia integrated_lora_sim.jl
```

### CSMA/CA比較シミュレーション
```bash
julia csma_ca_sim.jl
```

## 出力ファイル

### `result_integrated/`
- `integrated_sync_log_*.csv`: 同期ログ（成功/失敗、検出時刻、誤差、SNR）
- `integrated_sync_accuracy_*.csv`: 同期精度サマリー（統計量、相関係数）
- `integrated_tx_log_*.csv`: 送信ログ（全パケット詳細）
- `integrated_summary_*.csv`: 端末ごとの集計
- `integrated_result.png`: 可視化グラフ
- `integrated_term1_power_*.csv`: 端末1の受信電力データ

### `result_csma_ca/`
- `csma_ca_summary_*.csv`: CSMA/CA結果サマリー

## パフォーマンス

### 最新結果（50端末、8チャネル、10分間）

**Integrated LoRa Sim:**
- 同期成功率: **100%**
- PER: **4.77%**
- 総パケット数: 482
- 成功: 459, 衝突: 23

**CSMA/CA Sim:**
- PER: **4-5%**
- 総パケット数: 450-485

→ **Integrated LoRa SimはCSMA/CAと同等の性能を達成**

## 主要パラメータ

### 物理層
- 5G同期信号: 3.7GHz, 43dBm, 3.6MHz BW
- LoRaデータ: 920MHz, 13dBm, SF10, 125kHz BW

### MAC層
- スロット長: 400ms
- チャネル数: 8（AS923準拠）
- Duty Cycle: 1%
- バックオフ: 1-5スロット（ランダム）
- CS閾値: -80dBm

### シミュレーション
- 端末数: 50
- エリア: 500m × 500m
- 起動時刻: ポアソン分布（平均15秒、最大30秒）
- パケット生成: ポアソン分布（平均60秒間隔）
- シミュレーション時間: 10分

## 研究的特徴

1. **現実的なIoTシナリオ**: ポアソン分布による起動・送信タイミング
2. **高精度同期評価**: 理想値との誤差分析、相関分析
3. **改良型MAC**: スロット境界バックオフ + チャネル再選択
4. **Out-of-Band同期**: 同期とデータを異なる周波数帯で分離
5. **CSMA/CA比較**: 同一条件での性能比較が可能
