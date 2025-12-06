# Power_Peak: Integrated LoRa & 5G Synchronization Simulation

このプロジェクトは、LoRa通信と5GのOut-of-Band同期（SSB信号を利用した同期）を統合したシミュレーション環境です。

## ディレクトリ構成

```
Power_Peak/
├── integrated_lora_sim.jl    # 【メイン】統合シミュレーション実行スクリプト
├── main_simulation.jl        # (旧) 開発用メインスクリプト
├── modules/                  # 機能別モジュール群
│   ├── signal_generation.jl    # 信号生成 (5G SSB, LoRaチャープ)
│   ├── path_loss.jl            # パスロス計算 (Log-Distanceモデル)
│   ├── shadowing.jl            # シャドウイング計算 (対数正規分布)
│   ├── noise_generation.jl     # ノイズ生成 (AWGN)
│   ├── terminal_deployment.jl  # 端末配置 (ランダム, 固定, 円形など)
│   ├── lora_airtime.jl         # [NEW] LoRa Time on Air (ToA) 計算
│   └── local_clock.jl          # (未使用) ローカルクロック管理
├── result_integrated/        # シミュレーション結果出力先 (CSV, PNG)
└── archive/                  # 過去のスクリプトやバックアップ
```

## 主要ファイルの説明

### 1. `integrated_lora_sim.jl`
シミュレーションのメインエントリポイントです。以下の処理を一括で行います：
1.  **パラメータ設定**: 5G同期信号(3.7GHz)とLoRaデータ通信(920MHz)のパラメータ定義
    *   **ToA自動計算**: SF（拡散率）とペイロード長からパケット送信時間（Time on Air）を自動計算
2.  **端末配置**: 指定されたエリア内に端末を配置
3.  **Phase 1 (同期)**: 5G基地局からのSSB信号を受信し、各端末の同期タイミング（スロット開始時刻）を決定
4.  **Phase 2 (スケジューリング)**: 同期タイミングに基づき、各端末の送信スロットを割り当て
5.  **Phase 3 (MAC層)**: CSMA/CA（キャリアセンス）を用いたデータ送信シミュレーション
6.  **結果出力**: ログ(CSV)と可視化グラフ(PNG)を `result_integrated/` に保存

### 2. `modules/` (モジュール群)
機能ごとにファイルを分割し、再利用性を高めています。
*   **`lora_airtime.jl`**: LoRaの物理層パラメータ（SF, BW, CR等）に基づいて、正確なパケット送信時間（ToA）を計算します。
*   **`terminal_deployment.jl`**: 端末の配置ロジック。「円形ランダム配置(random_fixed)」などが定義されています。
*   **`path_loss.jl`**: 距離と周波数に応じた電波減衰を計算します（Log-Distance Path Loss Model）。

## 実行方法

```bash
julia integrated_lora_sim.jl
```

実行すると、コンソールにシミュレーション経過と**計算されたToA**が表示され、`result_integrated/` フォルダに結果ファイルが生成されます。

## シミュレーションの流れ

1.  **初期設定 & ToA計算**
    *   SF (Spreading Factor) とペイロードサイズから、パケット送信時間を自動計算
2.  **5G同期信号の受信** (3.7GHz, 43dBm)
    *   端末はランダムな時刻に起動 (`Startup Time`)
    *   起動後、最初に来たビーコン信号を検出して同期
3.  **LoRaデータ通信** (920MHz, 13dBm)
    *   同期したタイミングでスロットを合わせ、データ送信
    *   キャリアセンスで衝突回避
    *   正確なToAに基づいたチャネル占有時間のシミュレーション
