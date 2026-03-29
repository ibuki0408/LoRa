# Mathematical Theory for LoRaWAN Analysis

このディレクトリには、LoRaWANネットワークにおけるキャリアセンス（CS）およびパケット衝突の理論モデルを計算するためのスクリプト `theory_ppp.jl` が含まれています。以下に、使用されている主要な数学的モデルと詳細な確率式をまとめます。

## 1. キャリアセンス (CS/LBT) モデル

### 1.1 物理リンクモデル (Path Loss & Shadowing)
端末 $A$ から距離 $d$ [m] 離れた端末 $B$ への受信電力 $P_{rx}(d)$ [dBm] は、以下の対数距離パスロスモデルとシャドウイング成分 $S$ によって規定されます。

$$P_{rx}(d) = P_{tx} - \left( PL_{ref} + 10 \alpha \log_{10}(d) \right) - S$$

- $P_{tx}$ [dBm]: 送信電力 (Default: 13.0)
- $PL_{ref}$ [dB]: 基準パスロス (1m地点)
- $\alpha$: パスロス指数 (Default: 2.7)
- $S \sim \mathcal{N}(0, \sigma^2)$: シャドウイング (Default: $\sigma = 8.0$ dB)

### 1.2 実効センス確率 (Average Sensing Probability) $\eta_{avg}$
円形エリア（半径 $R_{max}$）内に一様に分布する2つのランダムな点 $A, B$ について、お互いにキャリアセンス閾値 $\gamma_{th}$ を超える電力を受信し、送信を抑制し合える平均確率 $\eta_{avg}$ を数値積分で求めます。

この確率は、点 $A$ の位置 $(r, 0)$ を固定した際の「境界内でのセンス領域の重なり」をエリア全体で平均化したものです。

$$\eta_{avg} = \iint_{Area} f(A) dA \iint_{Area} f(B) dB \cdot P(P_{rx}(\|A-B\|) > \gamma_{th})$$

離散形式による実装式（2次元数値積分）:
$$\eta_{avg} = \sum_{r} w_r \sum_{\rho} w_\rho \sum_{\theta} w_\theta \cdot Q\left( \frac{\gamma_{th} - \bar{P}_{rx}(dist(r, \rho, \theta))}{\sigma} \right)$$

- $dist(r, \rho, \theta) = \sqrt{r^2 + \rho^2 - 2r\rho\cos(\theta)}$
- $w_r, w_\rho$: 半径方向の重み ($2r/R^2$)
- $Q(x)$: 標準正規分布の右側累積確率関数

### 1.3 Slotted CS 競合モデル (Contention Factor $M$)
1つのタイムスロット・1つのチャネルにおいて、ある端末のセンス範囲内で同時に送信を試みる「ライバル端末数」の期待値を $M$ とします。

$$M = \frac{N \times \eta_{avg} \times T_{slot}}{T_{interval} \times N_{channels}}$$

- $N$: 全端末数
- $T_{slot}$: 競合ウィンドウ長 (Default: 100ms)
- $T_{interval}$: 平均送信周期 (Default: 30.0s)
- $N_{channels}$: 利用可能チャネル数

この $M$ 台による競合で、ジッタにより最小の遅延を引いた1台だけが送信でき、残りがブロックされると仮定すると、送信成功（LBT通過）確率 $P_{LBT}$ は以下になります。

$$P_{LBT} = \frac{1 - e^{-M}}{M}$$

---

## 2. 衝突・成功モデル (Poisson-Discrete Model)

### 2.1 実効負荷 $G$ (Offered Load)
解析対象の端末にとって衝突の脅威となる「自分を検知できない（＝隠れ端末）」かつ「現時点でアクティブな」端末の密度を $G_{hidden}$ とします。

$$G_{hidden} = \frac{N \times (1 - \eta_{avg}) \times T_{vuln}}{T_{interval} \times N_{channels}}$$

- 同期スロットの場合: $T_{vuln} = T_{slot} = 100$ms
- 非同期（Pure ALOHA）の場合: $T_{vuln} = 2 \times T_{airtime}$ [s]

### 2.2 キャプチャ効果 (Capture Effect) $P_{cap}(k)$
$k$ 台の干渉端末が存在する状況下で、希望信号の電力が全干渉電力の和よりも SIR 閾値以上に高い場合に受信成功とみなします。

$$P_{success} = \sum_{k=0}^{\infty} \frac{G^k e^{-G}}{k!} \cdot P_{cap}(k)$$

- $P_{cap}(0) = 1.0$ (干渉なし)
- $P_{cap}(k \ge 1)$: 実装では $k$ に依存する減衰モデル $\frac{C}{k^{1.3}}$ 等で近似

---

## 3. 指標の定義とシミュレータとの対応

### 3.1 Retry-aware PER (`retry_per`)
シミュレータの **Original PER** と比較すべき指標です。LBT（キャリアセンス）による一時的な送信待機（バックオフ）を「失敗」に含めず、隠れ端末との不可避な衝突のみをカウントします。
$$PER_{retry} = 1 - P_{success}(G_{hidden})$$

### 3.2 One-shot PER (`one_shot_per`)
再送を想定せず、LBTでのブロックも失敗とみなす、最も厳しい通信環境の指標です。
$$PER_{one\_shot} = 1 - (P_{LBT} \times P_{success}(G_{total\_active}))$$

### 3.3 CS Block Rate
端末が送信を試みた際、CS（キャリアセンス）により送信を控える確率です。
$$CS\_Block\_Rate = 1 - P_{LBT}$$
