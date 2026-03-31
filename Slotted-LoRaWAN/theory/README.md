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
セル内にランダムに配置された2端末が、お互いを「検知できる（＝LBTで避け合える）」平均的な確率です。単純な距離だけでなく、シャドウイングによる不規則な変動を含めてエリア全体で平均化しています（**空間的結合係数**とも呼べます）。

$$\eta_{avg} = \frac{1}{A^2} \iint_{Area} \iint_{Area} P(P_{rx}(\|A-B\|) > \gamma_{th}) \, dA \, dB$$

> [!TIP]
> **直感的な意味**: セル内の全端末のうち、自分が平均して何％の送信を「事前に耳を澄ませて回避できるか」を示します。$\eta_{avg} = 0.1$ なら、周囲の10%の端末とは互いに譲り合えますが、残りの90%は「隠れ端末」となり衝突の危険があります。

### 1.3 競合数 $M$ (Expected Contenders)
ある端末が送信しようとした際、そのセンス範囲内で「同じスロット・同じチャネル」を使って送信を試みているライバル端末の平均的な数です。

$$M = \frac{N \times \eta_{avg} \times T_{slot}}{T_{interval} \times N_{channels}}$$

- **$N \times \eta_{avg}$**: 自分のセンス範囲内に存在する平均端末数
- **$T_{slot} / T_{interval}$**: 1つのスロットに送信が重なる確率

### 1.4 LBT通過確率 (Transmission Probability) $P_{LBT}$
競合する $M$ 台のライバルのうち、自分が送信権を得られる確率です。

$$P_{LBT} = \frac{1 - e^{-M}}{M}$$

> [!IMPORTANT]
> **数式の由来（勝者1名モデル）**
> 1. 自分を含めて $k+1$ 台の端末が同時に送信を試みるとき、ジッタ（微小な待機時間）の差で「最初にチャネルが空いていると判定した1台」だけが送信でき、残りの $k$ 台はブロックされると仮定します。
> 2. このとき、自分がその「勝者の1台」に選ばれる確率は $1/(k+1)$ です。
> 3. ライバル数 $k$ がポアソン分布 $\text{Poi}(M)$ に従うとき、その期待値は $\sum_{k=0}^{\infty} \frac{1}{k+1} \frac{M^k e^{-M}}{k!}$ となり、これを計算すると上記の $\frac{1 - e^{-M}}{M}$ が導かれます。


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
