"""
Simulation Part 2: Pure ALOHA vs Slotted ALOHA vs LBT (CSMA/CA)
LoRaWAN-like 920 MHz band environment

Output:
  - results/sim_part2_pdr.csv : PDR results per (protocol, N)
  - results/sim_part2_pdr.png : Comparison graph
"""

import numpy as np
import csv
import os

# ─────────────────────────────────────────────
#  0.  シミュレーションパラメータ（変更不要）
# ─────────────────────────────────────────────
SEED      = 2026
N_TRIALS  = 100
SIM_TIME  = 30 * 60          # 30分 [s]
N_LIST    = [100, 200, 300, 400, 500]

# True  : SINRが閾値以上なら受信成功
# False : 衝突したら全パケット失敗（理論式の値と一致するかを確認）
CAPTURE_EFFECT = True

R            = 1000.0         # セル半径 [m]
LAMBDA       = 1.0 / 30.0    # パケット生成レート [pkt/s]
PKT_DUR      = 0.070          # パケット長 [s] (70ms)
P_TX_DBM     = 13.0           # 送信電力 [dBm]
D0           = 1.0            # 基準距離 [m]
L0_DB        = 32.0           # 基準パスロス [dB]
D_BP         = 200.0          # ブレークポイント距離 [m]
SIGMA_SH     = 8.0            # シャドーイング標準偏差 [dB]
NF_DB        = 6.0            # ノイズフィギュア [dB]
BW_HZ        = 125e3          # 帯域幅 [Hz]
GAMMA_REQ_DB = 6.0            # SINR閾値 [dB]
N_DBM        = -174.0 + 10.0 * np.log10(BW_HZ) + NF_DB  # 熱雑音電力 [dBm]
SLOT_DUR     = 0.100          # スロット長 [s] (100ms)
CS_DUR       = 0.005          # CS観測時間 [s] (5ms)
P_CS_DBM     = -80.0          # CS閾値 [dBm]
CW_MAX       = 15
MAX_BACKOFF  = 3

#単位変換
def dbm_to_w(x_dbm):
    return 10.0 ** ((x_dbm - 30.0) / 10.0)

def db_to_lin(x_db):
    return 10.0 ** (x_db / 10.0)

GAMMA_REQ_LIN = db_to_lin(GAMMA_REQ_DB)
N_W           = dbm_to_w(N_DBM)
P_CS_W        = dbm_to_w(P_CS_DBM)

# ─────────────────────────────────────────────
#  1.  端末配置（セル内に一様ランダム）
# ─────────────────────────────────────────────
def place_devices(n, rng):
    """半径Rの円内に一様ランダムに配置。(x, y, d) を返す。"""
    r   = R * np.sqrt(rng.uniform(0, 1, n))
    phi = rng.uniform(0, 2 * np.pi, n)
    x   = r * np.cos(phi)
    y   = r * np.sin(phi)
    d   = np.maximum(np.sqrt(x**2 + y**2), 1.0)
    return x, y, d

# ─────────────────────────────────────────────
#  2.  伝搬モデル（2勾配パスロス + シャドーイング）
# ─────────────────────────────────────────────
def path_loss_db(d):
    """2勾配パスロスモデル [dB]。d<=D_BP: 20log, d>D_BP: 30log。"""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    return np.where(
        d <= D_BP,
        L0_DB + 20.0 * np.log10(d / D0),
        L0_DB + 20.0 * np.log10(D_BP / D0) + 30.0 * np.log10(d / D_BP)
    )

def rx_power_dbm(d, shadow_db):
    """受信電力 [dBm] = 送信電力 - パスロス + シャドーイング"""
    return P_TX_DBM - path_loss_db(d) + shadow_db

# ─────────────────────────────────────────────
#  3.  トラフィック生成（ポアソン過程）
# ─────────────────────────────────────────────
def generate_packets(n_dev, rng):
    """各端末のパケット到着時刻リストを返す（ポアソン分布）。"""
    packets = []
    for _ in range(n_dev):
        n_pkts = rng.poisson(LAMBDA * SIM_TIME)
        if n_pkts == 0:
            packets.append(np.array([]))
        else:
            packets.append(np.sort(rng.uniform(0, SIM_TIME, n_pkts)))
    return packets

# ─────────────────────────────────────────────
#  4a. Pure ALOHA シミュレーション
# ─────────────────────────────────────────────
def sim_pure_aloha(n_dev, rng, capture=None):
    """Pure ALOHA: チャネル確認なし、即時送信。"""
    cap      = CAPTURE_EFFECT if capture is None else capture
    _x, _y, d = place_devices(n_dev, rng)
    p_rx_w   = dbm_to_w(rx_power_dbm(d, rng.normal(0, SIGMA_SH, n_dev)))
    packets  = generate_packets(n_dev, rng)

    # 全送信イベントを (開始時刻, 終了時刻, 端末ID) のリストに展開
    events = [(t, t + PKT_DUR, i)
              for i, times in enumerate(packets) for t in times]
    if not events:
        return 0, 0

    events.sort()
    starts = np.array([e[0] for e in events])
    ends   = np.array([e[1] for e in events])
    devs   = np.array([e[2] for e in events])
    n_succ = 0

    for idx, (t_s, t_e, dev) in enumerate(events):
        # 干渉しうる区間を二分探索で絞り込む
        i0 = np.searchsorted(starts, t_s - PKT_DUR, side='left')
        i1 = np.searchsorted(starts, t_e,           side='right')
        ifc_pwr    = 0.0
        has_ifc    = False
        for j in range(i0, i1):
            if j == idx: continue
            if starts[j] < t_e and ends[j] > t_s:
                has_ifc = True
                ifc_pwr += p_rx_w[devs[j]]

        if cap:
            if p_rx_w[dev] / (ifc_pwr + N_W) >= GAMMA_REQ_LIN:
                n_succ += 1
        else:
            if not has_ifc:
                n_succ += 1

    return len(events), n_succ

# ─────────────────────────────────────────────
#  4b. Slotted ALOHA シミュレーション
# ─────────────────────────────────────────────
def sim_slotted_aloha(n_dev, rng, capture=None):
    """Slotted ALOHA: スロット境界でのみ送信。同スロット内の複数送信は衝突。"""
    cap      = CAPTURE_EFFECT if capture is None else capture
    _x, _y, d = place_devices(n_dev, rng)
    p_rx_w   = dbm_to_w(rx_power_dbm(d, rng.normal(0, SIGMA_SH, n_dev)))
    packets  = generate_packets(n_dev, rng)

    # スロットインデックスでグループ化
    slot_map = {}
    n_tx = 0
    for i, times in enumerate(packets):
        for t in times:
            slot_map.setdefault(int(np.ceil(t / SLOT_DUR)), []).append(i)
            n_tx += 1

    n_succ = 0
    for devs in slot_map.values():
        if len(devs) == 1:
            dev = devs[0]
            if cap:
                if p_rx_w[dev] / N_W >= GAMMA_REQ_LIN:
                    n_succ += 1
            else:
                n_succ += 1  # 衝突なし → 成功
        else:
            if cap:
                for dev in devs:
                    ifc = sum(p_rx_w[o] for o in devs if o != dev)
                    if p_rx_w[dev] / (ifc + N_W) >= GAMMA_REQ_LIN:
                        n_succ += 1
            # capture=False の場合: 衝突 → 全パケット失敗（何もしない）

    return n_tx, n_succ

# ─────────────────────────────────────────────
#  4c. LBT (CSMA/CA) シミュレーション
# ─────────────────────────────────────────────
def sim_lbt(n_dev, rng, capture=None):
    """
    LBT: 送信前に5msのCS（キャリアセンス）観測
    CS判定は端末間の伝搬損失を使用
    ヒープによる2段階イベント処理:
      ev_type=0 (arrival) : CS観測開始 → クリアなら5ms後にcs_endをスケジュール
      ev_type=1 (cs_end)  : 5ms後に再確認 → クリアなら送信開始
    SINRは全イベント終了後に tx_log から一括計算
    """
    import heapq

    cap = CAPTURE_EFFECT if capture is None else capture
    x_pos, y_pos, d = place_devices(n_dev, rng)
    p_rx_w = dbm_to_w(rx_power_dbm(d, rng.normal(0, SIGMA_SH, n_dev)))
    packets = generate_packets(n_dev, rng)

    # 端末ペアのシャドーイングを事前生成（固定行列）
    # sh_peer[i, j] = 端末j → 端末i のシャドーイング [dB]（シミュ中固定）
    sh_peer = rng.normal(0, SIGMA_SH, (n_dev, n_dev))

    # ヒープ要素: (時刻, ev_type, gen_t, 端末ID, バックオフ回数)
    heap = []
    n_tx_attempt = 0
    for i, times in enumerate(packets):
        for t in times:
            heapq.heappush(heap, (t, 0, t, i, 0))
            n_tx_attempt += 1

    ongoing = []  # (tx_end_time, tx_device_idx) 現在送信中のリスト
    tx_log  = []  # (tx_start, tx_end, dev) SINR計算用ログ

    def _sensed_power(sensing_dev):
        """sensing_dev から見た受信電力合計（端末間パスロス + 固定シャドーイング）。"""
        total = 0.0
        for (_, tx_dev) in ongoing:
            dx   = x_pos[sensing_dev] - x_pos[tx_dev]
            dy   = y_pos[sensing_dev] - y_pos[tx_dev]
            d_ij = max(np.sqrt(dx**2 + dy**2), 1.0)
            p_dbm = P_TX_DBM - path_loss_db(d_ij)[0] + sh_peer[sensing_dev, tx_dev]
            total += dbm_to_w(p_dbm)
        return total

    while heap:
        t, ev_type, gen_t, dev, bo_cnt = heapq.heappop(heap)
        ongoing = [(e, d_) for (e, d_) in ongoing if e > t]
        sp = _sensed_power(dev)

        if ev_type == 0:
            # CS観測開始: クリアなら5ms後にcs_endをスケジュール
            if sp < P_CS_W:
                heapq.heappush(heap, (t + CS_DUR, 1, gen_t, dev, bo_cnt))
            elif bo_cnt < MAX_BACKOFF:
                wait = rng.integers(0, CW_MAX + 1) * SLOT_DUR
                heapq.heappush(heap, (t + wait, 0, gen_t, dev, bo_cnt + 1))

        else:
            # CS観測終了: 5ms継続してクリアなら送信開始
            if sp < P_CS_W:
                tx_end = t + PKT_DUR
                ongoing.append((tx_end, dev))
                tx_log.append((t, tx_end, dev))
            elif bo_cnt < MAX_BACKOFF:
                wait = rng.integers(0, CW_MAX + 1) * SLOT_DUR
                heapq.heappush(heap, (t + wait, 0, gen_t, dev, bo_cnt + 1))

    # 全イベント終了後にSINRを一括計算 (Pure ALOHAと同じO(N log N)方式)
    if not tx_log:
        return n_tx_attempt, 0

    tx_log.sort()
    starts = np.array([e[0] for e in tx_log])
    ends   = np.array([e[1] for e in tx_log])
    devs   = np.array([e[2] for e in tx_log])
    n_succ = 0

    for idx, (t_s, t_e, dev) in enumerate(tx_log):
        i0 = np.searchsorted(starts, t_s - PKT_DUR, side='left')
        i1 = np.searchsorted(starts, t_e,           side='right')
        ifc_pwr = 0.0
        has_ifc = False
        for j in range(i0, i1):
            if j == idx: continue
            if starts[j] < t_e and ends[j] > t_s:
                has_ifc = True
                ifc_pwr += p_rx_w[devs[j]]

        if cap:
            if p_rx_w[dev] / (ifc_pwr + N_W) >= GAMMA_REQ_LIN:
                n_succ += 1
        else:
            if not has_ifc:
                n_succ += 1

    return n_tx_attempt, n_succ

# ─────────────────────────────────────────────
#  5.  モンテカルロ実行
# ─────────────────────────────────────────────
PROTOCOLS = {
    "Pure ALOHA":    sim_pure_aloha,
    "Slotted ALOHA": sim_slotted_aloha,
    "LBT":           sim_lbt,
}

def run_all():
    rng_master = np.random.default_rng(SEED)
    results = {p: {} for p in PROTOCOLS}

    for proto_name, sim_fn in PROTOCOLS.items():
        print(proto_name, flush=True)
        for N in N_LIST:
            pdrs = []
            for _ in range(N_TRIALS):
                seed = rng_master.integers(0, 2**31)
                n_tx, n_succ = sim_fn(N, np.random.default_rng(seed))
                pdrs.append(n_succ / n_tx if n_tx > 0 else 0.0)
            results[proto_name][N] = float(np.mean(pdrs))
            print(f"  N={N}  PDR={results[proto_name][N]:.4f}", flush=True)

    return results

# ─────────────────────────────────────────────
#  6.  結果出力: CSV
# ─────────────────────────────────────────────
def save_csv(results, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        header = [
            "N", "Offered_Traffic_G",
            "Pure ALOHA (Sim)", "Pure ALOHA (Theory)",
            "Slotted ALOHA (Sim)", "Slotted ALOHA (Theory)",
            "LBT (Sim)",
        ]
        writer.writerow(header)
        for N in N_LIST:
            G_pkt  = N * LAMBDA * PKT_DUR   # Pure ALOHA 用 G
            G_slot = N * LAMBDA * SLOT_DUR  # Slotted ALOHA 用 G
            row = [
                N,
                f"{G_pkt:.6f}",
                f"{results['Pure ALOHA'][N]:.6f}",
                f"{np.exp(-2 * G_pkt):.6f}",          # 理論値: e^{-2G}
                f"{results['Slotted ALOHA'][N]:.6f}",
                f"{np.exp(-G_slot):.6f}",              # 理論値: e^{-G_slot}
                f"{results['LBT'][N]:.6f}",
            ]
            writer.writerow(row)

# ─────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    from datetime import datetime
    results = run_all()

    base = os.path.join(os.path.dirname(__file__), "results")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    save_csv(results, os.path.join(base, f"sim_part2_pdr_{timestamp}.csv"))
    # save_plot(results, os.path.join(base, f"sim_part2_pdr_{timestamp}.png"))  # グラフ出力（必要時にコメント解除）
