"""
Simulation Part 2: Pure ALOHA vs Slotted ALOHA vs LBT (CSMA/CA)
LoRaWAN-like 920 MHz band environment

【課題】このコードには TODO が複数あります。
      各 TODO を埋めてシミュレーションを完成させてください。
      完成したら実行して results/sim_part2_pdr.csv を確認してください。
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

def dbm_to_w(x_dbm):
    # TODO: dBm を W（ワット）に変換
    pass

def db_to_lin(x_db):
    return 10.0 ** (x_db / 10.0)

GAMMA_REQ_LIN = db_to_lin(GAMMA_REQ_DB)
N_W           = dbm_to_w(N_DBM)
P_CS_W        = dbm_to_w(P_CS_DBM)

# ─────────────────────────────────────────────
#  1.  端末配置
# ─────────────────────────────────────────────
def place_devices(n, rng):
    """半径Rの円内に n 台の端末を一様ランダムに配置。戻り値: (x, y, d)"""
    # TODO: r（動径）と phi（方位角）を生成せよ
    r   = None
    phi = None
    # TODO: (r, phi) から直交座標 (x, y) を計算せよ
    x   = None
    y   = None
    # TODO: BSからの距離 d を計算せよ（最小値 1.0 m を保証）
    d   = None
    return x, y, d

# ─────────────────────────────────────────────
#  2.  伝搬モデル
# ─────────────────────────────────────────────
def path_loss_db(d):
    """2勾配パスロスモデル [dB]"""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    return np.where(
        d <= D_BP,
        # TODO: d <= D_BP のときのパスロス式を書け
        0.0,
        # TODO: d >  D_BP のときのパスロス式を書け
        0.0,
    )

def rx_power_dbm(d, shadow_db):
    """受信電力 [dBm] = 送信電力 - パスロス + シャドーイング"""
    # TODO: 受信電力を計算して return せよ
    pass

# ─────────────────────────────────────────────
#  3.  トラフィック生成（変更不要）
# ─────────────────────────────────────────────
def generate_packets(n_dev, rng):
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
        # 干渉しうる区間を二分探索で絞り込む（変更不要）
        i0 = np.searchsorted(starts, t_s - PKT_DUR, side='left')
        i1 = np.searchsorted(starts, t_e,           side='right')

        ifc_pwr = 0.0
        has_ifc = False
        for j in range(i0, i1):
            if j == idx:
                continue
            # TODO: 送信 j と idx が時間的に重なるか判定し、
            #       重なる場合は has_ifc = True にして ifc_pwr に加算

        # TODO: 成功判定を書け
        #   cap=True  : SINR = p_rx_w[dev] / (ifc_pwr + N_W) >= GAMMA_REQ_LIN なら成功
        #   cap=False : has_ifc が False のときのみ成功

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

    # スロットインデックスでグループ化（変更不要）
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
            # TODO: 単独送信の成功判定を書け
            #   cap=True  : SINR = p_rx_w[dev] / N_W >= GAMMA_REQ_LIN なら成功
            #   cap=False : 衝突なしなので必ず成功
        else:
            # TODO: 衝突時の処理を書け
            #   cap=True  : 各端末について他の端末の受信電力を干渉として SINR 判定
            #               ifc = 同スロット内の自分以外の p_rx_w の合計
            #   cap=False : 全パケット失敗（何もしなくてよい）
            pass

    return n_tx, n_succ

# ─────────────────────────────────────────────
#  4c. LBT (CSMA/CA) シミュレーション
# ─────────────────────────────────────────────
def sim_lbt(n_dev, rng, capture=None):
    """
    LBT: 送信前に CS_DUR=5ms のキャリアセンス観測を行う。
    CS判定は端末間の伝搬損失を使用。

    ヒープの要素: (時刻, ev_type, gen_t, 端末ID, バックオフ回数)
      ev_type=0 (arrival) : CS観測開始
      ev_type=1 (cs_end)  : CS観測終了 → クリアなら送信

    SINR は全イベント終了後に tx_log から一括計算（変更不要）。
    """
    import heapq

    cap = CAPTURE_EFFECT if capture is None else capture
    x_pos, y_pos, d = place_devices(n_dev, rng)
    p_rx_w  = dbm_to_w(rx_power_dbm(d, rng.normal(0, SIGMA_SH, n_dev)))
    packets = generate_packets(n_dev, rng)

    # 端末ペアのシャドーイングを事前生成（固定行列）
    # sh_peer[i, j] = 端末j → 端末i のシャドーイング [dB]（シミュ中固定）
    sh_peer = rng.normal(0, SIGMA_SH, (n_dev, n_dev))

    heap = []
    n_tx_attempt = 0
    for i, times in enumerate(packets):
        for t in times:
            heapq.heappush(heap, (t, 0, t, i, 0))
            n_tx_attempt += 1

    ongoing = []  # (tx_end_time, tx_device_idx)
    tx_log  = []  # (tx_start, tx_end, dev)

    def _sensed_power(sensing_dev):
        """sensing_dev から見たチャネルの受信電力合計（端末間パスロス + 固定シャドーイング）。"""
        total = 0.0
        for (_, tx_dev) in ongoing:
            dx    = x_pos[sensing_dev] - x_pos[tx_dev]
            dy    = y_pos[sensing_dev] - y_pos[tx_dev]
            d_ij  = max(np.sqrt(dx**2 + dy**2), 1.0)
            p_dbm = P_TX_DBM - path_loss_db(d_ij)[0] + sh_peer[sensing_dev, tx_dev]
            total += dbm_to_w(p_dbm)
        return total

    while heap:
        t, ev_type, gen_t, dev, bo_cnt = heapq.heappop(heap)
        ongoing = [(e, d_) for (e, d_) in ongoing if e > t]
        sp = _sensed_power(dev)

        if ev_type == 0:
            # TODO: CS観測開始の処理を書け
            #   sp < P_CS_W（クリア）なら:
            #     → (t + CS_DUR) に ev_type=1 のイベントをヒープへ push
            #   sp >= P_CS_W（ビジー）かつ bo_cnt < MAX_BACKOFF なら:
            #     → wait = rng.integers(0, CW_MAX + 1) * SLOT_DUR だけ待って
            #       ev_type=0 のイベントをヒープへ push（bo_cnt + 1）
            pass

        else:  # ev_type == 1
            # TODO: CS観測終了の処理を書け
            #   sp < P_CS_W（まだクリア）なら:
            #     → 送信開始: tx_end = t + PKT_DUR
            #     → ongoing と tx_log に追加
            #   sp >= P_CS_W（観測中に誰かが送信開始）かつ bo_cnt < MAX_BACKOFF なら:
            #     → バックオフして ev_type=0 を push（bo_cnt + 1）
            pass

    # 全イベント終了後にSINRを一括計算（変更不要）
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
#  5.  モンテカルロ実行（変更不要）
# ─────────────────────────────────────────────
PROTOCOLS = {
    "Pure ALOHA (LBTなし)": sim_pure_aloha,
    "Pure ALOHA (LBTあり)": sim_lbt,
    "Slotted ALOHA":        sim_slotted_aloha,
}

def run_all():
    rng_master = np.random.default_rng(SEED)
    results = {p: {} for p in PROTOCOLS}

    for proto_name, sim_fn in PROTOCOLS.items():
        for N in N_LIST:
            pdrs = []
            for _ in range(N_TRIALS):
                seed = rng_master.integers(0, 2**31)
                n_tx, n_succ = sim_fn(N, np.random.default_rng(seed))
                pdrs.append(n_succ / n_tx if n_tx > 0 else 0.0)
            results[proto_name][N] = float(np.mean(pdrs))

    return results


# ─────────────────────────────────────────────
#  6.  CSV出力
# ─────────────────────────────────────────────
def save_csv(results, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "N", "Offered_Traffic_G",
            "Pure ALOHA (LBTなし) (Sim)", "Pure ALOHA Theory",
            "Pure ALOHA (LBTあり) (Sim)",
            "Slotted ALOHA (Sim)", "Slotted ALOHA Theory",
        ])
        for N in N_LIST:
            G_pkt  = N * LAMBDA * PKT_DUR
            G_slot = N * LAMBDA * SLOT_DUR
            row = [
                N,
                f"{G_pkt:.6f}",
                f"{results['Pure ALOHA (LBTなし)'][N]:.6f}",
                # TODO: Pure ALOHA の理論PDR を書け
                f"{0.0:.6f}",
                f"{results['Pure ALOHA (LBTあり)'][N]:.6f}",
                f"{results['Slotted ALOHA'][N]:.6f}",
                # TODO: Slotted ALOHA の理論PDR を書け
                f"{0.0:.6f}",
            ]
            writer.writerow(row)
    print(f"\nCSV saved → {path}")

# ─────────────────────────────────────────────
#  Main（変更不要）
# ─────────────────────────────────────────────
if __name__ == "__main__":
    print("シミュレーション開始")
    print(f"  熱雑音電力  : {N_DBM:.2f} dBm")
    print(f"  SINR閾値    : {GAMMA_REQ_DB:.1f} dB")
    print(f"  CS閾値      : {P_CS_DBM:.1f} dBm")
    print(f"  試行回数    : {N_TRIALS}")
    print(f"  シミュ時間  : {SIM_TIME/60:.0f} 分")

    results = run_all()

    base = os.path.join(os.path.dirname(__file__), "results")
    save_csv(results, os.path.join(base, "sim_part2_pdr.csv"))

    print("\n完了。")
