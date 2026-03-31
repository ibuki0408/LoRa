"""
Simulation Part 2: Pure ALOHA vs Slotted ALOHA vs LBT (CSMA/CA)
LoRaWAN-like 920 MHz band environment

Output:
  - results/sim_part2_pdr.csv   : PDR results per (protocol, N)
  - results/sim_part2_pdr.png   : Comparison graph
"""

import numpy as np
import csv
import os
import matplotlib
matplotlib.use("Agg")   # headless (no display required)
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────
#  0.  Simulation parameters
# ─────────────────────────────────────────────
SEED           = 2026
N_TRIALS       = 100          # Monte-Carlo trials
SIM_TIME       = 30 * 60      # 30 min [s]
N_LIST         = [100, 200, 300, 400, 500]

# ── Capture effect toggle ────────────────────────────────────────────────────
# True  : SINR-based success (realistic)
# False : any collision → all fail (matches classical ALOHA theory)
CAPTURE_EFFECT = True

# ── LBT carrier-sense mode ────────────────────────────────────────────────────
# True  : Device-to-device path loss (like Prop.jl) → hidden-node effects emerge
# False : BS-centric view ("神の視点") → all devices visible, no hidden nodes
LBT_PEER_SENSING = True

# Spatial
R             = 1000.0       # cell radius [m]

# Traffic
LAMBDA        = 1.0 / 30.0   # packet rate per device [pkt/s]
PKT_DUR       = 0.070        # packet duration [s]  (70 ms)

# PHY
P_TX_DBM      = 13.0         # transmit power [dBm]
D0            = 1.0          # reference distance [m]
L0_DB         = 32.0         # path-loss at d0 [dB]
D_BP          = 200.0        # breakpoint distance [m]
SIGMA_SH      = 8.0          # shadowing std [dB]
NF_DB         = 6.0          # noise figure [dB]
BW_HZ         = 125e3        # bandwidth [Hz]
GAMMA_REQ_DB  = 6.0          # SINR threshold [dB]

# Noise power [dBm]
N_DBM = -174.0 + 10.0 * np.log10(BW_HZ) + NF_DB  # ≈ -120.0 dBm

# MAC – Slotted ALOHA
SLOT_DUR      = 0.100        # 100 ms

# MAC – LBT
CS_DUR        = 0.005        # carrier-sense observation time [s]  (5 ms)
P_CS_DBM      = -80.0        # carrier-sense threshold [dBm]
CW_MAX        = 15           # contention window  (0 … CW_MAX)
MAX_BACKOFF   = 3            # max backoff attempts before drop

# ─────────────────────────────────────────────
#  Helpers: unit conversion
# ─────────────────────────────────────────────
def dbm_to_w(x_dbm):
    return 10.0 ** ((x_dbm - 30.0) / 10.0)

def db_to_lin(x_db):
    return 10.0 ** (x_db / 10.0)

GAMMA_REQ_LIN = db_to_lin(GAMMA_REQ_DB)
N_W           = dbm_to_w(N_DBM)
P_CS_W        = dbm_to_w(P_CS_DBM)

# ─────────────────────────────────────────────
#  1.  Spatial model
# ─────────────────────────────────────────────
def place_devices(n, rng):
    """Uniform random placement inside circle of radius R.
    Returns (x, y, d) where d = distance to BS.
    """
    r   = R * np.sqrt(rng.uniform(0, 1, n))
    phi = rng.uniform(0, 2 * np.pi, n)
    x   = r * np.cos(phi)
    y   = r * np.sin(phi)
    d   = np.sqrt(x**2 + y**2)               # distance to BS
    d   = np.maximum(d, 1.0)                  # avoid d < 1 m
    return x, y, d

# ─────────────────────────────────────────────
#  2.  Channel model
# ─────────────────────────────────────────────
def path_loss_db(d):
    """Two-slope path-loss model [dB]."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    L = np.where(
        d <= D_BP,
        L0_DB + 20.0 * np.log10(d / D0),
        L0_DB + 20.0 * np.log10(D_BP / D0) + 30.0 * np.log10(d / D_BP)
    )
    return L

def rx_power_dbm(d, shadow_db):
    """Received power [dBm] for a given distance and shadowing sample."""
    return P_TX_DBM - path_loss_db(d) + shadow_db

# ─────────────────────────────────────────────
#  3.  Traffic generation  (pure Poisson, event-based)
# ─────────────────────────────────────────────
def generate_packets(n_dev, rng):
    """
    Returns a list of length n_dev; each element is a sorted numpy array
    of packet-generation times [s] within [0, SIM_TIME].
    """
    packets = []
    for _ in range(n_dev):
        # Number of packets in SIM_TIME ~ Poisson(LAMBDA * SIM_TIME)
        n_pkts = rng.poisson(LAMBDA * SIM_TIME)
        if n_pkts == 0:
            packets.append(np.array([]))
        else:
            times = np.sort(rng.uniform(0, SIM_TIME, n_pkts))
            packets.append(times)
    return packets

# ─────────────────────────────────────────────
#  4.  Collision / SINR check
# ─────────────────────────────────────────────
def check_sinr(sig_power_w, interferers_w):
    """Return True if SINR >= threshold."""
    I   = np.sum(interferers_w)
    sinr = sig_power_w / (I + N_W)
    return sinr >= GAMMA_REQ_LIN

# ─────────────────────────────────────────────
#  5a. Pure ALOHA simulation
# ─────────────────────────────────────────────
def sim_pure_aloha(n_dev, rng, capture=None):
    """
    Pure ALOHA with O(N log N) interference calculation via searchsorted.
    capture: override CAPTURE_EFFECT global (None = use global).
    """
    cap = CAPTURE_EFFECT if capture is None else capture
    _x, _y, d = place_devices(n_dev, rng)
    shadows  = rng.normal(0, SIGMA_SH, n_dev)
    p_rx_dbm = rx_power_dbm(d, shadows)
    p_rx_w   = dbm_to_w(p_rx_dbm)

    packets  = generate_packets(n_dev, rng)

    events = []
    for i, times in enumerate(packets):
        for t in times:
            events.append((t, t + PKT_DUR, i))

    if len(events) == 0:
        return 0, 0

    events.sort(key=lambda e: e[0])
    n_tx   = len(events)
    n_succ = 0

    starts = np.array([e[0] for e in events])
    ends   = np.array([e[1] for e in events])
    devs   = np.array([e[2] for e in events])

    for idx, (t_s, t_e, dev) in enumerate(events):
        idx_start = np.searchsorted(starts, t_s - PKT_DUR, side='left')
        idx_end   = np.searchsorted(starts, t_e,           side='right')

        interferer_power = 0.0
        has_interferer   = False
        for jdx in range(idx_start, idx_end):
            if jdx == idx:
                continue
            if starts[jdx] < t_e and ends[jdx] > t_s:
                has_interferer = True
                interferer_power += p_rx_w[devs[jdx]]

        if not cap:
            # No capture: collision = fail
            if not has_interferer:
                n_succ += 1
        else:
            sinr = p_rx_w[dev] / (interferer_power + N_W)
            if sinr >= GAMMA_REQ_LIN:
                n_succ += 1

    return n_tx, n_succ

# ─────────────────────────────────────────────
#  5b. Slotted ALOHA simulation
# ─────────────────────────────────────────────
def sim_slotted_aloha(n_dev, rng, capture=None):
    """
    Slotted ALOHA: packet waits until the next slot boundary.
    capture: override CAPTURE_EFFECT global (None = use global).
    """
    cap = CAPTURE_EFFECT if capture is None else capture
    _x, _y, d = place_devices(n_dev, rng)
    shadows  = rng.normal(0, SIGMA_SH, n_dev)
    p_rx_dbm = rx_power_dbm(d, shadows)
    p_rx_w   = dbm_to_w(p_rx_dbm)

    packets  = generate_packets(n_dev, rng)

    slot_map = {}
    n_tx = 0
    for i, times in enumerate(packets):
        for t in times:
            slot_idx = int(np.ceil(t / SLOT_DUR))
            slot_map.setdefault(slot_idx, []).append(i)
            n_tx += 1

    n_succ = 0
    for slot_idx, devs in slot_map.items():
        if len(devs) == 1:
            dev = devs[0]
            if not cap:
                n_succ += 1          # no noise-only failure in ideal model
            else:
                sinr = p_rx_w[dev] / N_W
                if sinr >= GAMMA_REQ_LIN:
                    n_succ += 1
        else:
            if not cap:
                pass                 # collision → all fail
            else:
                for dev in devs:
                    ifc_power = sum(p_rx_w[o] for o in devs if o != dev)
                    sinr = p_rx_w[dev] / (ifc_power + N_W)
                    if sinr >= GAMMA_REQ_LIN:
                        n_succ += 1

    return n_tx, n_succ

# ─────────────────────────────────────────────
#  5c. LBT (CSMA/CA) simulation
# ─────────────────────────────────────────────
def sim_lbt(n_dev, rng, capture=None):
    """
    LBT (CSMA/CA) with two-stage heap events and deferred SINR evaluation.

    LBT_PEER_SENSING=True  : each device measures power received FROM nearby
                             transmitters using device-to-device path loss
                             (like Prop.jl). Hidden-node collisions can occur.
    LBT_PEER_SENSING=False : BS-centric 「神の視点」 – all transmissions visible.
    """
    cap        = CAPTURE_EFFECT if capture is None else capture
    peer_sense = LBT_PEER_SENSING

    x_pos, y_pos, d = place_devices(n_dev, rng)
    shadows  = rng.normal(0, SIGMA_SH, n_dev)
    p_rx_dbm = rx_power_dbm(d, shadows)        # device → BS received power
    p_rx_w   = dbm_to_w(p_rx_dbm)

    raw_packets = generate_packets(n_dev, rng)

    import heapq
    heap = []
    n_tx_attempt = 0

    for i, times in enumerate(raw_packets):
        for t in times:
            heapq.heappush(heap, (t, 0, t, i, 0))
            n_tx_attempt += 1

    # ongoing: (tx_end_time, tx_device_idx) – packets physically on air
    ongoing = []
    tx_log  = []   # (tx_start, tx_end, dev) – for deferred SINR

    def _sense(sensing_dev):
        """Total received power at sensing_dev from all ongoing transmitters."""
        if not peer_sense:
            # BS-centric: use BS-received power (no hidden-node effect)
            return sum(p_rx_w[tx_dev] for (_, tx_dev) in ongoing)
        # Peer-sensing: device-to-device path loss (like Prop.jl)
        total = 0.0
        for (_, tx_dev) in ongoing:
            dx   = x_pos[sensing_dev] - x_pos[tx_dev]
            dy   = y_pos[sensing_dev] - y_pos[tx_dev]
            d_ij = max(np.sqrt(dx*dx + dy*dy), 1.0)
            sh_ij      = rng.normal(0, SIGMA_SH)      # independent per pair
            p_peer_dbm = P_TX_DBM - path_loss_db(d_ij)[0] + sh_ij
            total += dbm_to_w(p_peer_dbm)
        return total

    while heap:
        t, ev_type, gen_t, dev, bo_cnt = heapq.heappop(heap)

        ongoing = [(e, d_) for (e, d_) in ongoing if e > t]
        sensed_power = _sense(dev)

        if ev_type == 0:
            # ── Observation start ─────────────────────────────────────────
            if sensed_power < P_CS_W:
                heapq.heappush(heap, (t + CS_DUR, 1, gen_t, dev, bo_cnt))
            else:
                if bo_cnt < MAX_BACKOFF:
                    wait = rng.integers(0, CW_MAX + 1) * SLOT_DUR
                    heapq.heappush(heap, (t + wait, 0, gen_t, dev, bo_cnt + 1))

        else:  # ev_type == 1
            # ── Observation end (5 ms later): final idle check ────────────
            if sensed_power < P_CS_W:
                tx_end = t + PKT_DUR
                ongoing.append((tx_end, dev))   # store dev index
                tx_log.append((t, tx_end, dev))
            else:
                if bo_cnt < MAX_BACKOFF:
                    wait = rng.integers(0, CW_MAX + 1) * SLOT_DUR
                    heapq.heappush(heap, (t + wait, 0, gen_t, dev, bo_cnt + 1))

    # ── Deferred SINR evaluation (same O(N log N) method as Pure ALOHA) ──
    if not tx_log:
        return n_tx_attempt, 0

    tx_log.sort(key=lambda x: x[0])
    starts = np.array([e[0] for e in tx_log])
    ends   = np.array([e[1] for e in tx_log])
    devs   = np.array([e[2] for e in tx_log])

    n_succ = 0
    for idx, (t_s, t_e, dev) in enumerate(tx_log):
        idx_start = np.searchsorted(starts, t_s - PKT_DUR, side='left')
        idx_end   = np.searchsorted(starts, t_e,           side='right')

        interferer_power = 0.0
        has_interferer   = False
        for jdx in range(idx_start, idx_end):
            if jdx == idx:
                continue
            if starts[jdx] < t_e and ends[jdx] > t_s:
                has_interferer = True
                interferer_power += p_rx_w[devs[jdx]]

        if not cap:
            if not has_interferer:
                n_succ += 1
        else:
            sinr = p_rx_w[dev] / (interferer_power + N_W)
            if sinr >= GAMMA_REQ_LIN:
                n_succ += 1

    return n_tx_attempt, n_succ

# ─────────────────────────────────────────────
#  6.  Monte-Carlo runner
# ─────────────────────────────────────────────
PROTOCOLS = {
    "Pure ALOHA":    sim_pure_aloha,
    "Slotted ALOHA": sim_slotted_aloha,
    "LBT":           sim_lbt,
}

def run_all():
    """Run Monte-Carlo for both capture=True and capture=False modes."""
    rng_master = np.random.default_rng(SEED)

    # results[capture_flag][protocol][N] = mean PDR
    results = {True: {p: {} for p in PROTOCOLS},
               False: {p: {} for p in PROTOCOLS}}

    for cap in [True, False]:
        label = "capture ON" if cap else "capture OFF (no-capture)"
        print(f"\n{'='*50}")
        print(f"  Mode: {label}")
        print(f"{'='*50}")
        for proto_name, sim_fn in PROTOCOLS.items():
            print(f"\n  --- {proto_name} ---")
            for N in N_LIST:
                pdrs = []
                for trial in range(N_TRIALS):
                    trial_seed = rng_master.integers(0, 2**31)
                    trial_rng  = np.random.default_rng(trial_seed)
                    n_tx, n_succ = sim_fn(N, trial_rng, capture=cap)
                    pdrs.append(n_succ / n_tx if n_tx > 0 else 0.0)

                mean_pdr = float(np.mean(pdrs))
                results[cap][proto_name][N] = mean_pdr
                G = N * LAMBDA * PKT_DUR
                print(f"    N={N:3d}  G={G:.4f}  PDR={mean_pdr:.4f}")

    return results

# ─────────────────────────────────────────────
#  7.  Output: CSV + PNG
# ─────────────────────────────────────────────
def save_csv(results, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        # Columns: N, G, then each protocol × each mode
        header = ["N", "Offered_Traffic_G"]
        for p in PROTOCOLS:
            header += [f"{p} (capture)", f"{p} (no-capture)"]
        writer.writerow(header)
        for N in N_LIST:
            G   = N * LAMBDA * PKT_DUR
            row = [N, f"{G:.6f}"]
            for p in PROTOCOLS:
                row += [f"{results[True][p][N]:.6f}",
                        f"{results[False][p][N]:.6f}"]
            writer.writerow(row)
    print(f"\nCSV saved → {path}")

def save_plot(results, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    G_list = [N * LAMBDA * PKT_DUR for N in N_LIST]

    # Protocol colour map
    colors = {
        "Pure ALOHA":    "#E74C3C",
        "Slotted ALOHA": "#3498DB",
        "LBT":           "#2ECC71",
    }
    markers = {"Pure ALOHA": "o", "Slotted ALOHA": "s", "LBT": "^"}

    # ── Simulation curves ─────────────────────────────────────────────────
    for proto_name in PROTOCOLS:
        c  = colors[proto_name]
        mk = markers[proto_name]
        # capture ON  → solid
        pdr_on  = [results[True][proto_name][N]  for N in N_LIST]
        # capture OFF → dashed, lighter
        pdr_off = [results[False][proto_name][N] for N in N_LIST]
        for ax, xs in zip(axes, [N_LIST, G_list]):
            ax.plot(xs, pdr_on,  color=c, marker=mk, linestyle="-",
                    linewidth=2, markersize=7,
                    label=f"{proto_name} (capture)")
            ax.plot(xs, pdr_off, color=c, marker=mk, linestyle="--",
                    linewidth=1.5, markersize=5, alpha=0.6,
                    label=f"{proto_name} (no-capture)")

    # ── Theoretical curves (dotted) ───────────────────────────────────────
    G_fine = np.linspace(1e-3, max(G_list) * 1.3, 400)
    N_fine = G_fine / (LAMBDA * PKT_DUR)
    th_kw  = dict(linewidth=1.6, linestyle=":", alpha=0.9)
    for ax in axes:
        x = N_fine if ax is axes[0] else G_fine
        ax.plot(x, np.exp(-2 * G_fine), color="#E74C3C",
                label="Pure ALOHA (theory)", **th_kw)
        ax.plot(x, np.exp(-G_fine),     color="#3498DB",
                label="Slotted ALOHA (theory)", **th_kw)

    # ── Axes formatting ───────────────────────────────────────────────────
    axes[0].set_xlabel("Number of Devices $N$", fontsize=12)
    axes[0].set_ylabel("PDR", fontsize=12)
    axes[0].set_title("PDR vs Number of Devices", fontsize=13, fontweight="bold")
    axes[0].legend(fontsize=8, ncol=2, loc="upper right")
    axes[0].grid(True, alpha=0.35)
    axes[0].set_ylim(-0.05, 1.05)

    axes[1].set_xlabel("Offered Traffic $G$", fontsize=12)
    axes[1].set_ylabel("PDR", fontsize=12)
    axes[1].set_title("PDR vs Offered Traffic", fontsize=13, fontweight="bold")
    axes[1].legend(fontsize=8, ncol=2, loc="upper right")
    axes[1].grid(True, alpha=0.35)
    axes[1].set_ylim(-0.05, 1.05)

    plt.suptitle(
        "Pure ALOHA vs Slotted ALOHA vs LBT  |  LoRaWAN-like 920 MHz\n"
        "solid = capture ON,  dashed = capture OFF (no-capture),  dotted = classical theory",
        fontsize=12, fontweight="bold", y=1.03)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Graph saved → {path}")

# ─────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    print("Starting simulation…")
    print(f"  Noise power  : {N_DBM:.2f} dBm")
    print(f"  SINR thresh  : {GAMMA_REQ_DB:.1f} dB")
    print(f"  CS threshold : {P_CS_DBM:.1f} dBm")
    print(f"  Trials       : {N_TRIALS}")
    print(f"  Sim time     : {SIM_TIME/60:.0f} min")

    results = run_all()

    base = os.path.join(os.path.dirname(__file__), "..", "results")
    save_csv(results, os.path.join(base, "sim_part2_pdr.csv"))
    save_plot(results, os.path.join(base, "sim_part2_pdr.png"))

    print("\nDone.")
