import numpy as np
import matplotlib.pyplot as plt
import random
from dataclasses import dataclass
from typing import List, Optional

# --- シミュレーション諸元 (Parameters) ---
@dataclass
class SimulationConfig:
    R: float = 1000.0  # セル半径 [m]
    Ptx_dbm: float = 13.0  # 送信電力 [dBm]
    packet_len_bytes: int = 1500 # 今回は100msを固定で使用 (プロパティ参照)
    d0: float = 1.0  # パスロス基準距離 [m]
    L0: float = 32.45 # 自由空間パスロス換算 (2.4GHz帯など想定, d=1m)
    dbp: float = 200.0 # ブレークポイント距離 [m] (この距離以降で減衰が急激になる)
    sigma_sh: float = 8.0 # シャドーイングの標準偏差 [dB]
    B_hz: float = 125e3 # 帯域幅 [Hz] (LoRaなどの想定)
    NF_db: float = 6.0 # ノイズフィギュア [dB]
    Gamma_req_db: float = 6.0 # 受信成功に必要なSINR閾値 [dB]
    Pcs_dbm: float = -80.0 # キャリアセンス閾値 [dBm] (LBTで使用)
    CW: int = 15 # コンテンションウィンドウサイズ (LBTで使用)
    max_retrans: int = 3 # 最大再送回数 (LBTで使用)
    use_capture_effect: bool = True # 捕捉効果 (Capture Effect) を考慮するか
    
    @property
    def packet_time_s(self):
        # パケット長を時間[s]で定義 (課題の想定に合わせて100ms)
        return 0.1 
    
    @property
    def noise_pwr_dbm(self):
        # 熱雑音電力 = -174 dBm/Hz + 10*log10(B) + NF
        return -174 + 10 * np.log10(self.B_hz) + self.NF_db

config = SimulationConfig()

# --- 伝搬モデル (Propagation Model) ---
def calculate_pathloss(d: float):
    """
    2勾配パスロスモデル (Two-slope model) の計算
    d <= dbp: 自由空間に近い減衰 (20 log10 d)
    d > dbp: 地面反射などの影響で急峻な減衰 (40 log10 d)
    """
    if d <= config.dbp:
        L = config.L0 + 20 * np.log10(max(d, config.d0) / config.d0)
    else:
        # ブレークポイント以降は勾配が40になるように接続
        L = config.L0 + 20 * np.log10(config.dbp / config.d0) + 40 * np.log10(d / config.dbp)
    return L

def get_received_power(dist: float, shadowing: float):
    # 受信電力 [dBm] = 送信電力 - パスロス + シャドーイング
    return config.Ptx_dbm - calculate_pathloss(dist) + shadowing

# --- 端末クラス (Terminal Helper) ---
class Terminal:
    def __init__(self, id: int):
        self.id = id
        # 円形セル内に一様にランダム配置するための極座標変換
        theta = random.uniform(0, 2 * np.pi)
        r = config.R * np.sqrt(random.uniform(0, 1)) # sqrtを使うことで面積一様になる
        self.dist = max(config.d0, r)
        self.shadowing = np.random.normal(0, config.sigma_sh)
        self.rx_power = get_received_power(self.dist, self.shadowing)

# --- MAC Protocols ---

def simulate_pure_aloha(N_terminals: int, lam_per_terminal: float, sim_time: float):
    terminals = [Terminal(i) for i in range(N_terminals)]
    
    events = []
    for t in terminals:
        curr_time = random.expovariate(lam_per_terminal)
        while curr_time < sim_time:
            events.append({
                'start_time': curr_time,
                'end_time': curr_time + config.packet_time_s,
                'rx_power': t.rx_power,
                'id': t.id
            })
            curr_time += random.expovariate(lam_per_terminal)
    
    events.sort(key=lambda x: x['start_time'])
    
    success_count = 0
    total_count = len(events)
    if total_count == 0: return 1.0
    
    # Efficient SINR check: only check nearby packets
    starts = np.array([e['start_time'] for e in events])
    ends = np.array([e['end_time'] for e in events])
    powers = 10**(np.array([e['rx_power'] for e in events]) / 10)
    noise_pwr_linear = 10**(config.noise_pwr_dbm/10)
    threshold_linear = 10**(config.Gamma_req_db / 10)
    
    for i in range(total_count):
        e_start = starts[i]
        e_end = ends[i]
        e_power = powers[i]
        
        # Range of potential overlaps: [e_start - T, e_end]
        idx_start = np.searchsorted(starts, e_start - config.packet_time_s, side='left')
        idx_end = np.searchsorted(starts, e_end, side='right')
        
        # Calculate interference from overlapping packets
        # indices to check: [idx_start, idx_end)
        interference_linear = 0
        collision_detected = False
        for j in range(idx_start, idx_end):
            if i == j: continue
            # Overlap check (necessary because of side windows)
            if not (ends[j] <= e_start or starts[j] >= e_end):
                interference_linear += powers[j]
                collision_detected = True
        
        # Success check
        if config.use_capture_effect:
            # Reality (Capture Effect): Success if SINR is enough
            if e_power / (interference_linear + noise_pwr_linear) >= threshold_linear:
                success_count += 1
        else:
            # Ideal (Perfect Collision Model): Success only if NO collision and power >= noise
            if not collision_detected and (e_power / noise_pwr_linear >= threshold_linear):
                success_count += 1
            
    return success_count / total_count

def simulate_slotted_aloha(N_terminals: int, lam_per_terminal: float, sim_time: float):
    slot_time = config.packet_time_s
    n_slots = int(sim_time / slot_time)
    
    success_count = 0
    total_sent = 0
    
    terminals = [Terminal(i) for i in range(N_terminals)]
    prob_tx = 1 - np.exp(-lam_per_terminal * slot_time)
    
    powers = np.array([t.rx_power for t in terminals])
    powers_linear = 10**(powers / 10)
    noise_pwr_linear = 10**(config.noise_pwr_dbm/10)
    threshold_linear = 10**(config.Gamma_req_db / 10)
    
    for s in range(n_slots):
        # Determine which terminals transmit in this slot
        active_mask = np.random.random(N_terminals) < prob_tx
        active_indices = np.where(active_mask)[0]
        n_active = len(active_indices)
        
        if n_active == 0:
            continue
            
        total_sent += n_active
        if n_active == 1:
            # Single transmitter
            idx = active_indices[0]
            if powers_linear[idx] / noise_pwr_linear >= threshold_linear:
                success_count += 1
        else:
            # Multiple transmitters (Collision)
            if config.use_capture_effect:
                # Reality (Capture Effect): Check SINR for each
                sum_pwr = np.sum(powers_linear[active_indices])
                for idx in active_indices:
                    interference = sum_pwr - powers_linear[idx]
                    if powers_linear[idx] / (interference + noise_pwr_linear) >= threshold_linear:
                        success_count += 1
            else:
                # Ideal: All packets in collision fail
                pass
                    
    return success_count / total_sent if total_sent > 0 else 1.0

def simulate_lbt_csma(N_terminals: int, lam_per_terminal: float, sim_time: float):
    terminals = [Terminal(i) for i in range(N_terminals)]
    
    # List to store all ACTUAL transmissions: (start, end, rx_power, tid)
    tx_log = []
    
    import heapq
    pq = []
    # (time, type, terminal_id, retrans_count)
    for i in range(N_terminals):
        heapq.heappush(pq, (random.expovariate(lam_per_terminal), 'arrival', i, 0))
    
    # Track current active transmissions for CS
    active_txs = [] # (end_time, rx_power)
    
    total_packets_arrived = 0
    
    while pq:
        time, type, tid, rc = heapq.heappop(pq)
        if time > sim_time: break
        
        if type == 'arrival':
            if rc == 0: total_packets_arrived += 1
            
            # Carrier Sense: clean up expired ones
            active_txs = [tx for tx in active_txs if tx[0] > time]
            
            current_energy_linear = sum([10**(tx[1]/10) for tx in active_txs])
            
            if 10 * np.log10(current_energy_linear + 1e-20) < config.Pcs_dbm:
                # Channel Clear -> Transmit
                end_time = time + config.packet_time_s
                tx_log.append({
                    'start_time': time,
                    'end_time': end_time,
                    'rx_power': terminals[tid].rx_power,
                    'id': tid
                })
                active_txs.append((end_time, terminals[tid].rx_power))
                # Schedule next packet arrival after this one ends (non-slotted non-persistent)
                heapq.heappush(pq, (end_time + random.expovariate(lam_per_terminal), 'arrival', tid, 0))
            else:
                # Channel Busy -> Backoff
                if rc < config.max_retrans:
                    backoff_slots = random.randint(0, config.CW)
                    # Use a realistic backoff slot (e.g., 20ms or similar, let's stick with 1ms or slightly larger)
                    # Requirement says 100ms packet, let's use 5ms backoff slot
                    backoff_time = backoff_slots * 0.005 
                    heapq.heappush(pq, (time + backoff_time, 'arrival', tid, rc + 1))
                else:
                    # Drop packet, schedule next arrival
                    heapq.heappush(pq, (time + random.expovariate(lam_per_terminal), 'arrival', tid, 0))

    if not tx_log: return 0.0
    
    # SINR check for all transmissions in log (same efficient method as Pure Aloha)
    tx_log.sort(key=lambda x: x['start_time'])
    starts = np.array([e['start_time'] for e in tx_log])
    ends = np.array([e['end_time'] for e in tx_log])
    powers = 10**(np.array([e['rx_power'] for e in tx_log]) / 10)
    noise_pwr_linear = 10**(config.noise_pwr_dbm/10)
    threshold_linear = 10**(config.Gamma_req_db / 10)
    
    success_count = 0
    for i in range(len(tx_log)):
        e_start = starts[i]
        e_end = ends[i]
        
        idx_start = np.searchsorted(starts, e_start - config.packet_time_s, side='left')
        idx_end = np.searchsorted(starts, e_end, side='right')
        
        collision_detected = False
        interference = 0
        for j in range(idx_start, idx_end):
            if i == j: continue
            if not (ends[j] <= e_start or starts[j] >= e_end):
                interference += powers[j]
                collision_detected = True
        
        if config.use_capture_effect:
            if powers[i] / (interference + noise_pwr_linear) >= threshold_linear:
                success_count += 1
        else:
            if not collision_detected and (powers[i] / noise_pwr_linear >= threshold_linear):
                success_count += 1
            
    return success_count / total_packets_arrived if total_packets_arrived > 0 else 1.0

# --- Evaluation Loop ---
def run_simulation():
    # SET SEED
    np.random.seed(2026)
    random.seed(2026)
    
    N_list = [100, 200, 300, 400, 500]
    lam_per_terminal = 1/30
    sim_time = 1800.0 # 30 minutes
    trials = 100
    
    results_pure = []
    results_slotted = []
    results_lbt = []
    
    print(f"Starting simulation: trials={trials}, sim_time={sim_time}s")
    
    # --- Capture Effect Toggle ---
    config.use_capture_effect = False # Set to False for the requested test
    print(f"Capture Effect: {'Enabled' if config.use_capture_effect else 'Disabled'}")
    
    import time as timer_lib
    start_sim = timer_lib.time()
    
    for N in N_list:
        p_pdr, s_pdr, l_pdr = [], [], []
        print(f"N={N}...", end="", flush=True)
        for t in range(trials):
            p_pdr.append(simulate_pure_aloha(N, lam_per_terminal, sim_time))
            s_pdr.append(simulate_slotted_aloha(N, lam_per_terminal, sim_time))
            l_pdr.append(simulate_lbt_csma(N, lam_per_terminal, sim_time))
            if (t+1) % 20 == 0: print(".", end="", flush=True)
            
        results_pure.append(np.mean(p_pdr))
        results_slotted.append(np.mean(s_pdr))
        results_lbt.append(np.mean(l_pdr))
        print(" Done")
        
    end_sim = timer_lib.time()
    print(f"Simulation completed in {end_sim - start_sim:.2f} seconds.")
        
    # Theoretical Values
    G_list = np.array(N_list) * lam_per_terminal * config.packet_time_s
    PDR_pure_theory = np.exp(-2 * G_list)
    PDR_slotted_theory = np.exp(-G_list)
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(N_list, results_pure, 'o-', label='Pure ALOHA (Sim)')
    plt.plot(N_list, PDR_pure_theory, '--', label='Pure ALOHA (Theory)', alpha=0.7)
    plt.plot(N_list, results_slotted, 's-', label='Slotted ALOHA (Sim)')
    plt.plot(N_list, PDR_slotted_theory, '--', label='Slotted ALOHA (Theory)', alpha=0.7)
    plt.plot(N_list, results_lbt, '^-', label='LBT (Sim)')
    
    plt.xlabel('Number of Terminals N')
    plt.ylabel('PDR')
    plt.title(f'PDR Comparison (Trials={trials}, Time={sim_time}s)')
    plt.legend()
    plt.grid(True)
    # 保存パスをカレントディレクトリ配下に変更
    plt.savefig('pdr_comparison_final.png')
    plt.show()

if __name__ == "__main__":
    run_simulation()
