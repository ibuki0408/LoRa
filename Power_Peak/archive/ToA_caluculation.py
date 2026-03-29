import math

def calculate_lora_toa(sf, bw, pl, cr, preamble_len=8, ih=0, de=1):
    """
    LoRaの総送信時間 (ToA: Time on Air) を計算する関数。
    標準的なセンサーIoT環境 (SF>=11, BW=125kHz) を想定し、DE=1をデフォルト値としています。

    Args:
        sf (int): Spreading Factor (拡散率, 6-12).
        bw (float): Bandwidth (帯域幅, Hz). 例: 125000.0 (125kHz).
        pl (int): Payload Length (ペイロードバイト数).
        cr (int): Coding Rate (符号化率パラメータ, 1-4). 実際の符号化率は 4/(4+CR).
        preamble_len (int): Preamble Length (プリアンブルシンボル数). デフォルトはLoRaWAN標準の8.
        ih (int): Implicit Header Mode (暗黙的ヘッダ, 0=明示的/標準, 1=暗黙的).
        de (int): Low Data Rate Optimization (低データレート最適化, 0または1).
        
    Returns:
        float: 総送信時間 ToA (秒).
    """
    
    # 1. 1シンボルあたりの時間 (T_symb) の計算
    t_symb = (2**sf) / bw
    
    # 2. プリアンブル時間 (T_preamble) の計算
    t_preamble = (preamble_len + 4.25) * t_symb
    
    # 3. ペイロードシンボル数 (N_payload) の計算
    
    # 複雑な式の分子
    numerator = 8 * pl - 4 * sf + 28 + 16 * cr - 20 * ih
    
    # 複雑な式の分母
    denominator = 4 * (sf - 2 * de)
    
    # 天井関数内の計算結果 (切り上げ対象)
    ceil_input = numerator / denominator
    
    # N_payload の主要部
    # max(0, ...) で負の値を防ぎ、ceil(input) で切り上げ、(CR+4)を乗算
    n_payload_main = math.ceil(max(0, ceil_input)) * (cr + 4)
    
    n_payload = 8 + n_payload_main
    
    # 4. ペイロード時間 (T_payload) の計算
    t_payload = n_payload * t_symb
    
    # 5. 総送信時間 (ToA) の計算
    toa = t_preamble + t_payload
    
    return toa

# --- 使用例：最悪ケースの計算 ---
# 設定: SF=12, BW=125kHz, PL=51バイト, CR=4, IH=0, DE=1
SF_MAX = 10
BW_NARROW = 125000.0 # 125 kHz
PL_MAX = 10
CR_MAX = 1
IH_STANDARD = 0
DE_OPTIMIZED = 0 # SF>=11 and BW=125kHz のため

toa_max = calculate_lora_toa(
    sf=SF_MAX, 
    bw=BW_NARROW, 
    pl=PL_MAX, 
    cr=CR_MAX, 
    ih=IH_STANDARD, 
    de=DE_OPTIMIZED
)

# --- 使用例：最短ケースの計算 ---
# 設定: SF=7, BW=500kHz, PL=1バイト, CR=1, IH=0, DE=0 (SF<11のため)
SF_MIN = 7
BW_WIDE = 125000.0 # 500 kHz
PL_MIN = 10
CR_MIN = 4
DE_NORMAL = 0 # SF<11のため

toa_min = calculate_lora_toa(
    sf=SF_MIN, 
    bw=BW_WIDE, 
    pl=PL_MIN, 
    cr=CR_MIN,
    ih=IH_STANDARD,
    de=DE_NORMAL
)


print("--- LoRa ToA 計算結果 ---")
print(f"最悪ケース (SF={SF_MAX}, PL={PL_MAX}B, CR={CR_MAX}) の ToA: {toa_max:.4f} 秒")
print(f"最短ケース (SF={SF_MIN}, PL={PL_MIN}B, CR={CR_MIN}) の ToA: {toa_min:.4f} 秒")