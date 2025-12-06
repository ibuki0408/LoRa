# ===== SINR ベースの衝突判定モジュール =====

using Printf

"""
SF ごとの必要 SINR (dB)

# 出典
- M. Bor et al., "Do LoRa Low-Power Wide-Area Networks Scale?", MSWiM 2016
- Semtech SX1276 Datasheet, Rev. 5, 2016
"""
const REQUIRED_SNR_DB = Dict(
    7  => -7.5,
    8  => -10.0,
    9  => -12.5,
    10 => -15.0,
    11 => -17.5,
    12 => -20.0
)

# Co-channel SIR (Signal to Interference Ratio) Threshold (dB)
# Semtech SX1276 Datasheet: ~6 dB for all SFs
const REQUIRED_CO_CHANNEL_SIR_DB = 6.0

function calculate_db(signal_w, noise_w)
    if noise_w <= 0
        return Inf
    end
    return 10 * log10(signal_w / noise_w)
end

function detect_collisions_sinr(finished_tx, sf::Int, noise_power_dbm::Float64)
    required_snr_db = REQUIRED_SNR_DB[sf]
    required_sir_db = REQUIRED_CO_CHANNEL_SIR_DB
    
    success = 0
    collisions = 0
    
    noise_w = 10^(noise_power_dbm / 10) * 1e-3
    
    for i in 1:length(finished_tx)
        p1 = finished_tx[i]
        signal_w = 10^(p1.rx_power_at_gw / 10) * 1e-3
        
        # 1. SNR Check (Noise only)
        snr_db = calculate_db(signal_w, noise_w)
        if snr_db < required_snr_db
            p1.status = "Lost (SNR)"
            collisions += 1 # technically lost, but counting as failure
            continue
        end
        
        # 2. SIR Check (Interference)
        interference_w = 0.0
        for j in 1:length(finished_tx)
            if i == j continue end
            p2 = finished_tx[j]
            
            # 時間重複チェック (Any overlap)
            if max(p1.start_ms, p2.start_ms) < min(p1.end_ms, p2.end_ms)
                interference_w += 10^(p2.rx_power_at_gw / 10) * 1e-3
            end
        end
        
        if interference_w > 0
            sir_db = calculate_db(signal_w, interference_w)
            if sir_db < required_sir_db
                p1.status = "Collision (SIR)"
                collisions += 1
                continue
            end
        end
        
        # Both checks passed
        p1.status = "Success"
        success += 1
    end
    
    return (success, collisions)
end
