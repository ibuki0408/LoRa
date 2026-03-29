# ===== collision_detection.jl =====
# Provides collision detection functions for integrated_lora_sim.jl

using Statistics

"""
    detect_collisions_simple_overlap(finished_tx)

Detect collisions based on simple time overlap in the same channel.
Updates the `status` field of TransmissionRecord objects to "Collision" if overlapped.
Returns (success_count, collision_count).
"""
function detect_collisions_simple_overlap(finished_tx)
    # Reset status first (assume Success initially, verify collision)
    # However, records come in as "Success" from MAC phase.
    
    n = length(finished_tx)
    if n == 0
        return (0, 0)
    end

    # Time sort for efficiency
    sort!(finished_tx, by = x -> x.start_ms)
    
    collision_indices = Set{Int}()
    
    for i in 1:n
        tx1 = finished_tx[i]
        
        # Check against subsequent packets
        for j in (i+1):n
            tx2 = finished_tx[j]
            
            # Optimization: If tx2 starts after tx1 ends, no overlap with tx1,
            # and since sorted, no further packets will overlap tx1.
            if tx2.start_ms >= tx1.end_ms
                break
            end
            
            # Check channel
            if tx1.channel != tx2.channel
                continue
            end
            
            # Overlap detected
            push!(collision_indices, i)
            push!(collision_indices, j)
        end
    end
    
    # Update Status
    for idx in collision_indices
        finished_tx[idx].status = "Collision"
    end
    
    collision_count = length(collision_indices)
    success_count = n - collision_count
    
    return (success_count, collision_count)
end

"""
    detect_collisions_sinr(finished_tx, sf, noise_dbm)

Detect collisions based on SINR processing with Capture Effect.
LoRa Capture Effect: Stronger signal can survive if SIR > Threshold (approx 6dB).
Updates `status` to "Collision" if packet fails to survive.
Returns (success_count, collision_count).
"""
function detect_collisions_sinr(finished_tx, sf, noise_dbm)
    n = length(finished_tx)
    if n == 0
        return (0, 0)
    end
    
    # Sort by start time
    sort!(finished_tx, by = x -> x.start_ms)
    
    # Keep track of failures
    failed_indices = Set{Int}()
    
    # Required SNR for LoRa (approximate values)
    # SF7:-7.5, SF8:-10, SF9:-12.5, SF10:-15, SF11:-17.5, SF12:-20 (dB)
    required_snr_db = -15.0 # Default for SF10
    if sf == 7 required_snr_db = -7.5
    elseif sf == 8 required_snr_db = -10.0
    elseif sf == 9 required_snr_db = -12.5
    elseif sf == 10 required_snr_db = -15.0
    elseif sf == 11 required_snr_db = -17.5
    elseif sf == 12 required_snr_db = -20.0
    end
    
    # Co-channel rejection (SIR threshold)
    # Usually around 6dB for same SF
    required_sir_db = 6.0
    
    # Noise Power in W
    noise_w = 10^(noise_dbm/10) * 1e-3
    
    for i in 1:n
        target = finished_tx[i]
        
        # Identify Interferers
        interferers = []
        
        # 1. Search backwards (potentially overlapping)
        k = i - 1
        while k >= 1
            other = finished_tx[k]
            # Since sorted by start_ms, other.start_ms <= target.start_ms.
            # We need overlap: max(start1, start2) < min(end1, end2)
            # here start2=target.start_ms. So target.start_ms < min(end1, target.end_ms)
            # Which simplifies to target.start_ms < end1 (other.end_ms).
            
            if other.end_ms <= target.start_ms
                # No overlap with this packet, but there might be earlier packets that end VERY late.
                # Unlike strict time-window, we don't have a guarantee that earlier start means earlier end.
                # So we should check.
                # However, for simulation efficiency, usually we can limit search range or just scan.
                # Given N is likely small (< 10000?), scanning a bit is fine.
            elseif other.channel == target.channel
                push!(interferers, other)
            end
            
            # Optimization: if start time difference is huge (larger than max packet duration), we can stop?
            # But we don't know max packet duration easily here.
            # Let's iterate all previous for correctness or use a reasonable window if needed.
            # For N=1000, O(N^2) is fine.
            k -= 1
        end
        
        # 2. Search forwards
        k = i + 1
        while k <= n
            other = finished_tx[k]
            if other.start_ms >= target.end_ms
                break # No intersection possible forwards
            end
            
            if other.channel == target.channel
                push!(interferers, other)
            end
            k += 1
        end
        
        # Calculate Interference Power
        interference_w = 0.0
        inf_ids = Int[]
        inf_dists = Float64[]
        
        for inf in interferers
            p_mw = 10^(inf.rx_power_at_gw / 10)
            interference_w += p_mw * 1e-3
            
            # Calculate distance between target and interferer
            dist = sqrt((target.x - inf.x)^2 + (target.y - inf.y)^2)
            push!(inf_ids, inf.terminal_id)
            push!(inf_dists, dist)
        end
        
        # Target Signal Power
        signal_w = 10^(target.rx_power_at_gw / 10) * 1e-3
        
        # Calculate SINR and SNR
        sinr_db = -Inf
        if !isempty(interferers) && interference_w > 0
            sir_w = signal_w / interference_w
            sinr_db = 10 * log10(sir_w)
        end
        
        snr_w = signal_w / noise_w
        snr_db = 10 * log10(snr_w)
        
        # Store metrics in TransmissionRecord
        target.num_interferers = length(interferers)
        target.sinr_db = sinr_db
        target.snr_db = snr_db
        target.interferer_ids = join(string.(inf_ids), ";")
        target.interferer_distances = join([@sprintf("%.1f", d) for d in inf_dists], ";")
        
        # Determine success/failure and reason
        is_success = true
        failure_reason = "None"
        
        # SIR Check (Capture Effect)
        if !isempty(interferers)
            if interference_w > 0
                if sinr_db < required_sir_db
                    is_success = false
                    failure_reason = "SINR_Insufficient"
                end
            else
                 # Should not happen if interferers list is not empty
            end
        end
        
        # SNR Check (Sensitivity)
        if snr_db < required_snr_db
            is_success = false
            # SNR failure takes precedence if both fail
            if failure_reason == "SINR_Insufficient"
                failure_reason = "SINR_and_SNR_Insufficient"
            else
                failure_reason = "SNR_Insufficient"
            end
        end
        
        # Update status and failure reason
        if !is_success
            push!(failed_indices, i)
            target.status = "Collision"
            target.failure_reason = failure_reason
        else
            target.failure_reason = "None"
        end
    end
    
    collision_count = length(failed_indices)
    success_count = n - collision_count
    
    return (success_count, collision_count)
end

