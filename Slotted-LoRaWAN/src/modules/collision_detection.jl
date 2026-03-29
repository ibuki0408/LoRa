module CollisionDetection

using Statistics, Printf
using ..PathLoss # For SINR detection distance calculation

export detect_collisions_simple_overlap, detect_collisions_sinr

"""
    detect_collisions_simple_overlap(finished_tx)
"""
function detect_collisions_simple_overlap(finished_tx)
    n = length(finished_tx)
    if n == 0
        return (0, 0)
    end

    sort!(finished_tx, by = x -> x.start_ms)
    
    collision_indices = Set{Int}()
    
    for i in 1:n
        tx1 = finished_tx[i]
        for j in (i+1):n
            tx2 = finished_tx[j]
            if tx2.start_ms >= tx1.end_ms
                break
            end
            if tx1.channel != tx2.channel
                continue
            end
            push!(collision_indices, i)
            push!(collision_indices, j)
        end
    end
    
    for idx in collision_indices
        finished_tx[idx].status = "Collision"
    end
    
    collision_count = length(collision_indices)
    success_count = n - collision_count
    
    return (success_count, collision_count)
end

"""
    detect_collisions_sinr(finished_tx, sf, noise_dbm)
"""
function detect_collisions_sinr(finished_tx, sf, noise_dbm)
    n = length(finished_tx)
    if n == 0
        return (0, 0)
    end
    
    sort!(finished_tx, by = x -> x.start_ms)
    failed_indices = Set{Int}()
    
    required_snr_db = -15.0
    if sf == 7 required_snr_db = -7.5
    elseif sf == 8 required_snr_db = -10.0
    elseif sf == 9 required_snr_db = -12.5
    elseif sf == 10 required_snr_db = -15.0
    elseif sf == 11 required_snr_db = -17.5
    elseif sf == 12 required_snr_db = -20.0
    end
    
    required_sir_db = 6.0
    noise_w = 10^(noise_dbm/10) * 1e-3
    
    for i in 1:n
        target = finished_tx[i]
        interferers = []
        
        k = i - 1
        while k >= 1
            other = finished_tx[k]
            if other.end_ms <= target.start_ms
                # Continue searching backwards as we don't have end-time sort guarantee
            elseif other.channel == target.channel
                push!(interferers, other)
            end
            k -= 1
        end
        
        k = i + 1
        while k <= n
            other = finished_tx[k]
            if other.start_ms >= target.end_ms
                break
            end
            if other.channel == target.channel
                push!(interferers, other)
            end
            k += 1
        end
        
        interference_w = 0.0
        inf_ids = Int[]
        inf_dists = Float64[]
        
        for inf in interferers
            p_mw = 10^(inf.rx_power_at_gw / 10)
            interference_w += p_mw * 1e-3
            dist = sqrt((target.x - inf.x)^2 + (target.y - inf.y)^2)
            push!(inf_ids, inf.terminal_id)
            push!(inf_dists, dist)
        end
        
        signal_w = 10^(target.rx_power_at_gw / 10) * 1e-3
        sinr_db = -Inf
        if !isempty(interferers) && interference_w > 0
            sinr_db = 10 * log10(signal_w / interference_w)
        end
        
        snr_db = 10 * log10(signal_w / noise_w)
        target.num_interferers = length(interferers)
        target.sinr_db = sinr_db
        target.snr_db = snr_db
        target.interferer_ids = join(string.(inf_ids), ";")
        target.interferer_distances = join([@sprintf("%.1f", d) for d in inf_dists], ";")
        
        is_success = true
        failure_reason = "None"
        
        if !isempty(interferers) && interference_w > 0 && sinr_db < required_sir_db
            is_success = false
            failure_reason = "SINR_Insufficient"
        end
        
        if snr_db < required_snr_db
            is_success = false
            failure_reason = (failure_reason == "SINR_Insufficient") ? "SINR_and_SNR_Insufficient" : "SNR_Insufficient"
        end
        
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

end # module
