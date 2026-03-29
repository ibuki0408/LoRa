module Utilities

using Statistics

export estimate_noise_floor_integrated, detect_crossings_integrated, debounce_integrated

"""
    estimate_noise_floor_integrated(rx_power, time_axis_ms; pre_signal_end_ms=9.0)

Estimate the noise floor from the early part of the received signal.
"""
function estimate_noise_floor_integrated(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}; pre_signal_end_ms::Float64=9.0)
    idxs = findall(t -> t < pre_signal_end_ms, time_axis_ms)
    if isempty(idxs)
        noise_power_w = median(rx_power)
        return 10 * log10(noise_power_w * 1000)
    end
    segment = rx_power[idxs]
    noise_power_w = median(segment)
    return 10 * log10(noise_power_w * 1000)
end

"""
    detect_crossings_integrated(rx_power, time_axis_ms, threshold_w)

Detect threshold crossings with linear interpolation for high precision.
"""
function detect_crossings_integrated(rx_power::Vector{Float64}, time_axis_ms::Vector{Float64}, threshold_w::Float64)
    crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[],
        :end_times => Float64[],
        :peak_times => Float64[],
        :peak_powers_w => Float64[],
        :peak_powers_dbm => Float64[]
    )
    
    is_above_threshold = false
    segment_start_idx = 0
    segment_peak_idx = 0
    segment_peak_power = 0.0
    
    for i in 2:length(rx_power)
        p_prev = rx_power[i-1]
        p_curr = rx_power[i]
        t_prev = time_axis_ms[i-1]
        t_curr = time_axis_ms[i]
        
        if !is_above_threshold && p_prev <= threshold_w && p_curr > threshold_w
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            
            push!(crossings[:start_times], crossing_time)
            is_above_threshold = true
            segment_start_idx = i
            segment_peak_idx = i
            segment_peak_power = p_curr
        end
        
        if is_above_threshold
            if p_curr > segment_peak_power
                segment_peak_idx = i
                segment_peak_power = p_curr
            end
        end
        
        if is_above_threshold && p_prev > threshold_w && p_curr <= threshold_w
            frac = (threshold_w - p_prev) / (p_curr - p_prev + 1e-20)
            crossing_time = t_prev + frac * (t_curr - t_prev)
            
            push!(crossings[:end_times], crossing_time)
            
            peak_time = time_axis_ms[segment_peak_idx]
            push!(crossings[:peak_times], peak_time)
            push!(crossings[:peak_powers_w], segment_peak_power)
            push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
            
            is_above_threshold = false
        end
    end
    
    if is_above_threshold
        if !isempty(time_axis_ms)
            push!(crossings[:end_times], time_axis_ms[end])
            peak_time = time_axis_ms[segment_peak_idx]
            push!(crossings[:peak_times], peak_time)
            push!(crossings[:peak_powers_w], segment_peak_power)
            push!(crossings[:peak_powers_dbm], 10 * log10(segment_peak_power * 1000))
        end
    end
    
    return crossings
end

"""
    debounce_integrated(crossings, debounce_time_ms)

Filter out spurious crossings within a debounce time.
"""
function debounce_integrated(crossings::Dict{Symbol, Vector{Float64}}, debounce_time_ms::Float64)
    if isempty(crossings[:peak_times])
        return crossings
    end

    filtered_crossings = Dict{Symbol, Vector{Float64}}(
        :start_times => Float64[],
        :end_times => Float64[],
        :peak_times => Float64[],
        :peak_powers_w => Float64[],
        :peak_powers_dbm => Float64[]
    )

    last_accepted_peak_time_ms = -Inf

    for i in 1:length(crossings[:peak_times])
        current_peak_time = crossings[:peak_times][i]

        if current_peak_time >= last_accepted_peak_time_ms + debounce_time_ms
            push!(filtered_crossings[:start_times], crossings[:start_times][i])
            push!(filtered_crossings[:end_times], crossings[:end_times][i])
            push!(filtered_crossings[:peak_times], current_peak_time)
            push!(filtered_crossings[:peak_powers_w], crossings[:peak_powers_w][i])
            push!(filtered_crossings[:peak_powers_dbm], crossings[:peak_powers_dbm][i])

            last_accepted_peak_time_ms = current_peak_time
        end
    end

    return filtered_crossings
end

end # module
