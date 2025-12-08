module PacketGeneration

export generate_next_poisson_time

"""
    generate_next_poisson_time(current_time_ms, mean_interval_ms, drift_factor, next_available_time_ms)

Generate the next packet transmission time based on a Poisson process (exponentially distributed intervals),
while strictly enforcing a Duty Cycle constraint.

# Arguments
- `current_time_ms`: The current simulation time (or the time of the last event).
- `mean_interval_ms`: The mean interval for the Poisson process (1/lambda).
- `drift_factor`: The clock drift factor for the terminal (e.g., 1.000005).
- `next_available_time_ms`: The earliest time the terminal is allowed to transmit due to Duty Cycle.

# Returns
- `actual_next_tx`: The calculated next transmission time (Float64).
"""
function generate_next_poisson_time(current_time_ms::Float64, mean_interval_ms::Float64, drift_factor::Float64, next_available_time_ms::Float64)
    # 1. Generate Poisson interval (Exponential Distribution)
    # interval = -mean * ln(rand)
    # Apply clock drift to the interval measurement
    interval = -mean_interval_ms * log(rand()) * drift_factor
    
    # 2. Calculate target time based on the random interval
    target_time = current_time_ms + interval
    
    # 3. Enforce Duty Cycle Constraint (Clamping)
    # If the target time is earlier than the allowed DC time, we must wait.
    # This effectively creates a "Poisson process with dead time".
    actual_next_tx = max(target_time, next_available_time_ms)
    
    return actual_next_tx
end

end
