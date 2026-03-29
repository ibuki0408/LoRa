module LocalClock

using Random

export ClockDriftParameters, calculate_drifted_time

struct ClockDriftParameters
    ppm::Float64
    factor::Float64
end

function calculate_drifted_time(absolute_time_ms::Float64, drift_params::ClockDriftParameters)
    return absolute_time_ms * drift_params.factor
end

end # module