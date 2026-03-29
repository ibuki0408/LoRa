module PacketGeneration

using Random, Statistics

export generate_poisson_arrival_times, generate_packet_sizes, analyze_traffic

"""
    generate_poisson_arrival_times(lambda_events_per_s, duration_s)
"""
function generate_poisson_arrival_times(lambda_events_per_s, duration_s)
    arrival_times = Float64[]
    current_time = 0.0
    while current_time < duration_s
        inter_arrival_time = -log(rand()) / lambda_events_per_s
        current_time += inter_arrival_time
        if current_time < duration_s
            push!(arrival_times, current_time * 1000.0) # Convert to ms
        end
    end
    return arrival_times
end

"""
    generate_packet_sizes(num_packets, mean_size_bytes)
"""
function generate_packet_sizes(num_packets, mean_size_bytes)
    return fill(Int(round(mean_size_bytes)), num_packets)
end

"""
    analyze_traffic(arrival_times_ms)
"""
function analyze_traffic(arrival_times_ms)
    if isempty(arrival_times_ms)
        println("トラフィックがありません。")
        return
    end
    n = length(arrival_times_ms)
    duration_s = maximum(arrival_times_ms) / 1000.0
    avg_rate = n / duration_s
    println("Traffic Analysis:")
    println("  Packets: $n")
    println("  Duration: $(round(duration_s, digits=2)) s")
    println("  Avg Rate: $(round(avg_rate, digits=2)) pkt/s")
end

end # module
