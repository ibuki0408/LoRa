include("src/Prop.jl")
using DSP, Statistics
using .Parameters
using .SignalGeneration

params = create_integrated_params()
sig_gen_params = SignalParameters(
    params.signal_duration_us * 1e-6,
    params.sync_center_freq_ghz * 1e9,
    params.signal_bw_mhz * 1e6,
    params.tx_sampling_rate_mhz * 1e6,
    params.tx_power_dbm
)

gen_start_ms = 0.0
gen_end_ms = 10.0
beacon_times = [1.0]

(sig_tx, _) = generate_sync_signal_segment(sig_gen_params, beacon_times, gen_start_ms, gen_end_ms)
rate_ratio = params.rx_sampling_rate_mhz / params.tx_sampling_rate_mhz
downsample_ratio = Int(round(params.tx_sampling_rate_mhz / params.rx_sampling_rate_mhz))

sig_rx_decimated = sig_tx[1:downsample_ratio:end]
sig_rx_resampled = resample(sig_tx, rate_ratio)

println("Original (TX) Peak Power: ", maximum(abs2.(sig_tx)))
println("Original (TX) Mean Power: ", mean(abs2.(sig_tx)))
println("Decimated (false) Peak Power: ", maximum(abs2.(sig_rx_decimated)))
println("Decimated (false) Mean Power: ", mean(abs2.(sig_rx_decimated)))
println("Resampled (true) Peak Power: ", maximum(abs2.(sig_rx_resampled)))
println("Resampled (true) Mean Power: ", mean(abs2.(sig_rx_resampled)))

