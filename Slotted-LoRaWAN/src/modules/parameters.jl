module Parameters

export IntegratedParameters, create_integrated_params

mutable struct IntegratedParameters
    # === 物理層 (PHY) ===
    signal_duration_us::Float64
    signal_bw_mhz::Float64
    terminal_bw_mhz::Float64
    tx_sampling_rate_mhz::Float64
    rx_sampling_rate_mhz::Float64
    tx_power_dbm::Float64
    noise_figure_db::Float64
    
    # === MAC層 ===
    num_terminals::Int
    area_size_m::Float64
    slot_length_ms::Float64
    packet_airtime_ms::Float64
    cs_threshold_dbm::Float64
    
    # === LoRa固有 ===
    spreading_factor::Int
    lora_payload_bytes::Int
    num_channels::Int
    
    # === シミュレーション制御 ===
    beacon_interval_ms::Float64
    simulation_duration_ms::Float64
    max_startup_delay_ms::Float64
    duty_cycle::Float64
    mean_event_interval_ms::Float64
    
    # === 環境モデル ===
    shadowing_enabled::Bool
    shadowing_std_db::Float64
    pass_loss_exp::Float64
    
    # === Out-of-band同期 ===
    sync_center_freq_ghz::Float64
    data_center_freq_ghz::Float64
    reference_path_loss_db::Float64
    
    # === 同期基地局位置 ===
    sync_bs_x_m::Float64
    sync_bs_y_m::Float64
    
    # === 同期検出 ===
    sync_observation_duration_ms::Float64
    gw_tx_power_dbm::Float64
    noise_floor_window_ms::Float64
    detection_margin_db::Float64
    min_samples::Int
    debounce_time_ms::Float64
    initial_window_duration_ms::Float64
    tx_jitter_max_ms::Float64
    
    # === 間欠受信 ===
    enable_intermittent_rx::Bool
    intermittent_window_ms::Float64
    initial_search_duration_ms::Float64
    
    # === 決定的ジッタ ===
    use_deterministic_jitter::Bool
    num_jitter_offsets::Int
    deterministic_jitter_random_ms::Float64
    
    collision_model::Symbol
    
    # === ACK/再送 (Optional in some scripts) ===
    enable_ack::Bool
    max_retries::Int
    ack_timeout_ms::Float64
    rx1_delay_ms::Float64
    backoff_base_ms::Float64
    
    # === LBT (Listen Before Talk) ===
    lbt_duration_ms::Float64
    lbt_sample_interval_ms::Float64
    
    # === 出力制御 ===
    enable_file_output::Bool
    enable_plot_output::Bool
    enable_detailed_logs::Bool
    force_async_mode::Bool
    max_buffer_size::Int
    
    # === PHY Toggles ===
    enable_carrier_sense::Bool
    enable_capture_effect::Bool
    target_sync_rate::Float64

    # === PHYモデル選択 ===
    use_waveform_phy::Bool  # true: 波形リサンプリングPHY / false: 電力モデルPHY（デフォルト）
end

function create_integrated_params()
    # デフォルト値の設定 (Prop.jl 基準)
    num_terminals = 100
    sf = 8
    beacon_interval_ms = 100.0
    slot_length_ms = 200.0
    mean_event_interval_ms = 30000.0
    simulation_duration_ms = 3600000.0
    
    tx_power_dbm = 13.0
    payload_bytes = 10
    num_channels = 8
    sync_freq_ghz = 3.7
    data_freq_ghz = 0.92
    
    area_size_m = 1000.0
    shadowing_std_db = 0.0
    pass_loss_exp = 2.7
    
    signal_duration_us = 66.67
    signal_bw_mhz = 1.905
    terminal_bw_mhz = 0.125
    tx_sampling_rate_mhz = 3.84
    rx_sampling_rate_mhz = 0.256
    noise_figure_db = 6.0
    
    sync_observation_duration_ms = 65000.0
    gw_tx_power_dbm = 40.0
    noise_floor_window_ms = 10.0
    detection_margin_db = 15.0
    min_samples = 3
    debounce_time_ms = 1.0
    initial_window_duration_ms = 400.0
    tx_jitter_max_ms = 10.0
    
    enable_intermittent_rx = true
    intermittent_window_ms = 6.0
    initial_search_duration_ms = 60.0
    use_deterministic_jitter = false
    collision_model = :sinr
    
    lbt_duration_ms = 5.0
    lbt_sample_interval_ms = 1.0
    cs_threshold_dbm = -80.0
    
    duty_cycle = 0.01
    max_startup_delay_ms = 60000.0
    force_async_mode = false
    max_buffer_size = 10
    
    enable_carrier_sense = true
    enable_capture_effect = true
    target_sync_rate = 1.0
    
    # ACK/再送デフォルト (Prop_Theory.jl 基準)
    enable_ack = false
    max_retries = 3
    ack_timeout_ms = 2000.0
    rx1_delay_ms = 1000.0
    backoff_base_ms = 1000.0

    return IntegratedParameters(
        signal_duration_us, signal_bw_mhz, terminal_bw_mhz,
        tx_sampling_rate_mhz, rx_sampling_rate_mhz, tx_power_dbm, noise_figure_db,
        
        num_terminals, area_size_m, slot_length_ms, 
        0.0, # packet_airtime_ms
        cs_threshold_dbm,
        
        sf, payload_bytes, num_channels,
        
        beacon_interval_ms, simulation_duration_ms, max_startup_delay_ms,
        duty_cycle, mean_event_interval_ms,
        
        true, # shadowing_enabled
        shadowing_std_db, pass_loss_exp,
        
        sync_freq_ghz, data_freq_ghz,
        20*log10(data_freq_ghz*1e9) - 147.55, # reference_path_loss_db
        
        0.0, 0.0, # sync_bs_x, sync_bs_y
        
        sync_observation_duration_ms, gw_tx_power_dbm, noise_floor_window_ms,
        detection_margin_db, min_samples, debounce_time_ms,
        initial_window_duration_ms, tx_jitter_max_ms,
        
        enable_intermittent_rx, intermittent_window_ms, initial_search_duration_ms,
        
        use_deterministic_jitter, 20, 5.0, # num_jitter_offsets, random_ms
        collision_model,
        
        enable_ack, max_retries, ack_timeout_ms, rx1_delay_ms, backoff_base_ms,
        
        lbt_duration_ms, lbt_sample_interval_ms,
        
        true, true, true, # enable_file, plot, detailed_logs
        force_async_mode, max_buffer_size,
        enable_carrier_sense, enable_capture_effect, target_sync_rate,
        true  # use_waveform_phy: デフォルトは電力モデル(false)
    )
end

end # module
