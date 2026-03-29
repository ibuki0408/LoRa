# ===== slot_simulation.jl =====
# 端末ごとの非同期送信ロジックモジュール
# ★★★ 2025/11/17 修正: Slot_id廃止、非同期起動対応 ★★★

using DataFrames, Dates, Random

"""
渡された「検出済み同期信号」に基づき、B案ロジックで送信イベントを生成する。
"""
function simulate_slot_transmission(
    valid_crossings::Dict{Symbol, Vector{Float64}}, # 既に窓内でフィルタリングされた信号
    sim_params::SimulationParameters, 
    terminal_clock::LocalClock,
    window_end_local_ms::Float64 # 受信窓が閉じるローカル時刻
)
    
    # パラメータ取得
    T_slot_ms = sim_params.slot_length_ms
    T_toa_ms = sim_params.lora_toa_ms
    T_sync_ms = sim_params.signal_interval_ms
    p_tx = sim_params.transmission_probability
    num_ch = sim_params.num_channels

    # 1. 基準時刻の取得
    # (Main側ですでにフィルタリングされているので、リストの先頭が「窓内の最初の信号」)
    all_local_peak_times = valid_crossings[:local_peak_times]
    
    if isempty(all_local_peak_times)
        return DataFrame()
    end
    
    t_first_local_ms = all_local_peak_times[1]
    
    # 2. スロット開始時刻の計算 (B案)
    # 窓の終わりまでに期待される最後のビーコン時刻
    num_beacons_after_first = floor((window_end_local_ms - t_first_local_ms) / T_sync_ms)
    t_slot_start_local_ms = t_first_local_ms + (num_beacons_after_first * T_sync_ms)

    # 3. 送信イベントDataFrame (slot_id を廃止 -> seq_no)
    tx_events = DataFrame(
        seq_no = Int[],                   # 端末ごとの通し番号
        target_tx_start_local_ms = Float64[], 
        actual_tx_start_global_ms = Float64[],
        actual_tx_end_global_ms = Float64[],
        channel_id = Int[]
    )

    # 4. 送信ループ
    current_local_target_ms = t_slot_start_local_ms
    seq_counter = 1

    # シミュレーション終了までループ
    while current_local_target_ms <= sim_params.total_duration_ms
        
        # トラフィック生成判定
        if rand() < p_tx
            selected_channel = rand(1:num_ch)
            
            target_tx_time_local_s = current_local_target_ms / 1000.0
            actual_tx_start_global_s = convert_to_global_time(terminal_clock, target_tx_time_local_s)
            
            # 送信終了時刻
            lora_toa_s = T_toa_ms / 1000.0
            actual_tx_end_global_s = actual_tx_start_global_s + lora_toa_s

            push!(tx_events, (
                seq_counter,
                current_local_target_ms,
                actual_tx_start_global_s * 1000.0,
                actual_tx_end_global_s * 1000.0,
                selected_channel
            ))
            seq_counter += 1
        end
        
        current_local_target_ms += T_slot_ms
    end
    
    return tx_events
end