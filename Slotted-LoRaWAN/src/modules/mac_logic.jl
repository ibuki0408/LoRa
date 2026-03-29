module MACLogic

using ..Parameters
using ..SimulationTypes
using ..TerminalDeployment

export schedule_tx_start!

"""
    schedule_tx_start!(event_queue, me, packet_id, base_time, params, terminal_sync_infos, synced_terminal_ids)

指定された base_time 以降の最良のタイミング（同期済みならスロット境界）で送信イベントをスケジュールする。
"""
function schedule_tx_start!(event_queue, me, packet_id, base_time, params, terminal_sync_infos, synced_terminal_ids)
    # --- 送信スケジュール (ジッタとスロット同期のみ) ---
    effective_base_time = base_time

    if me.terminal_id in synced_terminal_ids
        # 同期成功端末: スロット境界にスナップ
        start_ms = terminal_sync_infos[me.terminal_id]
        slot_len = params.slot_length_ms * me.clock_drift_factor
        
        diff_from_base = effective_base_time - start_ms
        num_slots = ceil(max(0.0, diff_from_base) / slot_len)
        tx_time_snapped = start_ms + num_slots * slot_len
        
        # ジッタ
        guard_band_ms = 10.0
        max_jitter_in_slot = max(0.0, params.slot_length_ms - params.packet_airtime_ms - guard_band_ms)
        jitter_limit = min(params.tx_jitter_max_ms, max_jitter_in_slot)
        jitter = rand() * jitter_limit
        tx_time_final = tx_time_snapped + jitter
    else
        # 非同期端末: バックオフ後の時刻を使用
        tx_time_final = effective_base_time
    end
    
    # シミュレーション期間内であればイベント登録
    if tx_time_final < params.simulation_duration_ms
        channel = rand(1:params.num_channels)
        new_cand = CandidateSlot(tx_time_final, me, 1, channel)
        new_event = MACEvent(tx_time_final, TX_START, me.terminal_id, packet_id, new_cand)
        
        # 時刻降順を維持して挿入
        idx = searchsortedfirst(event_queue, new_event, by=x->x.time_ms, rev=true)
        insert!(event_queue, idx, new_event)
    end
end

end # module
