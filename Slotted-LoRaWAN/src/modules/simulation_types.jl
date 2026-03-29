module SimulationTypes

export SyncResult, CandidateSlot, TransmissionRecord, MACEvent, EventType, TX_START, TX_END, ACK_CHECK, PACKET_GEN

# ==========================================
# 1. 同期結果構造体
# ==========================================
struct SyncResult
    terminal_id::Int
    success::Bool
    sync_time_ms::Float64
    sync_error_ms::Float64
end

# ==========================================
# 2. スロット候補構造体
# ==========================================
struct CandidateSlot
    time_global_ms::Float64  # グローバル時刻での開始
    terminal_node            # 端末オブジェクト
    slot_index::Int
    channel::Int             # チャネル番号（1-based）
end

# ==========================================
# 3. 送信記録構造体
# ==========================================
mutable struct TransmissionRecord
    terminal_id::Int
    start_ms::Float64
    end_ms::Float64
    tx_power_dbm::Float64
    x::Float64
    y::Float64
    status::String
    rx_power_at_gw::Float64
    channel::Int  # チャネル番号（1-based）
    
    # === ACK/再送および分析用フィールド ===
    retry_count::Int              # 再送回数
    is_retransmission::Bool       # 再送フラグ
    packet_id::Int                # パケットID (original_packet_id)
    failure_reason::String        # 詳細な失敗理由
    sinr_db::Float64             # 実際のSINR値 (dB)
    snr_db::Float64              # 実際のSNR値 (dB)
    num_interferers::Int         # 干渉パケット数
    cs_detected::Bool            # キャリアセンスで干渉検出したか
    interferer_ids::String       # 干渉した端末IDリスト (例: "1;23;45")
    interferer_distances::String # 干渉した端末との距離 [m] リスト (例: "120.5;340.2")
end

# ==========================================
# 4. イベント管理用
# ==========================================
@enum EventType begin
    TX_START   # 送信開始（LBT含む）
    TX_END     # 送信終了
    ACK_CHECK  # ACK確認（再送判断）
    PACKET_GEN # パケット生起（バッファへの追加）
end

struct MACEvent
    time_ms::Float64
    type::EventType
    terminal_id::Int
    packet_id::Int  # original_packet_id
    data::Any       # 必要に応じて追加情報を格納 (CandidateSlot, TransmissionRecordなど)
end

end # module
