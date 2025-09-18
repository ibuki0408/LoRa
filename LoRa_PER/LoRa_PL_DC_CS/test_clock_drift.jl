using Random, FFTW, Statistics, Printf, DelimitedFiles, LinearAlgebra, StatsBase

# クロックドリフトテスト用の簡易版
include("LoRa_SlottedALOHA_DC_CS_ClockDrift.jl")

println("=== クロックドリフトテスト開始 ===")

# テストパラメータ
sf = 7
bw = 125e3
num_devices = 2
sim_time = 1.0  # 短時間テスト
slot_duration = 0.1
dc_limit = 0.01
dc_window = 3600.0
cs_duration = 0.005
backoff_base = 0.1
backoff_max = 10.0
clock_accuracy = :standard
fixed_snr = -10.0  # 固定SNRでテスト
use_shadowing = false  # シャドウイング無効でテスト
rng = MersenneTwister(1234)

println("テストパラメータ:")
println("  端末数: $num_devices")
println("  シミュレーション時間: $sim_time 秒")
println("  スロット時間: $slot_duration 秒")
println("  固定SNR: $fixed_snr dB")
println("  クロック精度: $clock_accuracy")

# クロックドリフトテスト実行
per, positions, snrs_actual, stats, clock_drifts = slotted_aloha_dc_cs_clock_drift_per(
    sf, bw, num_devices;
    payload_range=(16,32),
    sim_time=sim_time,
    slot_duration=slot_duration,
    dc_limit=dc_limit,
    dc_window=dc_window,
    cs_duration=cs_duration,
    backoff_base=backoff_base,
    backoff_max=backoff_max,
    clock_accuracy=clock_accuracy,
    fixed_snr=fixed_snr,
    use_shadowing=use_shadowing,
    rng=rng
)

println("\n=== テスト結果 ===")
println("PER: $(per)")
println("総パケット数: $(stats[:total_packets])")
println("成功パケット数: $(stats[:success_packets])")
println("衝突パケット数: $(stats[:collision_packets])")
println("復調失敗パケット数: $(stats[:demod_fail_packets])")
println("キャリアセンスブロック数: $(stats[:cs_blocked_packets])")
println("バックオフブロック数: $(stats[:backoff_packets])")
println("クロックドリフトブロック数: $(stats[:clock_drift_packets])")

println("\n=== クロックドリフト情報 ===")
for (i, clock) in enumerate(clock_drifts)
    println("端末 $i:")
    println("  ドリフト: $(clock.drift_ppm) ppm")
    println("  蓄積ずれ: $(clock.accumulated_drift) 秒")
    println("  最後の更新時刻: $(clock.last_update_time) 秒")
end

println("\n=== テスト完了 ===")
