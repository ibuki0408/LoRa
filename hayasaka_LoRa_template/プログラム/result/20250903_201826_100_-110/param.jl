#1シミュレーション時間(日)
const SIM_PERIOD = 24 * 60.0 * 3
# const SIM_PERIOD = 24 * 60.0 * 1
# const SIM_PERIOD = 60.0
# const SIM_PERIOD_total = 24 * 60.0

#エリアサイズ(km)
const area_size = 0.5

#周波数幅(Hz)
const BW = 125000.0

#チャネル数
const channel_num = 4

#クラスタ数
const cluster_num = 4

#クラスタリング時間
# const clustering_time1 = 60.0 * 36
# const clustering_time2 = 60.0 * 72

#CDF計算時間
# const cdf_time1 = [10.0, 2800.0]
# const cdf_time2 = [3000.0, 5500.0]
# const cdf_time3 = [5501.0, 5505.0]
# const cdf_time1 = [5.0, 9.0]
# const cdf_time2 = [10.0, 2000.0]
# const cdf_time3 = [2200.0, 4000.0]
# const cdf_time1 = [10.0, 55.0]
# const cdf_time2 = [100.0, 1400.0]

#光速
const c = 3.0 * 10^8

#パスロス係数
const α = 4.0

#伝搬損失オフセット
const β = 9.5

#キャリアセンス閾値
const carrier_sense_threshold = -110
# const carrier_sense_threshold = -90
# const carrier_sense_threshold = -Inf
# const carrier_sense_threshold = Inf

#伝搬周波数係数
const γ = 4.5

#搬送波周波数
const f_c = 923.2

#シャドウイング標準偏差[dB]
const shadowing_standard_deviation = 3.48

#シャドウイング 素波の数
const SW_N = 500

#シャドウイング相関距離[m]
const dcor = 50.0

const noise_figure = 10 #雑音指数[dB]
const band_width = 125 * 10^3 #帯域幅
const noise_power_spectrum_density = -174 #雑音スペクトラム密度

const margin_db = 5


#送信電力(dBm)
const Tx_dB = 13

#SNRの閾値[SF, threshold]
const SNR_threshold_list = [[6, -5], [7, -7.5], [8, -10], [9, -12.5], [10, -15], [11, -17.5], [12, -20]]

#SIRの最悪値
# const SIR_threshold_worst = [[6, -11, -11, -11, -11, -11, -11], [-11, 6, -11, -11, -11, -11, -11], [-13, -13, 6, -13, -13, -13, -13], [-16, -16, -16, 6, -16, -16, -16], [-19, -19, -19, -19, 6, -19, -19], [-25, -25, -25, -25, -25, 6, -25], [-23, -23, -23, -23, -23, -23, 6]]

#SIRの閾値[干渉信号のSF(SF=6からSF=12)[参照信号のSF(SF=6からSF=12)]]
# const SIR_threshold_list = [[6, -11, -14, -17, -19, -22, -24], [-8, 6, -13, -17, -19, -22, -24], [-10, -11, 6, -16, -19, -22, -24], [-11, -13, -14, 6, -19, -22, -25], [-11, -14, -16, -17, 6, -22, -25], [-11, -14, -17, -19, -20, 6, -25], [-11, -14, -17, -20, -22, -23, 6]]
const SIR_threshold_list = [0 -8 -10 -11 -11 -11 -11; -11 0 -11 -13 -14 -14 -14; -14 -13 0 -14 -16 -17 -17; -17 -17 -16 0 -17 -19 -20; -19 -19 -19 -19 0 -20 -22; -22 -22 -22 -22 -22 0 -23; -24 -24 -24 -25 -25 -25 0]


#ENからのパケットに対するACK
ACK_EN_ONOFF = 0

#各パケット時間長(分)

#パケット送信周期(分)
const packet_period = 2.0

#サブフレーム数
const subframe_number = 4

#スロット
const slot_number = 150

#リソース数
# const number_PLIM_bits = Int(floor(log2(channel_num * slot_number)))

#リソース数
const resource_number = channel_num * slot_number

#ペイロードサイズ
const payload_bit = 160

#符号化率
const CR = 4 / 7

#クロックドリフト平均最小
const CD_mean_min = -1.91 * 10^(-3)

#クロックドリフト平均最大
const CD_mean_max = 0.28 * 10^(-3)

#クロックドリフト分散最小
const CD_variance_min = 9.59 * 10^(-11)

#クロックドリフト分散最大
const CD_variance_max = 3.19 * 10^(-10)
