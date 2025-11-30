# ===== collision_analyzer.jl =====
# 生成された送信イベントCSVを読み込み、衝突率を計算するツール

using CSV, DataFrames, Dates, Printf

function analyze_collisions(filepath::String)
    println("\n" * "="^60)
    println("衝突解析レポート")
    println("対象ファイル: $filepath")
    println("="^60)
    
    # CSV読み込み
    if !isfile(filepath)
        println("エラー: ファイルが見つかりません")
        return
    end
    
    df = CSV.read(filepath, DataFrame)
    total_packets = nrow(df)
    
    if total_packets == 0
        println("データがありません。")
        return
    end

    # 判定用フラグ配列 (trueなら衝突失敗)
    is_collided = falses(total_packets)
    
    # === 衝突判定ロジック (総当たり) ===
    # パケットAとパケットBについて：
    # 1. チャネルが同じ (channel_id)
    # 2. 時間が重なっている (Start_A < End_B  かつ  Start_B < End_A)
    # 両方満たせば衝突とする。

    # ※計算量削減のため、Start時間でソートしておく
    sort!(df, :actual_tx_start_global_ms)

    for i in 1:total_packets
        # 自分より後ろのパケットと比較
        for j in (i+1):total_packets
            
            # 時間最適化: 
            # パケットJの開始時刻が、パケットIの終了時刻を超えていたら、
            # それ以降のパケットKも全て重ならないのでループを抜ける
            if df[j, :actual_tx_start_global_ms] >= df[i, :actual_tx_end_global_ms]
                break
            end

            # 1. チャネルチェック
            if df[i, :channel_id] != df[j, :channel_id]
                continue # チャネルが違うのでセーフ
            end
            
            # 2. 時間重複チェック (ソート済み＆break条件があるので、ここは必ず重なっている)
            # 念のため厳密な条件: max(start_i, start_j) < min(end_i, end_j)
            
            start_i = df[i, :actual_tx_start_global_ms]
            end_i   = df[i, :actual_tx_end_global_ms]
            start_j = df[j, :actual_tx_start_global_ms]
            end_j   = df[j, :actual_tx_end_global_ms]
            
            if max(start_i, start_j) < min(end_i, end_j)
                # 衝突発生！ 両方を失敗にする
                is_collided[i] = true
                is_collided[j] = true
            end
        end
    end
    
    # === 結果集計 ===
    collision_count = count(is_collided)
    success_count = total_packets - collision_count
    per = collision_count / total_packets # Packet Error Rate
    throughput = success_count / total_packets

    println("解析結果:")
    println("• 総パケット数:   $total_packets")
    println("• 成功パケット数: $success_count")
    println("• 衝突パケット数: $collision_count")
    println("-"^30)
    println("• パケット衝突率 (PER): $(round(per * 100, digits=2)) %")
    println("• スループット (成功率): $(round(throughput * 100, digits=2)) %")
    println("="^60)
end

# ===== 自動実行ロジック =====
# results_multi_terminal フォルダ内の最新の all_tx_events_*.csv を探して実行する
function run_latest_analysis()
    target_dir = "results_multi_terminal"
    if !isdir(target_dir)
        println("ディレクトリ '$target_dir' が見つかりません。先にシミュレーションを実行してください。")
        return
    end

    files = readdir(target_dir, join=true)
    # ファイル名が "all_tx_events_" で始まり ".csv" で終わるものを抽出
    csv_files = filter(x -> occursin("all_tx_events_", x) && endswith(x, ".csv"), files)

    if isempty(csv_files)
        println("解析対象のCSVファイルが見つかりません。")
        return
    end

    # 最終更新日時またはファイル名でソートして最新を取得
    # (ファイル名にタイムスタンプが入っているので名前ソートでOK)
    latest_file = sort(csv_files)[end]
    
    analyze_collisions(latest_file)
end

# 実行
if abspath(PROGRAM_FILE) == @__FILE__
    run_latest_analysis()
end