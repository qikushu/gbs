# Total rice chromosome length of Nipponbare assembly in Mb 
chrMaxLength <- c("1"=43.270923,	"2"=35.937250,	"3"=36.413819,	"4"=35.502694,"5"=29.958434,	"6"=31.248787,	"7"=29.697621,	"8"=28.443022,	"9"=23.012720,	"10"=23.207287,	"11"=29.02110612, "12"=27.531856)

# Conversion rule of ABH genotypes and color code in output graphics.
# color_mapping <- c("A" = "orange", "B" = "#00bfff", "H" = "red", "-" = "gray")

# Read data
plantGenoList = readGGdata("../Downloads/ABH_22WRCFn_8_ver2.0.csv")
# plot ans save in PNG file
png("output.png", width = 2400, height = 2000)
plotHorizontal(plantGenoList=plantGenoList,chrMaxLength=chrMaxLength, color_mapping=color_mapping)
dev.off()




readGGdata = function(infile) {

    # Read data
    rawdf = read.csv(infile, head=T)

    # 系統名を取得 (Plant列)
    plantNames = rawdf$Plant[-c(1, 2)]

    # 1行目を評価して、空でない列番号を取得する
    notEmptyIndex = which(rawdf[1,] != "")
    df1 = rawdf[,notEmptyIndex]
    genotypeData = as.data.frame(t(df1))
    colnames(genotypeData) = append(c("Chr", "Pos"),plantNames)
    genotypeData$Pos= as.numeric(genotypeData$Pos)  # numeric型に変換
    genotypeData$Marker = rownames(genotypeData) 

    # 植物ごとにデータを分割
    forDeletion = c("Marker", "Chr", "Pos")
    plantGenoList <- lapply(plantNames, function(col_name) {
        cols_to_include <- c(forDeletion, col_name)
        df_subset <- genotypeData[cols_to_include]
        colnames(df_subset) <- cols_to_include
        return(df_subset)
    })
    # 名前付きリストで各植物体体のデータを保持
    names(plantGenoList) <- plantNames

    return(plantGenoList)
}


plotHorizontal = function(plantGenoList,chrMaxLength, color_mapping) {

    ### 染色体を横並びに並べるコード
    # 染色体同士の間隔
    chrInterval = 5
    # 植物体同士の間隔を指定。
    plantInterval = 2
    # マーカーの横棒の長さを指定。
    width = 12

    # 左端のスペースの間隔を指定。
    leftSpace = 10
    # 右端のスペースの間隔を指定。
    rightSpace = 10
    # 上端のスペースの間隔を指定。
    topSpace = 0

    # 名前付きリストで各植物体体のデータを保持
    plantNames = names(plantGenoList)
    # 個体数
    plantNum = length(plantNames)
    # 染色体の名前
    chrNames = unique(plantGenoList[[1]]$Chr)
    # 染色体の本数
    numChr = length(chrNames)
    # 全染色体の名前
    allChrNames = names(chrMaxLength)
    # 染色体の本数
    allChrNum = length(allChrNames)

    # xの最大幅
    maxXwidth = leftSpace + sum(chrMaxLength) + numChr * chrInterval + rightSpace + 10

    # yの最大幅
    maxYwidth =  (width + plantInterval ) * plantNum

    # プロット領域の設定
    plot(1, type="n", xlim=c(-20, maxXwidth), ylim=c(0,maxYwidth), xlab="X軸", ylab="Y軸")


    # j番目の植物体を処理する j=1に初期化
    j=1
    # chr描画のためのy方向のoffset
    yoff = 0;


    # 名前付きリスト`chromosome_data`に対してループ処理
    for (plantName in plantNames) {
    cat(plantName)
    # 植物体ごとのdfを取得
    df <- plantGenoList[[plantName]]
    colnames(df)[colnames(df) == plantName] <- "Genotype"
    # chr描画のためのx方向のoffset
    xoff = 0;
    # y_start座標の計算
    y_start = topSpace + yoff
    # y_end座標の計算
    y_end = topSpace + yoff + width
    # y_mid
    y_mid = (y_start + y_end) / 2
    # 系統名を記載
    text(x = -20, y = y_start, labels = plantName, cex = 1.5)


    for (chr in allChrNames) {

        # 染色体を描画
        chrEnd = chrMaxLength[chr] # 染色体の端
        

        
        # ユーザーデータが存在するかを確認
        if (any(df$Chr == chr)) {
            # 染色体ごとに抽出
            df2 = df[df$Chr==chr,]

            # 隣接する遺伝子型が同じ場合の処理
            df3 = flankGeno(df2, chr)

            # 新しい列 color を追加して色情報を格納
            df3$color <- color_mapping[df3$Genotype]


            # データフレームdf3を一行ずつ処理する for ループ
            for (i in 1:nrow(df3)) {
                row <- df3[i, ]  # 一行を抽出
                # xleft座標の配列
                xleft <- df3$MinPos + xoff
                # xright座標の配列
                xright <- df3$MaxPos + xoff
                # genotype color
                color = df3$color
                rect(xleft= xleft, ytop = y_start, xright=xright, ybottom=y_end, col=color, border=0)
            }
            # データフレームを一行ずつ処理する for ループ
            #for (i in 1:nrow(df2)) {
            #    row <- df2[i, ]  # 一行を抽出
                # x座標の配列
            #    x <- row$Pos + xoff
                # genotype color
            #    color = row$color
            #    segments(x0 = x, y0 = y_start, x1 = x, y1 = y_end, col = color)
            #}
        }
        rect(xleft= 0 + xoff, ytop = y_start, xright=chrEnd + xoff, ybottom=y_end, col=NA, border="grey")
        xoff = xoff + chrEnd + chrInterval
    }
    yoff = yoff + width + plantInterval

    j = j + 1;
    }
}


# 隣接する遺伝子型が同じ場合をまとめたデータフレームを返す
flankGeno = function(df2, chr) {
    df2$SequenceNumber <- createSequenceDataFrame(df2$Genotype)$SequenceNumber 
    # SequenceNumberの最大値
    maxLinkGroup = max(df2$SequenceNumber)
    
    result_df <- data.frame()  # 空のデータフレームを初期化

    for (i in 1:maxLinkGroup) {
        df3 = df2[df2$SequenceNumber == i,]
        # 遺伝子型
        gt = df3$Genotype[1]
        maxpos = max(df3$Pos)
        minpos = min(df3$Pos)
        result_df <- rbind(result_df, data.frame(Chr = chr, MinPos = minpos, MaxPos = maxpos, Genotype = gt))
    }

    return(result_df)
}


# 隣接するマーカーの遺伝子型が同一の場合、同一の遠し番号を割り振る
createSequenceDataFrame <- function(alphabet_vector) {
  # 通し番号を格納するためのベクトルを初期化
  sequence_number <- rep(NA, length(alphabet_vector))

  # 最初の要素に通し番号 1 を設定
  sequence_number[1] <- 1

  # 通し番号を割り振る
  for (i in 2:length(alphabet_vector)) {
    if (alphabet_vector[i] == alphabet_vector[i - 1]) {
      sequence_number[i] <- sequence_number[i - 1]
    } else {
      sequence_number[i] <- sequence_number[i - 1] + 1
    }
  }

  # 結果をデータフレームとして返す
  result_df <- data.frame(Alphabet = alphabet_vector, SequenceNumber = sequence_number)
  return(result_df)
}
