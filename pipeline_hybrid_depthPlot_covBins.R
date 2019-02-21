options(scipen=999)
require(modes)
library("ggplot2")
library("dplyr")
args <- commandArgs(TRUE)
outputPrefix <- args[1]
window <- args[2]
bedData <- read.table(paste(outputPrefix, window, "avgDepth-d.txt", sep="_"), header=T)
completeChromList <- unlist(list(unique(bedData$chrom)))
output <- data.frame()
covMeanFile <- paste(outputPrefix, window, "avgDepthSummary.txt", sep="_")
stepSize <- bedData[2,1]
covThresholds <- read.table(paste(outputPrefix, ".cov", sep=""), header=T)
lowCov <- covThresholds[1,2]
if (lowCov<3){
  lowCov <- 3
} else if (lowCov>10){
  lowCov <- 10
}
introgressCutoff <- 0.01

median <- median(bedData$meanValue)
write(paste("median avgDepth =", median, sep=" "), file=covMeanFile, append=T)
mean <- mean(bedData$meanValue)
write(paste("Mean avgDepth =", mean, sep=" "), file=covMeanFile, append=T)
maxMedian <- median(bedData$meanValue)
maxMean <- max(bedData$meanValue)
meanQuant <- quantile(bedData$meanValue, 0.99, na.rm=T)
bedData$log2 <- log2(bedData$meanValue/mean)
bedData$log2[bedData$log2<0] <- 0
logQuant <- quantile(bedData$log2, 0.99, na.rm=T)
bedData$meanSkewLim <- bedData$meanValue
if (maxMean>(meanQuant*4)){
  mean997 <- quantile(bedData$meanValue, 0.997, na.rm=T)
  bedData$meanSkewLim[bedData$meanSkewLim>mean997] <- mean997
}

split <- strsplit(as.character(bedData$chrom), "_")
bedData$species <- unlist(lapply(split, function(x) x[1]))
bedData$chrNum <- unlist(lapply(split, function(x) x[2]))

bedData_filtered <- filter(bedData, !grepl("chrmt|chrMT|micron", chrom))
median <- median(bedData$meanValue)
mean <- mean(bedData$meanValue)
medianFiltered <- median(bedData_filtered$meanValue)
meanFiltered <- mean(bedData_filtered$meanValue)

species <- c()
chrs <- c()
#Split chr name from spcecies name
for (i in 1:length(completeChromList)) {
  completeChrString <- toString(completeChromList[[i]])
  split <- unlist(strsplit(completeChrString, "_"))
  species <- c(species, split[1])
  chrs <- c(chrs, split[2])
}
uniSpecies <- unique(species)
uniChrs <- unique(chrs)
justChrList <- uniChrs[!grepl("chrmt|chrMT|chr2", uniChrs)]
if (length(justChrList)==16){
  if (justChrList[5]!="chrV") {
    justChrList <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")
  }
}

data <- bedData$log2
ampOut <- amps(data)
#plot(ampOut$Peaks[,1]~ampOut$Peaks[,2])
antimodes <- ampOut$Antimode[,1]
log2binInfo <- data.frame("binName"=character(), "log2BinStart"=numeric(), "log2BinEnd"=numeric(), "numWin"=numeric(), "perData"=numeric(),stringsAsFactors=FALSE)
k=1
for (i in 1:(length(antimodes)+1)) {
  if (i==1){
    binStart<-0
  }
  if (i==length(antimodes)+1) {
    binEnd <- max(data)
  } else {
    binEnd <- antimodes[i]
  }
  numWin <- length(which(data>=binStart & data<=binEnd))
  binName <- paste("bin",(k-1),sep="")
  log2binInfo[k,1] <- binName
  log2binInfo[k,2] <- binStart
  log2binInfo[k,3] <- binEnd
  log2binInfo[k,4] <- numWin
  log2binInfo[k,5] <- round(100*(numWin/length(data)),2)
  k <- k+1
  binStart <- binEnd
}
#head(log2binInfo)
numBinsWanted <- max(which(log2binInfo$numWin>=50))+1
log2binInfo$meanCovUpper <- (2^log2binInfo$log2BinEnd)*mean
log2BinsWanted <- log2binInfo[1:numBinsWanted,]
log2BinsToWrite <- log2binInfo[1:(numBinsWanted+1),]
log2BinsToWrite$binName[1] <- "belowThreshold"
log2BinsToWrite$binName[(numBinsWanted+1)] <- "aboveUpperBin"
log2BinsToWrite$log2BinEnd[(numBinsWanted+1)] <- log2binInfo$log2BinEnd[nrow(log2binInfo)]
log2BinsToWrite$numWin[(numBinsWanted+1)] <- sum(log2binInfo$numWin[(numBinsWanted+1):nrow(log2binInfo)])
log2BinsToWrite$perData[(numBinsWanted+1)] <- sum(log2binInfo$perData[(numBinsWanted+1):nrow(log2binInfo)])
log2BinsToWrite$meanCovUpper[(numBinsWanted+1)] <- log2binInfo$meanCovUpper[nrow(log2binInfo)]
#write.table(log2BinsToWrite, file=paste(outputPrefix, "_globalLog2BinsIncluded.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

data <- bedData$meanSkewLim
ampOut <- amps(data)
#plot(ampOut$Peaks[,1]~ampOut$Peaks[,2])
antimodes <- ampOut$Antimode[,1]
meanbinInfo <- data.frame("binName"=character(), "meanBinStart"=numeric(), "meanBinEnd"=numeric(), "numWin"=numeric(), "perData"=numeric(),stringsAsFactors=FALSE)
k=1
for (i in 1:(length(antimodes)+1)) {
  if (i==1){
    binStart<-0
  }
  if (i==length(antimodes)+1) {
    binEnd <- max(data)
  } else {
    binEnd <- antimodes[i]
  }
  numWin <- length(which(data>=binStart & data<=binEnd))
  binName <- paste("bin",(k-1),sep="")
  meanbinInfo[k,1] <- binName
  meanbinInfo[k,2] <- round(binStart, digits=3)
  meanbinInfo[k,3] <- round(binEnd, digits=3)
  meanbinInfo[k,4] <- numWin
  meanbinInfo[k,5] <- round(100*(numWin/length(data)),2)
  k <- k+1
  binStart <- binEnd
}
#head(meanbinInfo)
numBinsWanted <- max(which(meanbinInfo$numWin>=50))+1
#meanbinInfo$meanCovUpper <- (2^meanbinInfo$meanBinEnd)*mean
meanBinsWanted <- meanbinInfo[1:numBinsWanted,]
meanBinsToWrite <- meanbinInfo[1:(numBinsWanted+1),]
meanBinsToWrite$binName[1] <- "belowThreshold"
meanBinsToWrite$binName[(numBinsWanted+1)] <- "aboveUpperBin"
meanBinsToWrite$meanBinEnd[(numBinsWanted+1)] <- meanbinInfo$meanBinEnd[nrow(meanbinInfo)]
meanBinsToWrite$numWin[(numBinsWanted+1)] <- sum(meanbinInfo$numWin[(numBinsWanted+1):nrow(meanbinInfo)])
meanBinsToWrite$perData[(numBinsWanted+1)] <- sum(meanbinInfo$perData[(numBinsWanted+1):nrow(meanbinInfo)])
write.table(meanBinsToWrite, file=paste(outputPrefix, "_globalMeanBinsIncluded.txt", sep=""), row.names = F, sep="\t", quote = FALSE)

log2BinOut <- bedData
log2BinOut$covBin <- "aboveUpperBin"
for (i in 1:nrow(log2BinsWanted)){
  binName <- log2BinsWanted$binName[i]
  log2BinOut$covBin[which(log2BinOut$log2>=log2BinsWanted[i,2] & log2BinOut$log2<=log2BinsWanted[i,3])] <- binName
}
#write.table(log2BinOut, file=paste(outputPrefix, "_winAvgDepth_wLog2CovBins.txt", sep=""), row.names = F, quote = FALSE)

meanBinOut <- bedData
meanBinOut$covBin <- "aboveUpperBin"
for (i in 1:nrow(meanBinsWanted)){
  binName <- meanBinsWanted$binName[i]
  meanBinOut$covBin[which(meanBinOut$meanSkewLim>=meanBinsWanted[i,2] & meanBinOut$meanSkewLim<=meanBinsWanted[i,3])] <- binName
}
write.table(meanBinOut, file=paste(outputPrefix, "_winAvgDepth_wMeanCovBins.txt", sep=""), row.names = F, quote = FALSE)

log2LimitValue <- log2binInfo$meanCovUpper[nrow(log2BinsWanted)+1]
bedData$meanValueLog2Limited <- bedData$meanSkewLim
bedData$meanValueLog2Limited[bedData$meanValueLimited>log2LimitValue] <- log2LimitValue
log2LowerCutoff <- log2BinsWanted$log2BinEnd[1]
log2UpperCutoff <- log2BinsWanted$log2BinEnd[nrow(log2BinsWanted)]

meanLimitValue <- meanbinInfo$meanCovUpper[nrow(meanBinsWanted)+1]
bedData$meanValuemeanLimited <- bedData$meanSkewLim
bedData$meanValuemeanLimited[bedData$meanValueLimited>meanLimitValue] <- meanLimitValue
meanLowerCutoff <- meanBinsWanted$meanBinEnd[1]
meanUpperCutoff <- meanBinsWanted$meanBinEnd[nrow(meanBinsWanted)]

log2SpSum <- data.frame(species=character(), perAboveCutoff=character(), contributes=character(), stringsAsFactors = F)
binSum <- c()
k=1
for (spName in uniSpecies){
  spData <- bedData[which(bedData$species==spName),]
  nWinPer <- nrow(spData)*introgressCutoff
  cutoffData <- spData[which(spData$log2>=log2LowerCutoff),]
  percentage <- paste(round((nrow(cutoffData)/nrow(spData))*100, 2), "%" , sep="")
  log2SpSum[k,1] <- spName
  log2SpSum[k,2] <- percentage
  log2SpSum[k,3] <- nrow(cutoffData)>=nWinPer
  binString <- spName
  for (j in 1:nrow(log2BinsWanted)+1){
    if (j == nrow(log2BinsWanted)+1) {
      binData <- spData[which(spData$log2>=log2BinsWanted$log2BinEnd[j-1]),]
      binPer <- paste(round((nrow(binData)/nrow(spData))*100, 2), "%" , sep="")
      binString <- paste(binString, binPer, sep=",")
    } else {
      binData <- spData[which(spData$log2>=log2BinsWanted$log2BinStart[j] & spData$log2<=log2BinsWanted$log2BinEnd[j]),]
      binPer <- paste(round((nrow(binData)/nrow(spData))*100, 2), "%" , sep="")
      binString <- paste(binString, binPer, sep=",")
    }
  }
  if (k==1) {
    binSum <- binString
  } else {
    binSum <- c(binSum, binString)
  }
  k <- k+1
}
#print(binSum)
binHeader <- "Speices"
for (i in 1:nrow(log2BinsWanted)+1){
  if (i==nrow(log2BinsWanted)+1) {
    binName <- "avboveUpperBin"
    binHeader <- c(binHeader, binName)
  } else {
    binName <- paste("bin", i-1, ":", round(log2BinsWanted$meanCovUpper[i-1], 2), "X-", round(log2BinsWanted$meanCovUpper[i], 2), "X", sep="")
    binHeader <- c(binHeader, binName)
  }
}
binTemp <- data.frame(data=binSum)
binDF <- data.frame(do.call('rbind', strsplit(as.character(binTemp$data),',',fixed=TRUE)))
colnames(binDF) <- binHeader
log2SpSum <- cbind(log2SpSum, binDF[,2:ncol(binDF)])
#write.table(log2SpSum, file=paste(outputPrefix, "_log2CovBinsSummary.txt", sep=""), row.names = F, quote = FALSE)

meanSpSum <- data.frame(species=character(), perAboveCutoff=character(), contributes=character(), stringsAsFactors = F)
binSum <- c()
k=1
for (spName in uniSpecies){
  spData <- bedData[which(bedData$species==spName),]
  nWinPer <- nrow(spData)*introgressCutoff
  cutoffData <- spData[which(spData$meanSkewLim>=meanLowerCutoff),]
  percentage <- paste(round((nrow(cutoffData)/nrow(spData))*100, 2), "%" , sep="")
  meanSpSum[k,1] <- spName
  meanSpSum[k,2] <- percentage
  meanSpSum[k,3] <- nrow(cutoffData)>=nWinPer
  binString <- spName
  for (j in 1:nrow(meanBinsWanted)+1){
    if (j == nrow(meanBinsWanted)+1) {
      binData <- spData[which(spData$meanSkewLim>=meanBinsWanted$meanBinEnd[j-1]),]
      binPer <- paste(round((nrow(binData)/nrow(spData))*100, 2), "%" , sep="")
      binString <- paste(binString, binPer, sep=",")
    } else {
      binData <- spData[which(spData$meanSkewLim>=meanBinsWanted$meanBinStart[j] & spData$meanSkewLim<=meanBinsWanted$meanBinEnd[j]),]
      binPer <- paste(round((nrow(binData)/nrow(spData))*100, 2), "%" , sep="")
      binString <- paste(binString, binPer, sep=",")
    }
  }
  if (k==1) {
    binSum <- binString
  } else {
    binSum <- c(binSum, binString)
  }
  k <- k+1
}
#print(binSum)
binHeader <- "Speices"
for (i in 1:nrow(meanBinsWanted)+1){
  if (i==nrow(meanBinsWanted)+1) {
    binName <- "avboveUpperBin"
    binHeader <- c(binHeader, binName)
  } else {
    binName <- paste("bin", i-1, ":", round(meanBinsWanted$meanBinStart[i], 2), "X-", round(meanBinsWanted$meanBinEnd[i], 2), "X", sep="")
    binHeader <- c(binHeader, binName)
  }
}
binTemp <- data.frame(data=binSum)
binDF <- data.frame(do.call('rbind', strsplit(as.character(binTemp$data),',',fixed=TRUE)))
colnames(binDF) <- binHeader
meanSpSum <- cbind(meanSpSum, binDF[,2:ncol(binDF)])
write.table(meanSpSum, file=paste(outputPrefix, "_meanCovBinsSummary.txt", sep=""), row.names = F, quote = FALSE)

# sigSpecies <- spSum$species[which(spSum$contributes==T)]
# sigSpecies <- factor(sigSpecies, levels=sigSpecies)
# sigSpBins <- data.frame(species=character(), binName=character(), log2binStart=numeric(), log2binEnd=numeric())
# sigData <- data.frame()
# sigDataOut <- data.frame()
# sigSpColNames <- colnames(binOut)
# binOutTemp <- binOut
# for (spc in sigSpecies) {
#   sigData <- rbind(sigData, bedData[which(bedData$species==spc),])
#   spcWin <- binOutTemp[which(binOutTemp$species==spc),]
#   sigDataOut <- rbind(sigDataOut, spcWin)
#   spLog2 <- spcWin$log2
#   ampOut <- amps(spLog2)
#   antimodes <- ampOut$Antimode[,1]
#   binInfo <- data.frame("species"=character(), "binName"=character(), "log2BinStart"=numeric(), "log2BinEnd"=numeric(), "numWin"=numeric(), "perData"=numeric(),stringsAsFactors=FALSE)
#   k=1
#   for (i in 1:(length(antimodes)+1)) {
#     if (i==1){
#       binStart<-0
#     }
#     if (i==length(antimodes)+1) {
#       binEnd <- max(spLog2)
#     } else {
#       binEnd <- antimodes[i]
#     }
#     numWin <- length(which(spLog2>=binStart & spLog2<=binEnd))
#     binName <- paste("bin",(k-1),sep="")
#     binInfo[k,1] <- spc
#     binInfo[k,2] <- binName
#     binInfo[k,3] <- binStart
#     binInfo[k,4] <- binEnd
#     binInfo[k,5] <- numWin
#     binInfo[k,6] <- round(100*(numWin/length(spLog2)),2)
#     k <- k+1
#     binStart <- binEnd
#   }
#   #head(binInfo)
#   numBinsWanted <- max(which(binInfo$numWin>=500))+1
#   if (numBinsWanted=="-Inf") {
#     numBinsWanted <- 4
#   }
#   spBinsWanted <- binInfo[1:numBinsWanted,]
#   sigSpBins <- rbind(sigSpBins, spBinsWanted[,1:4])
#   newColNum <- ncol(sigDataOut)+1
#   sigDataOut[newColNum] <- "NA"
#   binOutTemp[newColNum] <- "NA"
#   for (i in 1:(nrow(spBinsWanted)+1)){
#     if (i==(nrow(spBinsWanted)+1)){
#       binName <- "aboveUpperBin"
#       sigDataOut[which(sigDataOut$log2>=spBinsWanted[(i-1),4] & sigDataOut$species==spc),newColNum] <- binName
#     } else {
#       binName <- spBinsWanted$binName[i]
#       sigDataOut[which(sigDataOut$log2>=spBinsWanted[i,3] & sigDataOut$log2<=spBinsWanted[i,4] & sigDataOut$species==spc),newColNum] <- binName
#     }
#   }
#   sigSpColNames <- c(sigSpColNames, paste(spc, "covBin", sep="-"))
# }
# colnames(sigDataOut) <- sigSpColNames
# sigSpBins$meanCovUpper <- (2^sigSpBins$log2BinEnd)*mean
## Can add back once I get the bins working for each spcies
#write.table(sigDataOut, file=paste(outputPrefix, "_contributingSpecies_covBins.txt", sep=""), row.names = F)
#write.table(sigSpBins, file=paste(outputPrefix, "_contributingSpecies_covBinsInfo.txt", sep=""), row.names = F)
speciesColors <- c("Scer"="red", "Spar"="#FFDB00", "Smik"="#49FF00", "Skud"="#00FF92", "SkudZP"="#00FF92", "SkudIFO"="#00FF92", "Sarb"="#0092FF", "Suva"="#4900FF", "Seub"="#FF00DB")
fillScale <- scale_fill_manual(values=speciesColors, name="Species")
colorScale <- scale_color_manual(name="Species", values = colors, breaks=uniSpecies)
theme <- theme_bw()+theme(legend.text=element_text(face="italic"))

binAvgLines <- geom_vline(xintercept = meanBinsWanted$meanBinEnd, color="gray60", linetype="dotdash")
binLogLines <- geom_vline(xintercept =log2BinsWanted$log2BinEnd, color="gray30", linetype="dotdash")
#sigBinLogLines <- geom_vline(data=sigSpBins, aes(xintercept = log2BinEnd, color=species), linetype="dotdash")
hBinAvgLines <- geom_hline(yintercept = meanBinsWanted$meanBinEnd, color="gray60", linetype="dotdash")
hBinLogLines <- geom_hline(yintercept =log2BinsWanted$log2BinEnd, color="gray30", linetype="dotdash")
#hSigBinLogLines <- geom_hline(data=sigSpBins, aes(yintercept =log2BinEnd, color=species), linetype="dotdash")

pdf(paste(outputPrefix, "_covDistPlots.pdf", sep=""))

ggplot(bedData, aes(x=meanSkewLim))+geom_density(aes(y=..scaled..))+ggtitle(paste(outputPrefix, "Windowed Mean Coverage Density Distribution", sep=" ")) + coord_flip() + theme + binAvgLines
ggplot(bedData, aes(x=log2))+geom_density(aes(y=..scaled..))+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Density Distribution", sep=" ")) + coord_flip() + theme + binLogLines

ggplot(bedData, aes(x=meanSkewLim, fill=species))+geom_density(aes(y=..scaled..), alpha=.75)+fillScale+ggtitle(paste(outputPrefix, "Windowed Mean Coverage Species Density Distribution", sep=" ")) + coord_flip()+theme + binAvgLines
ggplot(bedData, aes(x=log2, fill=species))+geom_density(alpha=.75)+fillScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Species Density Distribution", sep=" ")) + coord_flip()+theme + binLogLines
ggplot(bedData, aes(x=species, y=meanSkewLim, fill=species))+geom_violin(scale="width")+fillScale+ggtitle(paste(outputPrefix, "Windowed Mean Coverage Violin Plot", sep=" "))+scale_x_discrete(limits = uniSpecies)+theme + hBinAvgLines
ggplot(bedData, aes(x=species, y=log2, fill=species))+geom_violin(scale="width")+fillScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Violin Plot", sep=" "))+scale_x_discrete(limits = uniSpecies)+theme + hBinLogLines

dev.off()

#ggplot(sigData, aes(x=log2, fill=species))+geom_density(alpha=.75)+fillScale+colorScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Species Density Distribution", sep=" ")) + coord_flip()+theme + binLogLines
#ggplot(sigData, aes(x=species, y=log2, fill=species))+geom_violin()+fillScale+colorScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Violin Plot", sep=" "))+scale_x_discrete(limits = sigSpecies)+theme  + hBinLogLines


plot <- ggplot()
totalFill <- c()
totalFillMedian <- c()
totalFillMean <- c()
totalFillMaxMedian <- c()
totalFillMaxMean <- c()
chrTotalFill <- c()
chrTotalFillMean <- c()
totalPoint <- c()
totalPointMedian <- c()
totalPointMean <- c()
chrTotalPoint <- c()
chrTotalPointMean <- c()
xaxis <- scale_x_continuous()
vertLines <- geom_vline()
speciesBreaks <- c()
speciesLabel <- c()
labelPos <- c()
chrBreaks <- c()
chrLabel <- c()
byChrLabel <- c()
byChrBreaks <- c()
byChrLabelBreaks <- c()
lineBreaks <- c()
sigSpecies <- c()
line <- geom_abline(intercept=0, slope=0)

chrs <- justChrList
alphaValue <- 1
median_value <- median(bedData$meanSkewLim)
#Reorder by chr
chrOrdered <- data.frame()
for (j in 1:length(justChrList)){
  chr_name = toString(chrs[[j]])
  if (j %% 2 == 0)  {
    alphaValue <- 0.99
  } else{
    alphaValue <- 1
  }
  for (k in 1:length(uniSpecies)) {
    speciesName <- toString(uniSpecies[[k]])
    spcChrName <- paste(speciesName, chr_name, sep="_")
    if (is.element(spcChrName, bedData$chrom)){
      #print(spcChrName)
      chr <- subset(bedData, chrom==spcChrName) 
      chr <- cbind(abs_count = chr$meanSkewLim/median_value, chr)
      chr <- cbind(alpha = alphaValue, chr)
      chr <- cbind(spc_Chr_Pos= seq(0, (length(chr$chrom)-1)*stepSize, stepSize), chr)
      chr <- cbind(species_name = speciesName, chr)
      chr <- cbind(chrName = chr_name, chr)
      chrOrdered <- rbind(chrOrdered, chr)      
    }
  }
}
chrOrdered <- cbind(Chr_Pos = seq(0, (length(chrOrdered$chrom)-1)*stepSize, stepSize), chrOrdered)
chrOrdered$abs_count[chrOrdered$abs_count==Inf] <- NA

legendValues = c()

pdf(NULL)
if (length(uniSpecies)>1) {
  colorList <- palette(rainbow(length(uniSpecies)))
  colorList <- palette(rainbow(length(uniSpecies)))
} else {
  colorList <- c("blue")
}
#Assign colors for each species
colors <- c()
for (k in 1:length(uniSpecies)){
  speciesName <- uniSpecies[[k]]
  color <- colorList[k]
  colors[speciesName] <- color
}

sigTable <- c()
#Cycle through Species and Chr
#pull out important species by if the spcies mean cov is greater than the overall mean
for (k in 1:length(uniSpecies)) {
  sigMarker = 0
  speciesName <- toString(uniSpecies[[k]])
  spc <- subset(chrOrdered, grepl(speciesName, chrom))
  subLen <- length(spc[,1]) 
  labelPos <- append(labelPos, spc[1,"Genome_Pos"])
  speciesBreaks <- append(speciesBreaks, spc[1,"Genome_Pos"])
  spcMaxMedian <- median(spc$meanSkewLim)
  spcMaxMean <- mean(spc$meanSkewLim)
  spc99Mean <- quantile(spc$meanSkewLim, 0.99)[[1]]
  if (is.na(spc99Mean)){
    spc99Mean <- 0
    spcMaxMedian <- 0
    spcMaxMean <- 0
  }
  sigSpecies <- c(sigSpecies, speciesName)
  #speciesTable <- cbind(abs_count = spc$meanSkewLim/median_value, spc)
  speciesTable <- cbind(spc_Chr_Pos= seq(0, (length(spc$chrom)-1)*stepSize, stepSize), spc)
  speciesTable <- cbind(species_name = speciesName, speciesTable)
  sigTable <- rbind(sigTable, speciesTable)
  if (spcMaxMedian>median) {
    write(paste(speciesName, "Max median avgDepth =", spcMaxMedian, sep=" "), file=covMeanFile, append=T)   
    median_value <- median(spc$meanSkewLim)
  }
  if (spcMaxMedian>maxMedian){
    maxMedian <- spcMaxMedian
  }
  if (spcMaxMean>mean) {
    write(paste(speciesName, "Max mean avgDepth =", spcMaxMean, sep=" "), file=covMeanFile, append=T)
  }
  if (spcMaxMean>maxMean){
    maxMean <- spcMaxMean
    meanQuant <- quantile(spc$meanSkewLim, 0.99, na.rm=T)
    logQuant <- quantile(spc$log2, 0.99, na.rm=T)
  }
  speciesLabel <- append(speciesLabel, speciesName)
  legendValues = c(legendValues, paste(speciesName, colorList[k]), sep="=")
  #go through by chr 
  counter=1
  for (j in 1:length(justChrList)){
    chr_name = toString(chrs[[j]])
    spcChrName <- paste(speciesName, chr_name, sep="_")
    if (k==1) {
      byChrLabel <- append(byChrLabel, chr_name)
    }    
    if (is.element(spcChrName, spc$chrom)){
      chr <- subset(spc, chrom==spcChrName)   
      chr.out <- data.frame("speciesChrName"=spcChrName, "covMean"=mean(chr$meanSkewLim), "covMedian"=median(chr$meanSkewLim), "relativeMean"=mean(chr$relativeMean))
      output <- rbind(output, chr.out)
      chrBreaks <- append(chrBreaks, chr[1,"Genome_Pos"])
      if (counter %% 2 == 0)  {
        alphaValue <- 0.9
      } else{
        alphaValue <- 1
      }
      chrPlotFill <- geom_ribbon(data=chr, aes(x=Chr_Pos, ymin=0, ymax=, fill=species_name))
      chrTotalFill <- c(chrTotalFill, chrPlotFill)
      chrPlotFillMean <- geom_ribbon(data=chr, aes(x=Chr_Pos, ymin=0, ymax=log2, fill=species_name))
      chrTotalFillMean <- c(chrTotalFillMean, chrPlotFillMean)
      counter = counter + 1      
    } else {
      #     print("Doesn't exist")
      counter = counter + 1
    }
  }
}
write(paste("Max median avgDepth =", maxMedian, sep=" "), file=covMeanFile, append=T)
write(paste("Max mean avgDepth =", maxMean, sep=" "), file=covMeanFile, append=T)
write(paste("99% maxAvgDepth =", meanQuant, sep=" "), file=covMeanFile, append=T)

for (i in 1:length(byChrLabel)) {
  chr.Name <- byChrLabel[i]
  chr <- subset(chrOrdered,  chrName==chr.Name)
  byChrBreaks <- append(byChrBreaks, chr[1,"Chr_Pos"])
  chrEnd <- chr[length(chr[,1]),1]
  midpoint <- chr[1,"Chr_Pos"]+((chrEnd-chr[1,"Chr_Pos"])/2)
  byChrLabelBreaks <- c(byChrLabelBreaks, midpoint)
}

chrOrdered$meanValueLimited <- chrOrdered$meanSkewLim
chrOrdered$meanValueLimited[chrOrdered$meanValueLimited>meanQuant] <- max(meanQuant)

totalFill <- geom_ribbon(data=chrOrdered, aes(x=Genome_Pos, ymin=0, ymax=meanSkewLim, fill=species_name))
totalFillLog <- geom_ribbon(data=chrOrdered, aes(x=Genome_Pos, ymin=0, ymax=log2, fill=species_name))
totalPoint <- geom_point(data=chrOrdered, aes(x=Genome_Pos, y=meanSkewLim, colour=species_name))
totalPointLog <- geom_point(data=chrOrdered, aes(x=Genome_Pos, y=log2, colour=species_name))
chrTotalPoint <- geom_point(data=chrOrdered, aes(x=Chr_Pos, y=meanSkewLim, colour=species_name))
chrTotalPointMean <- geom_point(data=chrOrdered, aes(x=Chr_Pos, y=log2, colour=species_name))

avgPlotLine <- geom_line(data=bedData, aes(x=Genome_Pos, y=meanSkewLim), alpha=0.5)
line <- geom_abline(intercept=0, slope=0)
endPos <- bedData[length(bedData[,1]),1]
xaxis <- scale_x_continuous(breaks=speciesBreaks, labels=speciesLabel, name="Genome Position", limits=c(0,endPos))
if (length(byChrLabel)==length(byChrLabelBreaks)){
  chrXaxis <- scale_x_continuous(breaks=byChrLabelBreaks, labels=byChrLabel, name="Genome Position", limits=c(0,endPos))
} else {
  chrXaxis <- scale_x_continuous(name="Genome Position", limits=c(0,endPos)) 
}
vertLines <- geom_vline(xintercept = speciesBreaks)
chrVertLines <- geom_vline(xintercept = byChrBreaks)
chrLines <- geom_segment(aes(x = chrBreaks, xend = chrBreaks, y = 0, yend=max(bedData$meanSkewLim)), alpha=0.8, linetype="dotted")
chrLegend <- scale_fill_manual(name="Species", values = speciesColors)
chrPointLegend <- scale_color_manual(name="Species", values=speciesColors, breaks=uniSpecies)
covline10 <- geom_hline(yintercept=10, colour="grey", linetype=2)
covline3 <- geom_hline(yintercept=3, colour="grey", linetype=2)
covline <- geom_hline(yintercept=lowCov, colour="black", linetype=2)

pdf(paste(outputPrefix, "_", window, "_depth.pdf", sep=''), width=14)

yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,NA))
plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage w/ coverage cutoffs", sep=" "))
plot(plot+plotTitle+xaxis+yaxis+totalFill+vertLines+chrLines+line+ theme_classic()+chrLegend+covline10+covline3+covline)
plot(plot+plotTitle+xaxis+yaxis+totalPoint+vertLines+chrLines+line+ theme_classic()+chrPointLegend+scale_alpha(guide = 'none')+covline10+covline3+covline)
plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage with mean Cov Bins", sep=" "))
plot(plot+plotTitle+xaxis+yaxis+totalFill+vertLines+chrLines+line+ theme_classic()+chrLegend+hBinAvgLines)
plot(plot+plotTitle+xaxis+yaxis+totalPoint+vertLines+chrLines+line+ theme_classic()+chrPointLegend+scale_alpha(guide = 'none')+hBinAvgLines)
chrLines <- geom_segment(aes(x = chrBreaks, xend = chrBreaks, y =0, yend=(max(logQuant)/4)), alpha=0.5, linetype="dotted")
yaxis <- scale_y_continuous(name="log2(avg/whole genome avg)", limits = c(-(max(logQuant)/100),max(logQuant)))
plotTitle <- ggtitle(paste(outputPrefix, "log2 Mean Avg depth of coverage", sep=" "))
plot(plot+plotTitle+xaxis+yaxis+totalFillLog+chrLines+vertLines+line+ theme_classic()+chrLegend+hBinLogLines)
plot(plot+plotTitle+xaxis+yaxis+totalPointLog+chrLines+vertLines+line+ theme_classic()+chrPointLegend+hBinLogLines)

#plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage sorted by Chr", sep=" "))
#plot(plot+plotTitle+chrLegend+chrXaxis+yaxis+chrTotalFill+chrVertLines+line+ theme_classic()+covline)
#plot(plot+plotTitle+chrXaxis+yaxis+chrTotalPoint+chrVertLines+line+chrPointLegend+theme_classic()+scale_alpha(guide = 'none')+covline)

#outputFileName <- paste(outputPrefix, window, "speciesChromosomeCoverage.txt", sep="_")
#write.table(output, file=outputFileName, row.names=F, sep = "\t")

#Work with just the species who's mean depth is above the absolute mean (Should I be picking these by median?)

sigTable <- cbind(sig_Gen_Pos = seq(0, (length(sigTable$chrom)-1)*stepSize, stepSize), sigTable)
sigTable$meanValueLimited <- sigTable$meanSkewLim
sigTable$meanValueLimited[sigTable$meanValueLimited>meanQuant] <- max(meanQuant)
#sigOutputFileName <- paste(outputPrefix, window, "sigSpeciesAvgCoverage.txt", sep="_")
#write.table(sigTable, file=sigOutputFileName, row.names=F, sep = "\t")

sigSpeciesColors <- c()
for (i in sigSpecies){
  index <- which(uniSpecies %in% i)
  sigSpeciesColors <- c(sigSpeciesColors, colorList[index])
}

xaxis <- scale_x_continuous()
vertLines <- geom_vline()
speciesBreaks <- c()
speciesLabel <- c()
labelPos <- c()
chrBreaks <- c()
byChrLabelBreaks <- c()
byChrBreaks <- c()
byChrLabel <- c()
lineBreaks <- c()
facetBreaks <- c()
sigFill <- c()
sigFillMean <- c()
sigPoint <- c()
sigPointMean <- c()
sigChrFill <- c()
sigChrFillMean <- c()
sigChrPoint <- c()
sigChrPointMean <- c()
sigPointGrid <- c()

alphaValue <- 1
#Reorder by chr
sigChrOrdered <- data.frame()
for (j in 1:(length(justChrList))){
  if (j %% 2 == 0)  {
    alphaValue <- 0.99
  } else{
    alphaValue <- 1
  }
  chr_name = toString(chrs[[j]])
  for (k in 1:length(uniSpecies)) {
    speciesName <- toString(uniSpecies[[k]])
    spcChrName <- paste(speciesName, chr_name, sep="_")
    if (is.element(spcChrName, sigTable$chrom)){
      #print(spcChrName)
      chr <- subset(sigTable, chrom==spcChrName)
      chr <- cbind(alpha = alphaValue, chr)
      chr <- cbind(Chr_Pos = seq(0, (length(chr$chrom)-1)*stepSize, stepSize), chr)
      chr <- cbind(chrName = chr_name, chr)
      sigChrOrdered <- rbind(sigChrOrdered, chr)      
    }
  }
}
sigChrOrdered <- cbind(Sig_Chr_Pos = seq(0, (length(sigChrOrdered$chrom)-1)*stepSize, stepSize), sigChrOrdered)

labelPos <- c()
speciesBreaks <- c()
speciesLabel <- c()
legendValues = c()
sigFill <- c()
sigFillMean <- c()
sigPoint <- c()
sigPointMean <- c()
sigChrFill <- c()
sigChrFillMean <- c()
sigChrPoint <- c()
sigChrPointMean <- c()
facetBreaks <- list()
byChrBreaks <- c()
byChrLabel <- c()
chrBreaks <- c()
facetDF <- data.frame()
#Go through just the species with coverage
for (k in 1:length(sigSpecies)) {
  sigMarker = 0
  speciesName <- toString(sigSpecies[[k]])
  spc <- subset(sigChrOrdered, grepl(speciesName, chrom))
  scp <- cbind(speciesChrPos = seq(0, (length(spc$chrom)-1)*stepSize, stepSize), spc)
  subLen <- length(spc[,1]) 
  labelPos <- append(labelPos, spc[1,"sig_Gen_Pos"])
  speciesBreaks <- append(speciesBreaks, spc[1,"sig_Gen_Pos"])
  speciesLabel <- append(speciesLabel, speciesName)
  legendValues = c(legendValues, speciesName=colorList[k])
  counter=1 
  for (j in 1:(length(justChrList))){
    chr_name = toString(chrs[[j]])
    spcChrName <- paste(speciesName, chr_name, sep="_")
    if (is.element(spcChrName, spc$chrom)){
      chr <- subset(spc, chrom==spcChrName)
      spcFacetBreaks <- data.frame("Breaks"= chr[1,"spc_Chr_Pos"], "species_name"= speciesName)     
      facetDF <- rbind(facetDF, spcFacetBreaks)
      if (k==1) {
        #print(chr[1,1])
        byChrBreaks <- append(byChrBreaks, chr[1,"Sig_Chr_Pos"])  
        byChrLabel <- append(byChrLabel, chr_name)
      }    
      chrBreaks <- append(chrBreaks, chr[1,"sig_Gen_Pos"])
      if (counter %% 2 == 0)  {
        alphaValue <- 0.95
      } else{
        alphaValue <- 1
      }
      sigPlotFill <- geom_ribbon(data=chr, aes(x=sig_Gen_Pos, ymin=0, ymax=meanValueLimited, fill=species_name))
      sigFill <- c(sigFill, sigPlotFill)
      sigPlotPoint <- geom_point(data=chr, aes(x=sig_Gen_Pos, y=meanValueLimited, colour=species_name))
      sigPoint <- c(sigPoint, sigPlotPoint)
      sigChrPlot <- geom_ribbon(data=chr, aes(x=Sig_Chr_Pos, ymin=0, ymax=meanValueLimited, fill=species_name))
      sigChrFill <- c(sigChrFill, sigChrPlot)
      sigChrPlotPoint <- geom_point(data=chr, aes(x=Sig_Chr_Pos, y=meanValueLimited, colour=species_name))
      sigChrPoint <- c(sigChrPoint, sigChrPlotPoint)
      chrEnd <- chr[length(chr[,1]),1]
      midpoint <- chr[1,1]+((chrEnd-chr[1,1])/2)
    } else {
      #     print("Doesn't exist")    
    }
    counter = counter + 1
  }
}

byChrLabelBreaks <- c()
for (i in 1:length(byChrLabel)) {
  chr.Name <- byChrLabel[i]
  chr <- subset(sigChrOrdered,  chr.Name==chrName)
  chrEnd <- chr[length(chr[,1]),1]
  midpoint <- chr[1,1]+((chrEnd-chr[1,1])/2)
  byChrLabelBreaks <- c(byChrLabelBreaks, midpoint)
}

avgPlotLine <- geom_line(data=sigTable, aes(x=sig_Gen_Pos, y=meanSkewLim), alpha=0.5)

line <- geom_abline(intercept=0, slope=0)
endPos <- sigTable[length(sigTable[,1]),1]
xaxis <- scale_x_continuous(breaks=speciesBreaks, labels=speciesLabel, name="Genome Position", limits=c(0,endPos))
chrXaxis <- scale_x_continuous(breaks=byChrLabelBreaks, labels=byChrLabel, name="Genome Position", limits=c(0,endPos))
vertLines <- geom_vline(xintercept = speciesBreaks)
chrVertLines <- geom_vline(xintercept = byChrBreaks)
chrLines <- geom_segment(aes(x = chrBreaks, xend = chrBreaks, y = 0, yend=quantile(sigChrOrdered$meanSkewLim, 0.99, na.rm=T)), alpha=0.5, linetype="dotted")
chrLegend <- scale_fill_manual(name="Species", values = legendValues)
chrPointLegend <- scale_color_manual(name="Species", values=colorList, breaks=sigSpecies)
chrColors <-palette(rainbow(length(justChrList)))
chrColors <-palette(rainbow(length(justChrList)))
yaxis <- scale_y_continuous(name="Average Depth")
yaxisFacet <- scale_y_continuous(name="Abs Count")
facetSubset <- subset(facetDF, species_name==sigSpecies[length(sigSpecies)])
facetLabelBreaks <- c()
for (i in 1:length(byChrLabel)) {
  chrStart <- facetSubset[i, "Breaks"]
  if (i==length(byChrLabel)){
    chrEnd <- 12000000
  }
  else {
    chrEnd <- facetSubset[i+1, "Breaks"] 
  }
  midpoint <- chrStart+((chrEnd-chrStart)/2)
  facetLabelBreaks <- c(facetLabelBreaks, midpoint)
}
xaxisFacet <- scale_x_continuous(breaks=facetLabelBreaks, labels=byChrLabel, name="Genome Position", limits=c(0,sigChrOrdered[length(sigChrOrdered$chrom),"spc_Chr_Pos"]))
facetVertLines <- geom_vline(aes(xintercept = Breaks), facetDF, linetype=2)

#pdf(paste(outputPrefix, "_choosenSpcies.pdf", sep=''), width=14)

#if (!is.na(max(sigChrOrdered$abs_count))) {
#  ggplot(sigChrOrdered, aes(spc_Chr_Pos, abs_count, colour = chrName))+geom_point()+facet_grid(species_name ~ .)+theme_classic()+line+scale_colour_manual(values=chrColors)+xaxisFacet+facetVertLines+yaxisFacet+plotTitle+covline10+covline3+covline
#}
plotTitle <- ggtitle(paste(outputPrefix, " Depth w/ cov cutoff", sep=" "))
ggplot(sigChrOrdered, aes(spc_Chr_Pos, meanValueLimited, colour = species_name))+geom_point()+facet_grid(species_name ~ .)+theme_classic()+line+scale_colour_manual(values=speciesColors, name="Species")+covline10+covline3+covline+facetVertLines+yaxis+plotTitle+xaxisFacet+scale_alpha(guide = 'none')
ggplot()+geom_ribbon(data=sigChrOrdered, aes(x=spc_Chr_Pos, ymin=0, ymax=meanValueLimited, fill = species_name))+theme_classic()+line+scale_fill_manual(values=speciesColors, name="Species")+facetVertLines+yaxis+plotTitle+facet_grid(species_name ~ .)+xaxisFacet+covline10+covline3+covline
plotTitle <- ggtitle(paste(outputPrefix, " Depth w/ cov bins", sep=" "))
ggplot(sigChrOrdered, aes(spc_Chr_Pos, meanValueLimited, colour = species_name))+geom_point()+facet_grid(species_name ~ .)+theme_classic()+line+scale_colour_manual(values=speciesColors, name="Species")+facetVertLines+yaxis+plotTitle+xaxisFacet+scale_alpha(guide = 'none')+ hBinAvgLines
ggplot()+geom_ribbon(data=sigChrOrdered, aes(x=spc_Chr_Pos, ymin=0, ymax=meanValueLimited, fill = species_name))+theme_classic()+line+scale_fill_manual(values=speciesColors, name="Species")+facetVertLines+yaxis+plotTitle+facet_grid(species_name ~ .)+xaxisFacet+ hBinAvgLines
yaxis <- scale_y_continuous(name="Log2 Average Depth")
plotTitle <- ggtitle(paste(outputPrefix, " Log2 Depth w/ cov bins", sep=" "))
ggplot(sigChrOrdered, aes(spc_Chr_Pos, log2, colour = species_name))+geom_point()+facet_grid(species_name ~ .)+theme_classic()+line+scale_colour_manual(values=speciesColors, name="Species")+facetVertLines+yaxis+plotTitle+xaxisFacet+scale_alpha(guide = 'none')+ hBinLogLines
ggplot()+geom_ribbon(data=sigChrOrdered, aes(x=spc_Chr_Pos, ymin=0, ymax=log2, fill = species_name))+theme_classic()+line+scale_fill_manual(values=speciesColors, name="Species")+facetVertLines+yaxis+plotTitle+facet_grid(species_name ~ .)+xaxisFacet+ hBinLogLines

dev.off()

# yaxisFacet <- scale_y_continuous(name="Avg Depth")
# ggplot(sigChrOrdered, aes(spc_Chr_Pos, meanValue, colour = species_name))+geom_point()+facet_grid(species_name ~ .)+theme_classic()+line+scale_colour_manual(values=speciesColors, name="Species")+covline10+covline3+covline+facetVertLines+yaxis+plotTitle+xaxisFacet+scale_alpha(guide = 'none')+ hBinAvgLines
# ggplot()+geom_ribbon(data=sigChrOrdered, aes(x=spc_Chr_Pos, ymin=0, ymax=meanValue, fill = species_name))+theme_classic()+line+scale_fill_manual(values=speciesColors, name="Species")+facetVertLines+yaxis+plotTitle+facet_grid(species_name ~ .)+xaxisFacet+covline10+covline3+covline+ hBinAvgLines
# ggplot()+geom_ribbon(data=sigChrOrdered, aes(x=spc_Chr_Pos, ymin=0, ymax=meanValue, fill = chrName))+theme_classic()+line+scale_fill_manual(values=chrColors, name="Chromosome")+facetVertLines+yaxis+plotTitle+facet_grid(species_name ~ .)+xaxisFacet+covline10+covline3+covline+ hBinAvgLines
# xaxisFacet2 <- scale_x_continuous(breaks=seq(0,endPos,100000), labels=NULL, name="Chromosome Position")
# ggplot(sigChrOrdered, aes(Chr_Pos, meanValue, colour = species_name))+geom_point()+facet_grid(species_name ~ chrName,  scales = "free", space = "free")+theme_classic()+line+scale_colour_manual(values=speciesColors, name="Species")+yaxis+plotTitle+scale_alpha(guide = 'none')+xaxisFacet2+covline10+covline3+covline+ hBinAvgLines
# ggplot()+geom_ribbon(data=sigChrOrdered, aes(x=Chr_Pos, ymin=0, ymax=meanValue, fill = species_name))+facet_grid(species_name ~ chrName,  scales = "free", space = "free")+theme_classic()+line+scale_fill_manual(values=speciesColors, name="Species")+yaxis+plotTitle+scale_alpha(guide = 'none')+xaxisFacet2+covline10+covline3+covline
# ggplot()+geom_ribbon(data=sigChrOrdered, aes(x=Chr_Pos, ymin=0, ymax=meanValue, fill = chrName))+facet_grid(species_name ~ chrName,  scales = "free", space = "free")+theme_classic()+line+scale_colour_manual(values=chrColors, name="Chromosome")+yaxis+plotTitle+scale_alpha(guide = 'none')+xaxisFacet2+covline10+covline3+covline+ hBinAvgLines

# plot <- ggplot()
# plot(plot+plotTitle+xaxis+yaxis+sigFill+vertLines+chrLines+line+ theme_classic()+scale_fill_manual(values=colors, name="Species"))
# plot(plot+plotTitle+xaxis+yaxis+sigPoint+vertLines+chrLines+line+ theme_classic()+scale_colour_manual(values=colors, name="Species"))
# plot(plot+plotTitle+chrXaxis+yaxis+chrVertLines+line+scale_fill_manual(values=colors, name="Species")+ theme_classic()+sigChrFill)
# plot(plot+plotTitle+chrXaxis+yaxis+chrVertLines+line+scale_colour_manual(values=colors, name="Species")+theme_classic()+sigChrPoint)
