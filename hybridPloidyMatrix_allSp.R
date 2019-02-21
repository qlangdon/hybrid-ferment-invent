options(stringAsFactors=FALSE)
args <- commandArgs(TRUE)
hybSpName <- args[1]
covBinInfo <- read.table(paste(hybSpName, "_covBinKey.txt", sep=""), header=T)
spChrWinKey <- read.table("ScerSkudSuvaSeub_10kbWinKey.txt", header=T)
strainList <- covBinInfo$strain
allSp <- unique(spChrWinKey$species)

outData <- data.frame(spChrPos=spChrWinKey$spChrPos)
colName <- c("spChrPos")
for (strainName in strainList) {
  strainCovInfoFull <- covBinInfo[which(covBinInfo$strain==strainName),2:ncol(covBinInfo)]
  strainBinKey <- data.frame()
  for (i in 1:ncol(strainCovInfoFull)) {
    binName <- strainCovInfoFull[,i]
    if (!is.na(binName)) {
      key <- data.frame(ploidy=i, binName=strainCovInfoFull[,i])
      strainBinKey <- rbind(strainBinKey, key)
    }
  }
  uniBins <- unique(strainBinKey$binName, rm.na=T)
  strainDepth <- read.table(paste(strainName, "_winAvgDepth_wMeanCovBins.txt", sep=""), header=T)
  uniStSp <- unique(strainDepth$species)
  strainPloidy <- c()
  for (spName in allSp) {
    spWinKey <- spChrWinKey[which(spChrWinKey$species==spName),]
    if (spName %in% uniStSp) {
      strainSpPloidy <- c()
      for (j in 1:nrow(spWinKey)) {
        chrName <- as.character(spWinKey$chrom[j])
        startPos <- spWinKey$chrStart[j]
        binInfo <- as.character(strainDepth$covBin[which(strainDepth$chrom==chrName & strainDepth$start==startPos)])
        if (binInfo %in% uniBins) {
          strainSpPloidy <- c(strainSpPloidy, strainBinKey$ploidy[which(strainBinKey$binName==binInfo)])
        } else {
          strainSpPloidy <- c(strainSpPloidy, 0)
        }
      }
    } else {
      strainSpPloidy <- rep.int(0, nrow(spWinKey))
    }
    strainPloidy <- c(strainPloidy, strainSpPloidy)
  }
  outData <- cbind(outData, strainPloidy)
  colName <- c(colName, strainName)
}
colnames(outData) <- colName
write.table(outData, paste(hybSpName, "winPloidy10kb_ScerSkudSuvaSeub.txt", sep="_"), row.names = F, quote=F)

