

args <- commandArgs(TRUE)
#regionsInfo <- read.table(args[1], check.names = F)
regionsInfo <- read.table("Skud_winPloidy10kb.txt", header=T)
#Use Info Table
chrLengths <- read.table("../SkudChrLengths_wBin.txt", header = T)

uniChr <- unique(chrLengths$chrName)

infoCol <- regionsInfo$spChrPos
wantedStrains <- data.frame(spChrPos=infoCol, AMH=regionsInfo[,grep(("AMH"),names(regionsInfo))], EP2=regionsInfo[,grep(("EP2"),names(regionsInfo))], FM64=regionsInfo[,grep(("FM64"),names(regionsInfo))], NT50=regionsInfo[,grep(("NT50"),names(regionsInfo))], yHQL557=regionsInfo[,grep(("yHQL557"),names(regionsInfo))], yHQL560=regionsInfo[,grep(("yHQL560"),names(regionsInfo))])
wData <- wantedStrains[apply(wantedStrains[,-1], 1, function(x) !all(x==0)),]
#chrInfo <- strsplit(as.character(wData$`regionsInfo$spChrPo`), "_")[[2]]
split <- strsplit(as.character(wData$spChrPo), "-")
wData$chrInfo <- unlist(lapply(split, function(x) x[1]))
wData$start <- as.numeric(unlist(lapply(split, function(x) x[2])))
wData$end <- as.numeric(unlist(lapply(split, function(x) x[3])))
chrSplit <- strsplit(as.character(wData$chrInfo ), "_")
wData$chrName <- unlist(lapply(chrSplit, function(x) x[2]))


outData <- data.frame(loci=as.character(), genomeStart=as.numeric(), genomeEnd=as.numeric(), strains=as.character())
for (c in 1:(length(uniChr))) {
  chr <- uniChr[c]
  chrLenInfo <- chrLengths[which(chrLengths$chrName==chr),] 
  chrBinS <- chrLenInfo$binStart4way
  chrData <- wData[which(wData$chrName==chr),]
  if (nrow(chrData)>0) {
    #print(paste(as.character(chr), nrow(chrData), sep=" "))
    chrOut <- data.frame(loci=as.character(), genomeStart=as.numeric(), genomeEnd=as.numeric(), strains=as.character())
    chrData$chrStart <- chrData$start-chrBinS
    chrData$chrEnd <- chrData$end-chrBinS
    for (i in 1:nrow(chrData)) {
      binData <- chrData[i,]
      colNums <- which(binData[1,2:(ncol(binData)-6)]>0)
      strainInfo <- "strains:"
      for (colNum in colNums) {
        strainInfo <- paste(strainInfo, colnames(binData[colNum+1]), sep=" ")
      }
      lociName <- paste(binData$chrInfo, "_", (binData$chrStart/1000), "kb-", (binData$chrEnd/1000), "kb", sep="")
      gStart <- binData$chrStart+chrLenInfo$globalStart
      gEnd <- binData$chrEnd+chrLenInfo$globalStart
      if (gEnd > chrLenInfo$globalEnd) {
        gEnd <- chrLenInfo$globalEnd
      }
      binInfo <- data.frame(lociName = lociName, genomeStart = gStart, genomeEnd = gEnd, strains = strainInfo)
      chrOut <- rbind(chrOut, binInfo)
    }
    outData <- rbind(outData, chrOut)
  }
  #print(paste("Output:", as.character(chr), nrow(chrOut), sep=" "))
}

write.table(outData, "Skud_introRegions.txt", row.names = F, quote=F)
