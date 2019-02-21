options(scipen=999)
args <- commandArgs(TRUE)

strainList <- read.table(args[1], header=F)
chrLengths <- read.table(args[2], header=T)
#chrSplit <- strsplit(as.character(chrLengths$chrName), "_")
#chrLengths$species <- unlist(lapply(chrSplit, function(x) x[1]))
uniSpecies <- unique(chrLengths$species)

hybName <- ""
#spDFs <- list()
for (i in 1:length(uniSpecies)) {
  hybName <- paste(hybName, "X", uniSpecies[i], sep="")
}
hybName <- substr(hybName, 2, nchar(hybName))

for (i in 1:length(uniSpecies)) {
  fileOutName <- paste(hybName, uniSpecies[i], "ploidyInfo.txt", sep="_")
  tempData <- data.frame("strainID"=as.character(), "chrName"=as.character(), "ploidy"=as.character(), "genomeStart"=as.numeric(), "genomeEnd"=as.numeric(), "chrStart"=as.numeric(), "chrEnd"=as.numeric(), "regionLength"=as.numeric(), "species"=as.character(), "chrNum"=as.character(), "spStart"=as.numeric(), "spEnd"=as.numeric())
  write.table(tempData, file=fileOutName, row.names = F, sep="\t", quote=F)
}
#spDFs[i] <- data.frame("strainID"=as.character(), "chrName"=as.character(), "ploidy"=as.character(), "genomeStart"=as.numeric(), "genomeEnd"=as.numeric(), "chrStart"=as.numeric(), "chrEnd"=as.numeric(), "regionLength"=as.numeric(), "species"=as.character, "chrNum"=as.character, "spStart"=as.numeric(), "spEnd"=as.numeric())

for (j in (1:nrow(strainList))) {
  strainName <- strainList[j,1]
  strainFileName <- paste(as.character(strainName), "_min10KB_ploidyPositions.txt", sep="")
  if (file.exists(strainFileName)) {
    strainData <- read.table(strainFileName, header=T)
    chrSplit <- strsplit(as.character(strainData$chrName), "_")
    strainData$species <- unlist(lapply(chrSplit, function(x) x[1]))
    strainData$chrNum <- unlist(lapply(chrSplit, function(x) x[2]))
    for (i in 1:length(uniSpecies)){
      fileOutName <- paste(hybName, uniSpecies[i], "ploidyInfo.txt", sep="_")
      speciesName <- uniSpecies[i]
      spChrLens <- chrLengths[which(chrLengths$species==speciesName),]
      spGlobalStart <- spChrLens$spStart[1]
      strainSpData <- strainData[which(strainData$species==speciesName),]
      strainSpData$spStart <- strainSpData$genomeStart-spGlobalStart
      strainSpData$spEnd <- strainSpData$genomeEnd-spGlobalStart
      write.table(strainSpData, file=fileOutName, col.names=F, row.names = F, sep="\t", quote=F, append=T)
      #spDFs[i] <- rbind(spDFs[i], strainSpData)
    }
  } else {
    print(paste("Cannot find", strainFileName, sep=" "))
  }
}



