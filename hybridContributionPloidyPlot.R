options(scipen=999)
library("ggplot2")

setwd("~/Documents/Hittinger_Lab/whole_genome/brewingProject/contributionPlots")
#args <- commandArgs(TRUE)
ploidyFileName <- "ScerSkudSuvaSeub_ploidyInfo.txt"
chrLengthsFile <- "ScerSkudSuvaSeub_binnedChrLengths.txt"
strainOrderFile <- "strainOrder_contributionPlot.txt"

chrLengths <- read.table(chrLengthsFile, header=T)
strainOrder <- read.table(strainOrderFile, header=T)

speciesColors <- c("Scer"="red", "Spar"="#FFDB00", "Smik"="#49FF00", "Skud"="#00FF92", "SkudZP"="#00FF92", "SkudIFO"="#00FF92", "Sarb"="#0092FF", "Suva"="#4900FF", "Seub"="#FF00DB")

ploidyTable <- read.table(ploidyFileName, header=T)
maxPloidy <- max(ploidyTable$ploidy)
ploidyTable$yStart <- NA
ploidyTable$yEnd <- NA
uniStrains <- unique(strainOrder$genomeName)
uniSpecies <- unique(ploidyTable$species)

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in uniStrains) {
  strainData <- ploidyTable[which(ploidyTable$strainID==strainName),]
  strainMaxPloidy <- max(strainData$ploidy)
  #print(strainName)
  spMaxY <- strainYstart+strainMaxPloidy
  spLabelY <- strainYstart+(strainMaxPloidy/2)
  yBreaks <- c(yBreaks, spMaxY)
  yLabels <- c(yLabels, spLabelY)
  strainData$yStart <- strainYstart
  strainData$yEnd <- strainData$ploidy+strainYstart
  #strainData$uniStrainID <- strainOrder$newName[which(strainOrder$genomeID==strainName)]
  #ploidyTable$yStart[which(ploidyTable$strainID==strainName)] <- strainYstart
  #ploidyTable$yEnd[which(ploidyTable$strainID==strainName)] <- ploidyTable$ploidy[which(ploidyTable$strainID==strainName)]+strainYstart
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}

xBreaks <- c()
start <- 0
for (i in 1:nrow(chrLengths)){
  chrEnd <- chrLengths$endWhenBin[i]
  halfPoint <- start+(chrEnd - start)/2
  xBreaks <- c(xBreaks, halfPoint)
  start <- chrEnd
}

xaxis <- scale_x_continuous(name = "Genome Position", limits = c(0,chrLengths$endWhenBin[length((chrLengths$endWhenBin))]), breaks = c(6125000, 17890000, 29315000, 41000000), labels = c("Scer", "Skud", "Suva","Seub"), expand = c(0, 0))
yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
chrLines <- geom_vline(xintercept = c(0,chrLengths$endWhenBin), linetype="dashed", alpha=0.4)
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
#colorScale <- scale_color_manual(name="Species", values = speciesColors, breaks=uniSpecies)

pdf("allHyrbids_strainClustered_Genome_Contributions.pdf", width=11)

plotTitle <- ggtitle("All hybrids included in brewing project Hybrid Genome Contributions")
ggplot()+geom_rect(data=plotTable, aes(xmin = newGenoStart, xmax = newGenoEnd, ymin =yStart, ymax=yEnd, fill=species))+fillScale+xaxis+yaxis+chrLines+strainLines+theme_classic()+plotTitle+theme(axis.text.y = element_blank())

dev.off()


####ploidy by Alpha

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in uniStrains) {
  strainData <- ploidyTable[which(ploidyTable$strainID==strainName),]
  strainMaxPloidy <- max(strainData$ploidy)
  #print(strainName)
  spMaxY <- strainYstart+1
  spLabelY <- strainYstart+0.5
  yBreaks <- c(yBreaks, spMaxY)
  yLabels <- c(yLabels, spLabelY)
  strainData$yStart <- strainYstart
  strainData$yEnd <- 1+strainYstart
  strainData$alpha <- 1
  if (strainMaxPloidy == 4){
    strainData$alpha[which(strainData$ploidy==4)] <- 1
    strainData$alpha[which(strainData$ploidy==3)] <- 0.88
    strainData$alpha[which(strainData$ploidy==2)] <- 0.76
    strainData$alpha[which(strainData$ploidy==1)] <- 0.64
  } else if (strainMaxPloidy==3) {
    strainData$alpha[which(strainData$ploidy==3)] <- 0.88
    strainData$alpha[which(strainData$ploidy==2)] <- 0.76
    strainData$alpha[which(strainData$ploidy==1)] <- 0.64
  } else if (strainMaxPloidy==2) {
    strainData$alpha[which(strainData$ploidy==2)] <- 0.76
    strainData$alpha[which(strainData$ploidy==1)] <- 0.64
  }
  #strainData$uniStrainID <- strainOrder$newName[which(strainOrder$genomeID==strainName)]
  #ploidyTable$yStart[which(ploidyTable$strainID==strainName)] <- strainYstart
  #ploidyTable$yEnd[which(ploidyTable$strainID==strainName)] <- ploidyTable$ploidy[which(ploidyTable$strainID==strainName)]+strainYstart
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}

xBreaks <- c()
start <- 0
for (i in 1:nrow(chrLengths)){
  chrEnd <- chrLengths$endWhenBin[i]
  halfPoint <- start+(chrEnd - start)/2
  xBreaks <- c(xBreaks, halfPoint)
  start <- chrEnd
}

xaxis <- scale_x_continuous(name = "Genome Position", limits = c(0,chrLengths$endWhenBin[length((chrLengths$endWhenBin))]), breaks = c(6125000, 17890000, 29315000, 41000000), labels = c("Scer", "Skud", "Suva","Seub"), expand = c(0, 0))
yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
chrLines <- geom_vline(xintercept = c(0,chrLengths$endWhenBin), linetype="dashed", alpha=0.3)
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
#colorScale <- scale_color_manual(name="Species", values = speciesColors, breaks=uniSpecies)

pdf("allHyrbids_ploidyAlpha_genomeContributions.pdf", width=11)

plotTitle <- ggtitle("All hybrids included in brewing project Hybrid Genome Contributions")
ggplot()+geom_rect(data=plotTable, aes(xmin = newGenoStart, xmax = newGenoEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=plotTable$alpha)+fillScale+xaxis+yaxis+chrLines+strainLines+theme_classic()+plotTitle+theme(axis.text.y = element_blank(), axis.text.x=element_text(face="italic"))

dev.off()

###Contirbution Summaries
ploidyTable$ploidyXregLen <- ploidyTable$ploidy * ploidyTable$regionLength
strainSummary <- data.frame()
for (strainName in uniStrains){
  tempTable <- data.frame(strainID=strainName)
  strainData <- ploidyTable[which(ploidyTable$strainID==strainName),]
  for (speciesName in uniSpecies) {
    spData <- strainData[which(strainData$species==speciesName),]
    spSum <- sum(spData$ploidyXregLen)
    tempTable <- cbind(tempTable, spSum)
  }
  colnames(tempTable) <- c("strainID", "ScerCont", "SkudZPCont", "SuvaCont", "SeubCont")
  totalSum <- sum(tempTable$ScerCont, tempTable$SkudZPCont, tempTable$SuvaCont, tempTable$SeubCont)
  tempTable$totalContent <- totalSum
  strainSummary <- rbind(strainSummary, tempTable)
}
strainSummary$ScerNorm <- strainSummary$ScerCont/strainSummary$totalContent
strainSummary$SkudNorm <- strainSummary$SkudZPCont/strainSummary$totalContent
strainSummary$SuvaNorm <- strainSummary$SuvaCont/strainSummary$totalContent
strainSummary$SeubNorm <- strainSummary$SeubCont/strainSummary$totalContent
strainSummary$ScerByGenome <- strainSummary$ScerCont/12000000
strainSummary$SkudByGenome <- strainSummary$SkudZPCont/12000000
strainSummary$SuvaByGenome <- strainSummary$SuvaCont/12000000
strainSummary$SeubByGenome <- strainSummary$SeubCont/12000000
write.table(format(strainSummary, scientific=FALSE), file="hybridStrainContributionSummary.txt", row.names=F, sep = "\t", quote = FALSE)


####For Species

xaxis <- scale_x_continuous(name = "Genome Position", limits = c(0,chrLengths$endWhenBin[length((chrLengths$endWhenBin))]), breaks = xBreaks, labels = chrLengths$chrName, expand = c(0, 0))
yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, labels=strainOrder$newName, expand = c(0, 0))
chrLines <- geom_vline(xintercept = c(0,chrLengths$endWhenBin), linetype="dashed")
strainLines <- geom_hline(yintercept = yBreaks)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
#colorScale <- scale_color_manual(name="Species", values = speciesColors, breaks=uniSpecies)

pdf(paste(as.character(uniSpecies), "renamedOtherHyb_Hybrid_Genome_Contributions.pdf", sep="_"), width=11)

plotTitle <- ggtitle(paste(as.character(uniSpecies), "other hybrids included in brewing project Hybrid Genome Contributions", sep=" "))
ggplot()+geom_rect(data=plotTable, aes(xmin = spStart, xmax = spEnd, ymin =yStart, ymax=yEnd, fill=species))+
  fillScale+xaxis+yaxis+chrLines+strainLines+theme_classic()+plotTitle

dev.off()
