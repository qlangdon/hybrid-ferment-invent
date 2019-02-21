options(scipen=999)
library("ggplot2")

speciesColors <- c("Scer"="red", "Spar"="#FFDB00", "Smik"="#49FF00", "Skud"="#00FF92", "SkudZP"="#00FF92", "SkudIFO"="#00FF92", "Sarb"="#0092FF", "Suva"="#4900FF", "Seub"="#FF00DB")

setwd("~/Documents/Hittinger_Lab/whole_genome/brewingProject/contributionPlots")
#args <- commandArgs(TRUE)
ploidyFileName <- "ScerSkudSuvaSeub_ploidyInfo.txt"
chrLengthsFile <- "ScerSuvaSeub_binnedChrLengths.txt"
chrLengths <- read.table(chrLengthsFile, header=T)

strainsWanted <- c("yHAB34_ScerXSuvaXSeub", "Muri_ScerXSuvaXSeub_krogerus")

ploidyTable <- read.table(ploidyFileName, header=T)
maxPloidy <- max(ploidyTable$ploidy)
ploidyTable$yStart <- NA
ploidyTable$yEnd <- NA

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in strainsWanted) {
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

xaxis <- scale_x_continuous(name = "Genome Position", breaks = c(6125000, 18130000, 29810000), labels = c("Scer", "Suva","Seub"), expand = c(0, 0))
#yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
#chrLines <- geom_vline(xintercept = c(0,chrLengths$endWhenBin), linetype="dashed", alpha=0.3)
#strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

spLines <- geom_vline(xintercept=c(12250000,23820000))
fillScale <- scale_fill_manual(values=speciesColors, name="Species", guide=FALSE)
#colorScale <- scale_color_manual(name="Species", values = speciesColors, breaks=uniSpecies)

plotTabMod <- plotTable
plotTabMod$strainID <- as.character(plotTabMod$strainID)
plotTabMod$strainID[which(plotTabMod$strainID=="yHAB34_ScerXSuvaXSeub")] <- "WLP351"
plotTabMod$strainID[which(plotTabMod$strainID=="Muri_ScerXSuvaXSeub_krogerus")] <- "Muri"

ggplot()+geom_rect(data=plotTabMod, aes(xmin = genomeStart, xmax = genomeEnd, ymin =0, ymax=ploidy, fill=species))+fillScale+facet_grid(strainID~.)+ylab("Ploidy")+spLines+xaxis+theme_classic()+theme(axis.text.x=element_text(face="italic"))
ggsave("Muri-WLP351.pdf", height=8, width=12)

pdf("allHyrbids_ploidyAlpha_genomeContributions.pdf", width=11)

plotTitle <- ggtitle("All hybrids included in brewing project Hybrid Genome Contributions")
ggplot()+geom_rect(data=plotTable, aes(xmin = genomeStart, xmax = genomeEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=plotTable$alpha)+fillScale+xaxis+yaxis+chrLines+strainLines+theme_classic()+plotTitle+theme(axis.text.y = element_blank(), axis.text.x=element_text(face="italic"))

dev.off()

