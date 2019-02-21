options(scipen=999)
library("ggplot2")

setwd("~/Documents/Hittinger_Lab/whole_genome/brewingProject/contributionPlots")
#args <- commandArgs(TRUE)
ploidyFileName <- "ScerSkudSuvaSeub_ploidyInfo.txt"

speciesColors <- c("Scer"="red", "Spar"="#FFDB00", "Smik"="#49FF00", "Skud"="#00FF92", "SkudZP"="#00FF92", "SkudIFO"="#00FF92", "Sarb"="#0092FF", "Suva"="#4900FF", "Seub"="#FF00DB")

ploidyTable <- read.table(ploidyFileName, header=T)
maxPloidy <- max(ploidyTable$ploidy)
ploidyTable$yStart <- NA
ploidyTable$yEnd <- NA
uniSpecies <- unique(ploidyTable$species)


###PAD1 & FDC1 ###
setwd("~/Documents/Hittinger_Lab/whole_genome/brewingProject/genesOfInterest/")
PAD1FDC1geneLoc <- read.table("../genesOfInterest/PAD1-FDC1_geneLoc.txt", header=T)
PAD1FDC1geneLoc$chrName <- PAD1FDC1geneLoc$chrom
PFchrLens <- read.table("../genesOfInterest/PAD1-FDC1_chrLens.txt", header=T)
PFends <- PFchrLens[which(PFchrLens$ends>0),]
PFends$half <- PFends$ends-100000
ScerHalf <- PFends$half[which(PFends$chrName=="Scer_chrIV")]
SkudHalf <- PFends$half[which(PFends$chrName=="SkudZP_chrIV")]
SuvaHalf <- PFends$half[which(PFends$chrName=="Suva_chrII")]
SeubHalf <- PFends$half[which(PFends$chrName=="Seub_chrXIII")]

##Subset#
SuvaXSeubList <- read.table("SuvaXSeub_wPAD1FDC1.txt", header=T)
ueStrains <- unique(SuvaXSeubList$genomeName)
PFwant <- PFends[which(PFends$chrName=="Suva_chrII" | PFends$chrName=="Seub_chrXIII"),]

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in ueStrains) {
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
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}
toPlotTemp <- plotTable[which(plotTable$chrName=="Suva_chrII" | plotTable$chrName=="Seub_chrXIII"),]
toPlotTemp$newStart <- toPlotTemp$chrStart
toPlotTemp$newStart[which(toPlotTemp$species=="Suva" & toPlotTemp$chrStart<SuvaHalf & toPlotTemp$chrEnd>SuvaHalf)] <- SuvaHalf+1
toPlotTemp$newStart[which(toPlotTemp$species=="Seub" & toPlotTemp$chrStart<SeubHalf & toPlotTemp$chrEnd>SeubHalf)] <- SeubHalf+1
toPlot <- rbind(toPlotTemp[which(toPlotTemp$species=="Suva" & toPlotTemp$newStart>SuvaHalf),], toPlotTemp[which(toPlotTemp$species=="Seub" & toPlotTemp$newStart>SeubHalf),])

yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Suva" | PAD1FDC1geneLoc$species=="Seub"),]
genePos$chr <- factor(as.character(genePos$chrom), levels=c("Suva_chrII", "Seub_chrXIII"))
PFwant$chr <- factor(as.character(PFwant$chrName), levels=c("Suva_chrII", "Seub_chrXIII"))
toPlot$chr <- factor(as.character(toPlot$chrName), levels=c("Suva_chrII", "Seub_chrXIII"))

ggplot()+geom_rect(data=toPlot, aes(xmin = newStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlot$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends))+geom_vline(data=genePos, aes(xintercept = start), linetype="dotted") +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())+geom_vline(data=PFwant, aes(xintercept = half), linetype="dashed", color="grey")

toPlotMod <- toPlot
toPlotMod$newStart[which(toPlotMod$strainID=="yHQL580_SuvaXSeub" & toPlotMod$newStart==940000)] <- 930000
toPlotMod <- toPlotMod[-which(toPlotMod$strainID=="yHQL561_SuvaXSeub" & toPlotMod$newStart==950000),]
toPlotMod <- toPlotMod[-which(toPlotMod$strainID=="yHQL559_SuvaXSeub" & toPlotMod$newStart==950000),]
toPlotMod <- toPlotMod[-which(toPlotMod$strainID=="yHQL572_SuvaXSeub" & toPlotMod$newStart==950000),]
toPlotMod <- toPlotMod[-which(toPlotMod$strainID=="CBS380_SuvaXSeub_okuno" & toPlotMod$newStart==950000),]
toPlotMod <- toPlotMod[-which(toPlotMod$strainID=="yHQL522_SuvaXSeub_plate11" & toPlotMod$newStart==950000),]
toPlotMod <- toPlotMod[-which(toPlotMod$strainID=="yHQL555_SuvaXSeub" & toPlotMod$newStart==950000),]
toPlotMod <- toPlotMod[-which(toPlotMod$strainID=="yHQL570_SuvaXSeub" & toPlotMod$newStart==950000),]

ggplot()+geom_rect(data=toPlotMod, aes(xmin = newStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlotMod$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends))+geom_vline(data=genePos, aes(xintercept = start), linetype="dotted") +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())+geom_vline(data=PFwant, aes(xintercept = half), linetype="dashed", color="grey")

ggsave("SuvaXSeub_endChrPAD1FDC1.jpg", height=4, width=1.5)

###ScerXSkud###
ScerXSkudList <- read.table("ScerXSkud_wPAD1FDC1.txt", header=T)
wStrains <- unique(ScerXSkudList$genomeName)
PFwant <- PFends[which(PFends$chrName=="Scer_chrIV" | PFends$chrName=="SkudZP_chrIV"),]

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in wStrains) {
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
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}
toPlotTemp <- plotTable[which(plotTable$chrName=="Scer_chrIV" | plotTable$chrName=="SkudZP_chrIV"),]
toPlotTemp$newStart <- toPlotTemp$chrStart
toPlotTemp$newStart[which(toPlotTemp$species=="Scer" & toPlotTemp$chrStart<ScerHalf & toPlotTemp$chrEnd>ScerHalf)] <- ScerHalf+1
toPlotTemp$newStart[which(toPlotTemp$species=="SkudZP" & toPlotTemp$chrStart<SkudHalf & toPlotTemp$chrEnd>SkudHalf)] <- SkudHalf+1
toPlot <- rbind(toPlotTemp[which(toPlotTemp$species=="Scer" & toPlotTemp$newStart>ScerHalf),], toPlotTemp[which(toPlotTemp$species=="SkudZP" & toPlotTemp$newStart>SkudHalf),])
#toPlot$chrTemp <- as.character(toPlot$chr)
#toPlot$chr[which(toPlot$chrTemp=="SkudZP_chrIV")] <- "Skud_chrIV"

yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)

genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Scer" | PAD1FDC1geneLoc$species=="Skud"),]
genePos$chr <- genePos$chrom
PFwant$chr <- factor(as.character(PFwant$chrName), levels=c("Scer_chrIV", "SkudZP_chrIV"))
toPlot$chr <- factor(as.character(toPlot$chrName), levels=c("Scer_chrIV", "SkudZP_chrIV"))

PFwant$chr <- as.character(PFwant$chrName)
PFwant$chr[which(PFwant$chr=="SkudZP_chrIV")] <- "Skud_chrIV"

ggplot()+geom_rect(data=toPlot, aes(xmin = newStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlot$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends), linetype="dotted")+geom_vline(data=genePos, aes(xintercept = start)) +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())+geom_vline(data=PFwant, aes(xintercept = half), linetype="dashed", color="grey")

toPlotMod <- toPlot
toPlotMod$newStart[which(toPlotMod$strainID=="FM1043_ScerXSkud" & toPlotMod$newStart==1520000)] <- 1480000

ggplot()+geom_rect(data=toPlotMod, aes(xmin = newStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlotMod$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends))+geom_vline(data=genePos, aes(xintercept = start), linetype="dotted") +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())+geom_vline(data=PFwant, aes(xintercept = half), linetype="dashed", color="grey")

ggsave("ScerXSkud_endChrPAD1FDC1.jpg", height=4, width=1.5)



##lager
ScerXSeubList <- read.table("ScerXSeub_wPAD1FDC1.txt", header=T)
lagStrains <- unique(ScerXSeubList$genomeName)
PFwant <- PFends[which(PFends$chrName=="Scer_chrIV" | PFends$chrName=="Seub_chrXIII"),]

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in lagStrains) {
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
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}
toPlotTemp <- plotTable[which(plotTable$chrName=="Scer_chrIV" | plotTable$chrName=="Seub_chrXIII"),]
toPlotTemp$newStart <- toPlotTemp$chrStart
toPlotTemp$newStart[which(toPlotTemp$species=="Scer" & toPlotTemp$chrStart<ScerHalf & toPlotTemp$chrEnd>ScerHalf)] <- ScerHalf+1
toPlotTemp$newStart[which(toPlotTemp$species=="Seub" & toPlotTemp$chrStart<SeubHalf & toPlotTemp$chrEnd>SeubHalf)] <- SeubHalf+1
toPlot <- rbind(toPlotTemp[which(toPlotTemp$species=="Scer" & toPlotTemp$newStart>ScerHalf),], toPlotTemp[which(toPlotTemp$species=="Seub" & toPlotTemp$newStart>SeubHalf),])

yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
#colorScale <- scale_color_manual(name="Species", values = speciesColors, breaks=uniSpecies)

genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Scer" | PAD1FDC1geneLoc$species=="Seub"),]
genePos$chr <- factor(as.character(genePos$chrom), levels=c("Scer_chrIV", "Seub_chrXIII"))
PFwant$chr <- factor(as.character(PFwant$chrName), levels=c("Scer_chrIV", "Seub_chrXIII"))
toPlot$chr <- factor(as.character(toPlot$chrName), levels=c("Scer_chrIV", "Seub_chrXIII"))

ggplot()+geom_rect(data=toPlot, aes(xmin = newStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlot$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends), linetype="dotted")+geom_vline(data=genePos, aes(xintercept = start)) +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())+geom_vline(data=PFwant, aes(xintercept = half), linetype="dashed", color="grey")

toPlotMod <- toPlot
toPlotMod <- toPlotMod[-which(toPlotMod$species=="Seub" & toPlotMod$newStart>=910000),]
toPlotMod$chrEnd[which(toPlotMod$chrEnd==1520000)] <- 1530000
toPlotMod$chrEnd[which(toPlotMod$chrEnd==870000)] <- 880000
toPlotMod$newStart[which(toPlotMod$strainID=="yHAB585_lager" & toPlotMod$newStart==1460000)] <- 1440001
toPlotMod$newStart[which(toPlotMod$strainID=="DBVPG6261_lager_Hewitt" & toPlotMod$newStart==870000)] <- 860001
toPlotMod <- rbind(toPlotMod, data.frame(strainID="yHAB589_lager", chrName="Scer_chrIV", ploidy=1, genomeStart=2890000, genomeEnd=2900000, chrStart=1510000, chrEnd=1520000, regionLength=10000, species="Scer", chrNum="chrIV", spStart=2890000, spEnd=2900000, newGenoStart=2890000, newGenoEnd=2900000, yStart=13, yEnd=14, alpha=0.76, newStart=1510000, chr="Scer_chrIV"))

ggplot()+geom_rect(data=toPlotMod, aes(xmin = newStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlotMod$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends))+geom_vline(data=genePos, aes(xintercept = start), linetype="dotted") +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())+geom_vline(data=PFwant, aes(xintercept = half), linetype="dashed", color="grey")


ggsave("ScerXSeub_endChrPAD1FDC1.jpg", height=8, width=1.5)


####Last 20Kb ###
PFends$last50 <- PFends$ends-50000
Scerlast50 <- PFends$last50[which(PFends$chrName=="Scer_chrIV")]
Skudlast50 <- PFends$last50[which(PFends$chrName=="SkudZP_chrIV")]
Suvalast50 <- PFends$last50[which(PFends$chrName=="Suva_chrII")]
Seublast50 <- PFends$last50[which(PFends$chrName=="Seub_chrXIII")]

##Subset#
SuvaXSeubList <- read.table("SuvaXSeub_wPAD1FDC1.txt", header=T)
ueStrains <- unique(SuvaXSeubList$genomeName)
PFwant <- PFends[which(PFends$chrName=="Suva_chrII" | PFends$chrName=="Seub_chrXIII"),]

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in ueStrains) {
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
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}
toPlotTemp <- plotTable[which(plotTable$chrName=="Suva_chrII" | plotTable$chrName=="Seub_chrXIII"),]
toPlotTemp$newStart <- toPlotTemp$chrStart
toPlotTemp$newStart[which(toPlotTemp$species=="Suva" & toPlotTemp$chrStart<Suvalast50 & toPlotTemp$chrEnd>Suvalast50)] <- Suvalast50+1
toPlotTemp$newStart[which(toPlotTemp$species=="Seub" & toPlotTemp$chrStart<Seublast50 & toPlotTemp$chrEnd>Seublast50)] <- Seublast50+1
toPlot <- rbind(toPlotTemp[which(toPlotTemp$species=="Suva" & toPlotTemp$newStart>Suvalast50),], toPlotTemp[which(toPlotTemp$species=="Seub" & toPlotTemp$newStart>Seublast50),])

yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
#colorScale <- scale_color_manual(name="Species", values = speciesColors, breaks=uniSpecies)

genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Suva" | PAD1FDC1geneLoc$species=="Seub"),]
genePos$chr <- factor(as.character(genePos$chrom), levels=c("Suva_chrII", "Seub_chrXIII"))
PFwant$chr <- factor(as.character(PFwant$chrName), levels=c("Suva_chrII", "Seub_chrXIII"))
toPlot$chr <- factor(as.character(toPlot$chrName), levels=c("Suva_chrII", "Seub_chrXIII"))

ggplot()+geom_rect(data=toPlot, aes(xmin = chrStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlot$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends), linetype="dotted")+geom_vline(data=genePos, aes(xintercept = start)) +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())

ggsave("SuvaXSeub_last50PAD1FDC1.jpg", height=4, width=1.5)

###ScerXSkud###
ScerXSkudList <- read.table("ScerXSkud_wPAD1FDC1.txt", header=T)
wStrains <- unique(ScerXSkudList$genomeName)
PFwant <- PFends[which(PFends$chrName=="Scer_chrIV" | PFends$chrName=="SkudZP_chrIV"),]

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in wStrains) {
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
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}
toPlotTemp <- plotTable[which(plotTable$chrName=="Scer_chrIV" | plotTable$chrName=="SkudZP_chrIV"),]
toPlotTemp$newStart <- toPlotTemp$chrStart
toPlotTemp$newStart[which(toPlotTemp$species=="Scer" & toPlotTemp$chrStart<Scerlast50 & toPlotTemp$chrEnd>Scerlast50)] <- Scerlast50+1
toPlotTemp$newStart[which(toPlotTemp$species=="SkudZP" & toPlotTemp$chrStart<Skudlast50 & toPlotTemp$chrEnd>Skudlast50)] <- Skudlast50+1
toPlot <- rbind(toPlotTemp[which(toPlotTemp$species=="Scer" & toPlotTemp$newStart>Scerlast50),], toPlotTemp[which(toPlotTemp$species=="SkudZP" & toPlotTemp$newStart>Skudlast50),])
#toPlot$chrTemp <- as.character(toPlot$chr)
#toPlot$chr[which(toPlot$chrTemp=="SkudZP_chrIV")] <- "Skud_chrIV"

yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)

genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Scer" | PAD1FDC1geneLoc$species=="Skud"),]
genePos$chr <- genePos$chrom
PFwant$chr <- factor(as.character(PFwant$chrName), levels=c("Scer_chrIV", "SkudZP_chrIV"))
toPlot$chr <- factor(as.character(toPlot$chrName), levels=c("Scer_chrIV", "SkudZP_chrIV"))

PFwant$chr <- as.character(PFwant$chrName)
PFwant$chr[which(PFwant$chr=="SkudZP_chrIV")] <- "Skud_chrIV"

ggplot()+geom_rect(data=toPlot, aes(xmin = chrStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlot$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends), linetype="dotted")+geom_vline(data=genePos, aes(xintercept = start)) +facet_grid(~chr, scales = "free")+theme(axis.text.x = element_blank())

ggsave("ScerXSkud_last50PAD1FDC1.jpg", height=4, width=1.5)



##lager
ScerXSeubList <- read.table("ScerXSeub_wPAD1FDC1.txt", header=T)
lagStrains <- unique(ScerXSeubList$genomeName)
PFwant <- PFends[which(PFends$chrName=="Scer_chrIV" | PFends$chrName=="Seub_chrXIII"),]

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in lagStrains) {
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
  strainYstart <- spMaxY
  plotTable <- rbind(plotTable, strainData)
}
toPlotTemp <- plotTable[which(plotTable$chrName=="Scer_chrIV" | plotTable$chrName=="Seub_chrXIII"),]
toPlotTemp$newStart <- toPlotTemp$chrStart
toPlotTemp$newStart[which(toPlotTemp$species=="Scer" & toPlotTemp$chrStart<Scerlast50 & toPlotTemp$chrEnd>Scerlast50)] <- Scerlast50+1
toPlotTemp$newStart[which(toPlotTemp$species=="Seub" & toPlotTemp$chrStart<Seublast50 & toPlotTemp$chrEnd>Seublast50)] <- Seublast50+1
toPlot <- rbind(toPlotTemp[which(toPlotTemp$species=="Scer" & toPlotTemp$newStart>Scerlast50),], toPlotTemp[which(toPlotTemp$species=="Seub" & toPlotTemp$newStart>Seublast50),])

yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
#colorScale <- scale_color_manual(name="Species", values = speciesColors, breaks=uniSpecies)

genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Scer" | PAD1FDC1geneLoc$species=="Seub"),]
genePos$chr <- factor(as.character(genePos$chrom), levels=c("Scer_chrIV", "Seub_chrXIII"))
PFwant$chr <- factor(as.character(PFwant$chrName), levels=c("Scer_chrIV", "Seub_chrXIII"))
toPlot$chr <- factor(as.character(toPlot$chrName), levels=c("Scer_chrIV", "Seub_chrXIII"))

ggplot()+geom_rect(data=toPlot, aes(xmin = newStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=toPlot$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(expand = c(0,0))+geom_vline(data=PFwant, aes(xintercept = ends), linetype="dotted")+geom_vline(data=genePos, aes(xintercept = start))+geom_vline(data=PFwant, aes(xintercept = last50), linetype="dashed", color="grey") +facet_grid(~chr, scales = "free")

+theme(axis.text.x = element_blank())

ggsave("ScerXSeub_last50PAD1FDC1.jpg", height=4, width=1.5)




ScerXSkudList <- read.table("ScerXSkud_wPAD1FDC1.txt", header=T)
ckStrains <- unique(ScerXSkudList$genomeName)

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in ckStrains) {
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


yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)
#colorScale <- scale_color_manual(name="Species", valcks = speciesColors, breaks=uniSpecies)

#plotTitle <- ggtitle("ScerXSkud with PAD1 FDC1")
#ggplot()+geom_rect(data=plotTable, aes(xmin = genomeStart, xmax = genomeEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=plotTable$alpha)+fillScale+xaxis+yaxis+chrLines+strainLines+theme_classic()+plotTitle+theme(axis.text.y = element_blank())

chrsWanted <- plotTable[which(plotTable$chrName=="Scer_chrIV" | plotTable$chrName=="SkudZP_chrIV"),]

ggplot()+geom_rect(data=chrsWanted, aes(xmin = chrStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=chrsWanted$alpha)+
  fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + 
  scale_x_continuous(limits = c(0, NA))+facet_grid(~chrName, scales = "free")

#PAD1FDC1geneLoc <- read.table("../genesOfInterest/PAD1-FDC1_geneLoc.txt", header=T)
#PAD1FDC1geneLoc$chrName <- PAD1FDC1geneLoc$chrom
#PFchrLens <- read.table("../genesOfInterest/PAD1-FDC1_chrLens.txt", header=T)
genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Scer" | PAD1FDC1geneLoc$species=="Skud"),]
chrPos <- PFchrLens[which(PFchrLens$chrName=="Scer_chrIV" | PFchrLens$chrName=="SkudZP_chrIV"),]
#chrsWanted$chr <- factor(as.character(chrsWanted$chrName), levels=c("Scer_chrIV", "Skud_chrIV"))
#chrWanted$chrName <-factor(chrWanted$chr, levels=c("Scer_chrIV", "Skud_chrIV"))
#genePos$chr <- factor(as.character(genePos$chrom), levels=c("Scer_chrIV", "Skud_chrIV"))
#chrPos$chr <- factor(as.character(chrPos$chrName), levels=c("Scer_chrIV", "Skud_chrIV"))
genePos$chrName <- as.character(genePos$chrom)
genePos$chrName[which(genePos$chrom=="Skud_chrIV")] <- "SkudZP_chrIV"

ggplot()+geom_rect(data=chrsWanted, aes(xmin = chrStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=chrsWanted$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(limits = c(0, NA), name="Chr Pos")+geom_vline(data=chrPos, aes(xintercept = ends), linetype="dotted")+geom_vline(data=genePos, aes(xintercept = start)) +facet_grid(~chrName, scales = "free")+theme(axis.text.x = element_blank())

ggsave("ScerXSkud_chrsPAD1FDC1.jpg", height=4, width=1.5)


###ScerXSeub###
ScerXSeubList <- read.table("ScerXSeub_wPAD1FDC1.txt", header=T)
ceStrains <- unique(ScerXSeubList$genomeName)

plotTable <- data.frame()
strainYstart <- 0
yLabels <- c()
yBreaks <- 0
for (strainName in ceStrains) {
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

yaxis <- scale_y_continuous(limits = c(0,max(plotTable$yEnd)), name="Strain", breaks=yLabels, expand = c(0, 0))
strainLines <- geom_hline(yintercept = yBreaks, alpha=0.4)

fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies, guide=FALSE)

chrsWanted <- plotTable[which(plotTable$chrName=="Scer_chrIV" | plotTable$chrName=="Seub_chrXIII"),]

ggplot()+geom_rect(data=chrsWanted, aes(xmin = chrStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=chrsWanted$alpha)+
  fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + 
  scale_x_continuous(limits = c(0, NA))+facet_grid(~chrName, scales = "free")

genePos <- PAD1FDC1geneLoc[which(PAD1FDC1geneLoc$species=="Scer" | PAD1FDC1geneLoc$species=="Seub"),]
chrPos <- PFchrLens[which(PFchrLens$chrName=="Scer_chrIV" | PFchrLens$chrName=="Seub_chrXIII"),]
genePos$chrName <- as.character(genePos$chrom)

ggplot()+geom_rect(data=chrsWanted, aes(xmin = chrStart, xmax = chrEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=chrsWanted$alpha)+fillScale+yaxis+strainLines+theme_classic()+theme(axis.text.y = element_blank()) + scale_x_continuous(limits = c(0, NA), name="Chr Pos")+geom_vline(data=chrPos, aes(xintercept = ends), linetype="dotted")+geom_vline(data=genePos, aes(xintercept = start)) +facet_grid(~chrName, scales = "free")+theme(axis.text.x = element_blank())

ggsave("ScerXSeub_chrsPAD1FDC1.jpg", height=8, width=1.5)


#dev.off()

####ChrXIII####
SeubChrXIII <- plotTable[which(plotTable$chrName=="Seub_chrXIII"),]
strains <- unique(plotTable$strainID)
lastXIII <- data.frame()
for (strainName in strains) {
  strainXIII <- SeubChrXIII[which(SeubChrXIII$strainID==strainName),]
  lastWins <- strainXIII[(nrow(strainXIII)-1):nrow(strainXIII),]
  lastXIII <- rbind(lastXIII, lastWins)
} 


ScerChrIV <- plotTable[which(plotTable$chrName=="Scer_chrIV"),]
strains <- unique(plotTable$strainID)
lastIV <- data.frame()
for (strainName in strains) {
  strainIV <- ScerChrIV[which(ScerChrIV$strainID==strainName),]
  lastWins <- strainIV[(nrow(strainIV)-1):nrow(strainIV),]
  lastIV <- rbind(lastIV, lastWins)
} 

SkudChrIV <- plotTable[which(plotTable$chrName=="SkudZP_chrIV"),]
strains <- unique(plotTable$strainID)
SkudLIV <- data.frame()
for (strainName in strains) {
  strainIV <- SkudChrIV[which(SkudChrIV$strainID==strainName),]
  lastWins <- strainIV[(nrow(strainIV)-1):nrow(strainIV),]
  SkudLIV <- rbind(SkudLIV, lastWins)
} 

SuvaChrII <- plotTable[which(plotTable$chrName=="Suva_chrII"),]
strains <- unique(plotTable$strainID)
lastII <- data.frame()
for (strainName in strains) {
  strainII <- SuvaChrII[which(SuvaChrII$strainID==strainName),]
  lastWins <- strainII[(nrow(strainII)-1):nrow(strainII),]
  lastII <- rbind(lastII, lastWins)
} 


finalSet <- data.frame()



WLP802 <- ploidyTable[which(ploidyTable$strainID=="yHCT131_lager"),]
chrXIII <- WLP802[which(WLP802$chrNum=="chrXIII"),]
ggplot(chrXIII)+geom_rect(aes(xmin=chrStart, xmax=chrEnd, ymin=0, ymax=ploidy, fill=species), alpha=0.7)+fillScale+theme_classic()+facet_grid(species~.)

+geom_rect(data=plotTable, aes(xmin = newGenoStart, xmax = newGenoEnd, ymin =yStart, ymax=yEnd, fill=species), alpha=plotTable$alpha)+fillScale+xaxis+yaxis+chrLines+strainLines+theme_classic()+plotTitle+theme(axis.text.y = element_blank())

chrSplit <- strsplit(as.character(chrLengths$chrName), "_")
chrLengths$species <- unlist(lapply(chrSplit, function(x) x[1]))


