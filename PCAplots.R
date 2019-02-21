options(scipen=999)
options(stringsAsFactors = F)
library(ggplot2)
setwd("~/Documents/Hittinger_Lab/whole_genome/brewingProject/PCAs")

colors <- c("Domesticated"="#d52900", "Hybrid"="#cc79a7", "Wild"="#0072b2")
fillScale <- scale_fill_manual(values=colors, name="Origin")

###Used
##Skud
popColors <- c("AsiaB"="#404040", "AsiaA"="#848484", "Europe"="#ACACAC", "EU"="#ACACAC", "Hybrids"="#CBCBCB", "Wine Hybrids"="#CBCBCB", "BeerHybrids"="#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))

SkudSum <- read.table("Skud_PCAsummary.txt", check.names = F, header=T)
SkudEUasiaA <- read.table("PCAscores_Skud_EUasiaA.txt", check.names = F)
SkudStrainsEUasiaA <-row.names(SkudEUasiaA)
split <- strsplit(SkudStrainsEUasiaA, "_")
SkudEUasiaA$strainID <- unlist(lapply(split, function(x) x[1]))
SkudEUasiaA$syn <- unlist(lapply(split, function(x) x[2]))
SkudEUasiaA$spIso <- unlist(lapply(split, function(x) x[3]))
SkudEUasiaA$pop <- unlist(lapply(split, function(x) x[4]))
SkudEUasiaA$spIso[which(SkudEUasiaA$strainID=="yHKS141")] <- "Skud"
SkudEUasiaA$pop[which(SkudEUasiaA$strainID=="yHKS141")] <- "EU"
SkudEUasiaA$syn[which(SkudEUasiaA$strainID=="yHKS141")] <- NA
SkudEUasiaA$origin <- "Hybrid"
SkudEUasiaA$origin[which(SkudEUasiaA$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")

SkudEUasiaA$popLab <- SkudEUasiaA$pop
SkudEUasiaA$popLab[which(SkudEUasiaA$pop=="Europe")] <- "Europe"
SkudEUasiaA$popLab[which(SkudEUasiaA$pop=="EU")] <- "Europe"
SkudEUasiaA$popLab[which(SkudEUasiaA$spIso=="CIDER")] <- "Hybrids"
SkudEUasiaA$popLab[which(SkudEUasiaA$spIso=="WINE")] <- "Hybrids"
SkudEUasiaA$popLab[which(SkudEUasiaA$spIso=="BEER")] <- "Hybrids"
uniPops <- unique(SkudEUasiaA$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SkudEUasiaA[which(SkudEUasiaA$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
#popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSkud"| popLabDF$pop=="lager-saaz"),]
mid <- min(SkudEUasiaA$PC1)+((max(SkudEUasiaA$PC1)-min(SkudEUasiaA$PC1))/2)
spLabDF <- data.frame(x=mid, y=min(SkudEUasiaA$PC2), lab="S. kudriavzevii")

####This is the best###
ggplot(SkudEUasiaA, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=1,height=5), hjust = "inward") + 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), fontface="italic", vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=4),axis.title=element_text(size=6),text = element_text(size =4))

ggsave("SkudEUasiaA_PC1-PC2_25-35.jpg", height=2.5, width=3.5)


###### Seub #########
SeubSum <- read.table("Seub_PCAsummary.txt", check.names = F, header=T)
popColors <- c("PA1"="#404040", "PA2"="#6D6D6D", "PB2"="#8A8A8A", "PB1"="#A2A2A2",  "Holarctic"="#B6B6B6",  "HOL"="#B6B6B6", "Complex Hybrids"="#C7C7C7", "complexHyb"="#C7C7C7", "SuvaXSeub"="#C7C7C7", "Lager-Saaz"="#D7D7D7", "Lager-Frohberg" = "#E6E6E6", "Lager" = "#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))
### All ###
SeubAll <- read.table("PCAscores_Seub_All.txt", check.names = F, sep="\t")
SeubStrainsAll <-row.names(SeubAll)
split <- strsplit(SeubStrainsAll, "_")
SeubAll$strainID <- unlist(lapply(split, function(x) x[1]))
SeubAll$syn <- unlist(lapply(split, function(x) x[2]))
SeubAll$spIso <- unlist(lapply(split, function(x) x[3]))
SeubAll$pop <- unlist(lapply(split, function(x) x[4]))
SeubAll$extra <- unlist(lapply(split, function(x) x[5]))
SeubAll$extra[which(SeubAll$extra=="1")] <- NA
SeubAll$extra[which(SeubAll$extra=="2")] <- NA
SeubAll$spIso[which(is.na(SeubAll$extra)==F)] <- SeubAll$pop[which(is.na(SeubAll$extra)==F)]
SeubAll$pop[which(is.na(SeubAll$extra)==F)] <- SeubAll$extra[which(is.na(SeubAll$extra)==F)]
SeubAll$origin <- "Hybrid"
SeubAll$origin[which(SeubAll$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")

SeubAll$popLab <- SeubAll$pop
SeubAll$popLab[which(SeubAll$pop=="lager-")] <- "Lager-Frohberg"
SeubAll$popLab[which(SeubAll$pop=="SuvaXSeub")] <- "complexHyb"
SeubAll$popLab[which(SeubAll$pop=="SuvaXSkudXSeub")] <- "complexHyb"
uniPops <- unique(SeubAll$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SeubAll[which(SeubAll$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSeub"| popLabDF$pop=="lager-saaz"),]
mid <- min(SeubAll$PC1)+((max(SeubAll$PC1)-min(SeubAll$PC1))/2)
spLabDF <- data.frame(x=mid, y=min(SeubAll$PC2), lab="S. eubayanus")
popLabTrim$pop[which(popLabTrim$pop=="Lager-Frohberg")] <- "Lager"
popLabTrim$pop[which(popLabTrim$pop=="complexHyb")] <- "Complex Hybrids"
popLabTrim$pop[which(popLabTrim$pop=="PB1")] <- "PopB-1"
popLabTrim$pop[which(popLabTrim$pop=="PB2")] <- "PopB-2"
popLabTrim$pop[which(popLabTrim$pop=="PA1")] <- "PopA-1"
popLabTrim$pop[which(popLabTrim$pop=="PA2")] <- "PopA-2"
popLabTrim$pop[which(popLabTrim$pop=="HOL")] <- "Holarctic"
popLabTrim <- popLabTrim[-which(popLabTrim$pop=="PopB-2" | popLabTrim$pop=="PopA-2" | popLabTrim$pop=="Complex Hybrids"),]
popLabTrim$pop[which(popLabTrim$pop=="Lager")] <- "Hybrids"
popLabTrim$pop[which(popLabTrim$pop=="PopA-1")] <- "PopA"
popLabTrim$pop[which(popLabTrim$pop=="PopB-1")] <- "PopB"

####This is the best###
ggplot(SeubAll, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabTrim, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=5,height=-15), hjust = "inward")+ 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), fontface="italic", vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=4),axis.title=element_text(size=6),text = element_text(size =4))

ggsave("SeubAll_PC1-PC2.jpg", height=2.5, width=3.5)


##Suva
SuvaSum <- read.table("Suva_PCAsummary.txt", check.names = F, header=T)
popColors <- c("Aust"="#404040", "SA-B"="#6D6D6D", "SA-A"="#C1C1C1", "Holarctic"="#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))
SuvaSoAmHol <- read.table("PCAscores_Suva_SoAmHolarc.txt", check.names = F, sep="\t")
SuvaStrainsSoAmHol <-row.names(SuvaSoAmHol)
split <- strsplit(SuvaStrainsSoAmHol, "_")
SuvaSoAmHol$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaSoAmHol$syn <- unlist(lapply(split, function(x) x[2]))
SuvaSoAmHol$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaSoAmHol$pop <- unlist(lapply(split, function(x) x[4]))
SuvaSoAmHol$origin <- "Hybrid"
SuvaSoAmHol$origin[which(SuvaSoAmHol$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")

SuvaSoAmHol$popLab <- SuvaSoAmHol$pop
SuvaSoAmHol$popLab[-which(SuvaSoAmHol$pop=="SA-A" | SuvaSoAmHol$pop=="SA-B")] <- "Holarctic"
uniPops <- unique(SuvaSoAmHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SuvaSoAmHol[which(SuvaSoAmHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
mid <- min(SuvaSoAmHol$PC1)+((max(SuvaSoAmHol$PC1)-min(SuvaSoAmHol$PC1))/2)
spLabDF <- data.frame(x=mid, y=min(SuvaSoAmHol$PC2), lab="S. uvarum")
popLabDF$pop[which(popLabDF$pop=="Holarctic")] <- "Holarctic\n& Hyrbids"

####This is the best###
ggplot(SuvaSoAmHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=5,height=-15), hjust = "inward")+ 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), fontface="italic", vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=4),axis.title=element_text(size=6),text = element_text(size =4))

ggsave("SuvaSoAmHol_PC1-PC2.jpg", height=2.5, width=3.5)



###Scer###
renameOri <- read.table("strainPopOriKey.txt", check.names = F, header = T)
ScerSum <- read.table("Scer_PCAsummary.txt", check.names = F, header=T)
popColors <- c("WildMisc"="#404040", "SakeAsia"="#6D6D6D", "Mosaic"="#8A8A8A", "BreadMixed"="#A2A2A2",  "Beer2"="#B6B6B6", "MedOak"="#C7C7C7", "Wine"="#D7D7D7", "AleBeer1" = "#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))
### All ###
Scer.galGonAll <- read.table("PCAscores_Scer_galGonAll.txt", check.names = F, sep="\t")
Scer.galGonAll.strains <-row.names(Scer.galGonAll)
split <- strsplit(Scer.galGonAll.strains, "_")
Scer.galGonAll$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.galGonAll$syn <- unlist(lapply(split, function(x) x[2]))
Scer.galGonAll$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.galGonAll$tempPop <- unlist(lapply(split, function(x) x[4]))
Scer.galGonAll$phase <- unlist(lapply(split, function(x) x[5]))
Scer.galGonAll$phase[which(Scer.galGonAll$pop==1)] <- 1
Scer.galGonAll$phase[which(Scer.galGonAll$pop==2)] <- 2
Scer.galGonAll$tempPop[which(Scer.galGonAll$tempPop==1)] <- "Unplaced"
Scer.galGonAll$tempPop[which(Scer.galGonAll$tempPop==2)] <- "Unplaced"
Scer.galGonAll$origin <- NA
Scer.galGonAll$pop <- NA
for (strainID in Scer.galGonAll.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  strainPop <- as.character(renameOri$pop[which(renameOri$rename==strainID)])
  Scer.galGonAll$origin[which(rownames(Scer.galGonAll)==strainID)] <- strainOri
  Scer.galGonAll$pop[which(rownames(Scer.galGonAll)==strainID)] <- strainPop
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")

uniPops <- unique(Scer.galGonAll$pop)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- Scer.galGonAll[which(Scer.galGonAll$pop==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF <- rbind(popLabDF, labTemp)
}
mid <- min(Scer.galGonAll$PC1)+((max(Scer.galGonAll$PC1)-min(Scer.galGonAll$PC1))/2)
spLabDF <- data.frame(x=mid, y=min(Scer.galGonAll$PC2), lab="S. cerevisiae")
popLabMod <- popLabDF
popLabMod$pop[which(popLabMod$pop=="AleBeer1")] <-"Ale/Beer1" 
popLabMod$pop[which(popLabMod$pop=="BreadMixed")] <-"Bread/Mixed" 
popLabMod$pop[which(popLabMod$pop=="WildMisc")] <-"Wild Misc" 
popLabMod$pop[which(popLabMod$pop=="SakeAsia")] <-"Sake/Asian" 
popLabMod <- popLabMod[-which(popLabMod$pop=="MedOak"),]

####This is the best###
ggplot(Scer.galGonAll, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=pop), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabMod, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=5,height=-15), hjust = "inward")+ 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), fontface="italic", vjust = "inward", size = 4) +
  theme(axis.text=element_text(size=6),axis.title=element_text(size=8),text = element_text(size =6), legend.text=element_text(size=7), legend.title=element_text(size=8))

ggsave("ScerAll_PC1-PC2.jpg", height=3, width=5.5)

###For Supplment
#### Beer1_lessFro ###
renameOri <- read.table("beer1strainPopKey_wheatStout.txt", check.names = F, header = T)
popColors <- c("aleMosaic"="#666666", "belGerAltKol"="#909090", "beerBelGerAltKol"="#909090", "US"="#AFAFAF", "beerUS"="#AFAFAF", "britIsles"="#C9C9C9", "beerBritishIsles"="#C9C9C9",  "lager-saaz" = "#DFDFDF",  "Lager-Saaz" = "#DFDFDF", "lager-frohberg"="#F2F2F2", "Lager-Frohberg"="#F2F2F2", "Wheat"="#666666", "Stout"="#666666")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))

Scer.Beer1_lessFro <- read.table("PCAscores_Scer_Beer1lF.txt", check.names = F, sep="\t")
Scer.Beer1_lessFro.strains <-row.names(Scer.Beer1_lessFro)
split <- strsplit(Scer.Beer1_lessFro.strains, "_")
Scer.Beer1_lessFro$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Beer1_lessFro$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Beer1_lessFro$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Beer1_lessFro$tempPop <- unlist(lapply(split, function(x) x[4]))
Scer.Beer1_lessFro$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Beer1_lessFro$phase[which(Scer.Beer1_lessFro$pop==1)] <- 1
Scer.Beer1_lessFro$phase[which(Scer.Beer1_lessFro$pop==2)] <- 2
Scer.Beer1_lessFro$tempPop[which(Scer.Beer1_lessFro$tempPop==1)] <- "Unplaced"
Scer.Beer1_lessFro$tempPop[which(Scer.Beer1_lessFro$tempPop==2)] <- "Unplaced"
Scer.Beer1_lessFro$origin <- NA
Scer.Beer1_lessFro$pop <- NA
for (strainID in Scer.Beer1_lessFro.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  strainPop <- as.character(renameOri$pop[which(renameOri$rename==strainID)])
  Scer.Beer1_lessFro$origin[which(rownames(Scer.Beer1_lessFro)==strainID)] <- strainOri
  Scer.Beer1_lessFro$pop[which(rownames(Scer.Beer1_lessFro)==strainID)] <- strainPop
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")

uniPops <- unique(Scer.Beer1_lessFro$pop)
uniPops <- uniPops[-which(uniPops=="aleMosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- Scer.Beer1_lessFro[which(Scer.Beer1_lessFro$pop==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF <- rbind(popLabDF, labTemp)
}
mid <- min(Scer.Beer1_lessFro$PC1)+((max(Scer.Beer1_lessFro$PC1)-min(Scer.Beer1_lessFro$PC1))/2)
spLabDF <- data.frame(x=mid, y=max(Scer.Beer1_lessFro$PC2), lab="Ale/Beer1 - less Frohberg")
popLabMod <- popLabDF
popLabMod$pop[which(popLabMod$pop=="lager-saaz")] <- "Lager-Saaz"
popLabMod$pop[which(popLabMod$pop=="lager-frohberg")] <- "Lager-Frohberg"
popLabMod$pop[which(popLabMod$pop=="belGerAltKol")] <- "Belgian/German/Alt-Kolsch"
popLabMod$pop[which(popLabMod$pop=="britIsles")] <- "British Isles"

####This is the best###
ggplot(Scer.Beer1_lessFro, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=pop), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabMod, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=5,height=5), hjust = "inward")+ 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), vjust = "inward", size = 4) +
  theme(axis.text=element_text(size=6),axis.title=element_text(size=8),text = element_text(size =6), legend.text=element_text(size=7), legend.title=element_text(size=8))+guides(fill=F)

ggsave("ScerBeer1lessFro_wheatStout_PC1-PC2.jpg", height=4, width=6)

#AllFrohberg
Scer.Beer1 <- read.table("PCAscores_Scer_Beer1.txt", check.names = F, sep="\t")
Scer.Beer1.strains <-row.names(Scer.Beer1)
split <- strsplit(Scer.Beer1.strains, "_")
Scer.Beer1$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Beer1$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Beer1$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Beer1$pop <- unlist(lapply(split, function(x) x[4]))
Scer.Beer1$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Beer1$phase[which(Scer.Beer1$pop==1)] <- 1
Scer.Beer1$phase[which(Scer.Beer1$pop==2)] <- 2
Scer.Beer1$pop[which(Scer.Beer1$pop==1)] <- "Unplaced"
Scer.Beer1$pop[which(Scer.Beer1$pop==2)] <- "Unplaced"
Scer.Beer1$origin <- NA
for (strainID in Scer.Beer1.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.Beer1$origin[which(rownames(Scer.Beer1)==strainID)] <- strainOri
  strainPop <- as.character(renameOri$pop[which(renameOri$rename==strainID)])
  Scer.Beer1$pop[which(rownames(Scer.Beer1)==strainID)] <- strainPop
  
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")

uniPops <- unique(Scer.Beer1$pop)
uniPops <- uniPops[-which(uniPops=="aleMosaic")]
#uniPops <- uniPops[-which(uniPops=="Unplaced")]
#uniPops <- uniPops[-which(uniPops=="beer")]
#uniPops <- uniPops[-which(uniPops=="beerWheat")]
#uniPops <- uniPops[-which(is.na(uniPops))]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- Scer.Beer1[which(Scer.Beer1$pop==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF <- rbind(popLabDF, labTemp)
}
mid <- min(Scer.Beer1$PC1)+((max(Scer.Beer1$PC1)-min(Scer.Beer1$PC1))/2)
spLabDF <- data.frame(x=mid, y=max(Scer.Beer1$PC2), lab="Ale/Beer1 - all Frohberg")
popLabMod <- popLabDF
popLabMod$pop[which(popLabMod$pop=="lager-saaz")] <- "Lager-Saaz"
popLabMod$pop[which(popLabMod$pop=="lager-frohberg")] <- "Lager-Frohberg"
popLabMod$pop[which(popLabMod$pop=="belGerAltKol")] <- "Belgian/German/Alt-Kolsch"
popLabMod$pop[which(popLabMod$pop=="britIsles")] <- "British Isles"
popLabMod$pop[which(popLabMod$pop=="beerUS")] <- "US"

####This is the best###
ggplot(Scer.Beer1, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=pop), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabMod, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=1,height=1), hjust = "inward")+ 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), vjust = "inward", size = 4) +
  theme(axis.text=element_text(size=6),axis.title=element_text(size=8),text = element_text(size =6), legend.text=element_text(size=7), legend.title=element_text(size=8))+guides(fill=F)

ggsave("ScerBeer1allFro_WheatStout_PC1-PC2.jpg", height=4, width=6)



##Skud Supp
popColors <- c("AsiaB"="#404040", "AsiaA"="#848484", "Europe"="#ACACAC", "EU"="#ACACAC", "Hybrids"="#CBCBCB", "Wine Hybrids"="#CBCBCB", "WineHybs"="#CBCBCB", "BeerHybrids"="#E6E6E6", "BeerHybs"="#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))

SkudSum <- read.table("Skud_PCAsummary.txt", check.names = F, header=T)
SkudEU <- read.table("PCAscores_Skud_EU.txt", check.names = F)
SkudStrainsEU <-row.names(SkudEU)
split <- strsplit(SkudStrainsEU, "_")
SkudEU$strainID <- unlist(lapply(split, function(x) x[1]))
SkudEU$syn <- unlist(lapply(split, function(x) x[2]))
SkudEU$spIso <- unlist(lapply(split, function(x) x[3]))
SkudEU$pop <- unlist(lapply(split, function(x) x[4]))
SkudEU$spIso[which(SkudEU$strainID=="yHKS141")] <- "Skud"
SkudEU$pop[which(SkudEU$strainID=="yHKS141")] <- "Europe"
SkudEU$syn[which(SkudEU$strainID=="yHKS141")] <- NA
SkudEU$origin <- "Hybrid"
SkudEU$origin[which(SkudEU$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")

SkudEU$popLab <- SkudEU$pop
SkudEU$popLab[which(SkudEU$pop=="Europe")] <- "EU"
SkudEU$popLab[which(SkudEU$spIso=="CIDER")] <- "WineHybs"
SkudEU$popLab[which(SkudEU$spIso=="WINE")] <- "WineHybs"
SkudEU$popLab[which(SkudEU$spIso=="BEER")] <- "BeerHybs"
uniPops <- unique(SkudEU$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SkudEU[which(SkudEU$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
mid <- min(SkudEU$PC1)+((max(SkudEU$PC1)-min(SkudEU$PC1))/2)
spLabDF <- data.frame(x=mid, y=min(SkudEU$PC2), lab="Europe and Hybrids")
popLabMod <- popLabDF
popLabMod$pop[which(popLabMod$pop=="BeerHybs")] <- "Beer Hybrids"
popLabMod$pop[which(popLabMod$pop=="WineHybs")] <- "Wine Hybrids"
popLabMod$pop[which(popLabMod$pop=="EU")] <- "Europe"


####This is the best###
ggplot(SkudEU, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabMod, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=1,height=5), hjust = "inward") + 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=4),axis.title=element_text(size=6),text = element_text(size =4))

ggsave("SkudEU_PC1-PC2.jpg", height=2.5, width=3.5)

###SkudEUintrogressed
SkudEUIntro <- read.table("PCAscores_Skud_IntroEU.txt", check.names = F)
SkudStrainsEUIntro <-row.names(SkudEUIntro)
split <- strsplit(SkudStrainsEUIntro, "_")
SkudEUIntro$strainID <- unlist(lapply(split, function(x) x[1]))
SkudEUIntro$syn <- unlist(lapply(split, function(x) x[2]))
SkudEUIntro$spIso <- unlist(lapply(split, function(x) x[3]))
SkudEUIntro$pop <- unlist(lapply(split, function(x) x[4]))
SkudEUIntro$spIso[which(SkudEUIntro$strainID=="yHKS141")] <- "Skud"
SkudEUIntro$pop[which(SkudEUIntro$strainID=="yHKS141")] <- "Europe"
SkudEUIntro$syn[which(SkudEUIntro$strainID=="yHKS141")] <- NA
SkudEUIntro$origin <- "Hybrid"
SkudEUIntro$origin[which(SkudEUIntro$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")

SkudEUIntro$popLab <- SkudEUIntro$pop
SkudEUIntro$popLab[which(SkudEUIntro$pop=="Europe")] <- "EU"
SkudEUIntro$popLab[which(SkudEUIntro$spIso=="CIDER")] <- "WineHybs"
SkudEUIntro$popLab[which(SkudEUIntro$spIso=="WINE")] <- "WineHybs"
SkudEUIntro$popLab[which(SkudEUIntro$spIso=="BEER")] <- "BeerHybs"
uniPops <- unique(SkudEUIntro$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SkudEUIntro[which(SkudEUIntro$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSuvaXSeub" | popLabDF=="ScerXSkudXSeub"),]
mid <- min(SkudEUIntro$PC1)+((max(SkudEUIntro$PC1)-min(SkudEUIntro$PC1))/2)
spLabDF <- data.frame(x=mid, y=min(SkudEUIntro$PC2), lab="Europe and Hybrids\nWith minor introgressed hybrids")
popLabMod <- popLabTrim
popLabMod$pop[which(popLabMod$pop=="BeerHybs")] <- "Beer Hybrids"
popLabMod$pop[which(popLabMod$pop=="WineHybs")] <- "Wine Hybrids"
popLabMod$pop[which(popLabMod$pop=="EU")] <- "Europe"

####This is the best###
ggplot(SkudEUIntro, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabMod, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=1,height=5), hjust = "inward") + 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=4),axis.title=element_text(size=6),text = element_text(size =4))

ggsave("SkudEUintro_PC1-PC2.jpg", height=2.5, width=3.5)


###SeubHolarcitc
SeubSum <- read.table("Seub_PCAsummary.txt", check.names = F, header=T)
popColors <- c("PA1"="#404040", "PA2"="#6D6D6D", "PB2"="#8A8A8A", "PB1"="#A2A2A2",  "Holarctic"="#B6B6B6",  "HOL"="#B6B6B6", "Complex Hybrids"="#C7C7C7", "complexHyb"="#C7C7C7", "SuvaXSeub"="#C7C7C7", "Lager-Saaz"="#D7D7D7", "Lager-Frohberg" = "#E6E6E6", "Lager" = "#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))

SeubHol <- read.table("PCAscores_Seub_Holarctic.txt", check.names = F, sep="\t")
SeubStrainsHol <-row.names(SeubHol)
split <- strsplit(SeubStrainsHol, "_")
SeubHol$strainID <- unlist(lapply(split, function(x) x[1]))
SeubHol$syn <- unlist(lapply(split, function(x) x[2]))
SeubHol$spIso <- unlist(lapply(split, function(x) x[3]))
SeubHol$pop <- unlist(lapply(split, function(x) x[4]))
SeubHol$extra <- unlist(lapply(split, function(x) x[5]))
SeubHol$extra[which(SeubHol$extra=="1")] <- NA
SeubHol$extra[which(SeubHol$extra=="2")] <- NA
SeubHol$spIso[which(is.na(SeubHol$extra)==F)] <- SeubHol$pop[which(is.na(SeubHol$extra)==F)]
SeubHol$pop[which(is.na(SeubHol$extra)==F)] <- SeubHol$extra[which(is.na(SeubHol$extra)==F)]
SeubHol$origin <- "Hybrid"
SeubHol$origin[which(SeubHol$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")

SeubHol$popLab <- SeubHol$pop
SeubHol$popLab[which(SeubHol$pop=="lager-")] <- "lager-frohberg"
SeubHol$popLab[which(SeubHol$pop=="SuvaXSeub")] <- "complexHyb"
SeubHol$popLab[which(SeubHol$pop=="ScerXSkudXSeub")] <- "complexHyb"
uniPops <- unique(SeubHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SeubHol[which(SeubHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
#popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSeub"),]
mid <- min(SeubHol$PC1)+((max(SeubHol$PC1)-min(SeubHol$PC1))/2)
spLabDF <- data.frame(x=min(SeubHol$PC1), y=min(SeubHol$PC2), lab="Holarctic and Hybrids")
popLabMod <- popLabDF
popLabMod$pop[which(popLabMod$pop=="lager-saaz")] <- "Lager-Saaz"
popLabMod$pop[which(popLabMod$pop=="lager-frohberg")] <- "Lager-Frohberg"
popLabMod$pop[which(popLabMod$pop=="HOL")] <- "Holarctic"
popLabMod$pop[which(popLabMod$pop=="complexHyb")] <- "Complex Hybrids"

####This is the best###
ggplot(SeubHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabMod, aes(label = pop, x=PC1, y=PC2), alpha=.5, size=3, position=position_jitter(width=1,height=5), vjust="inward") + 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), vjust = "inward", size = 4, hjust = "inward") +
  guides(fill=FALSE)+theme(axis.text=element_text(size=7),axis.title=element_text(size=9),text = element_text(size =8))

ggsave("SeubHolPCA.jpg", height=5, width=7)

### IntroHol ###
SeubIntroHol <- read.table("PCAscores_Seub_IntroHol.txt", check.names = F, sep="\t")
SeubStrainsIntroHol <-row.names(SeubIntroHol)
split <- strsplit(SeubStrainsIntroHol, "_")
SeubIntroHol$strainID <- unlist(lapply(split, function(x) x[1]))
SeubIntroHol$syn <- unlist(lapply(split, function(x) x[2]))
SeubIntroHol$spIso <- unlist(lapply(split, function(x) x[3]))
SeubIntroHol$pop <- unlist(lapply(split, function(x) x[4]))
SeubIntroHol$extra <- unlist(lapply(split, function(x) x[5]))
SeubIntroHol$extra[which(SeubIntroHol$extra=="1")] <- NA
SeubIntroHol$extra[which(SeubIntroHol$extra=="2")] <- NA
SeubIntroHol$spIso[which(is.na(SeubIntroHol$extra)==F)] <- SeubIntroHol$pop[which(is.na(SeubIntroHol$extra)==F)]
SeubIntroHol$pop[which(is.na(SeubIntroHol$extra)==F)] <- SeubIntroHol$extra[which(is.na(SeubIntroHol$extra)==F)]
SeubIntroHol$origin <- "Hybrid"
SeubIntroHol$origin[which(SeubIntroHol$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")

SeubIntroHol$popLab <- SeubIntroHol$pop
SeubIntroHol$popLab[which(SeubIntroHol$pop=="lager-")] <- "lager-frohberg"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="ScerXSkudXSuvaXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="ScerXSkudXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="ScerXSuvaXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="SuvaXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(is.na(SeubIntroHol$pop))] <- "complexHyb"
uniPops <- unique(SeubIntroHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SeubIntroHol[which(SeubIntroHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
#mid <- min(SeubIntroHol$PC1)+((max(SeubIntroHol$PC1)-min(SeubIntroHol$PC1))/2)
spLabDF <- data.frame(x=min(SeubIntroHol$PC1), y=min(SeubIntroHol$PC2), lab="Holarctic and Hybrids\nWith minor introgressed hybrids")
popLabMod <- popLabDF
popLabMod$pop[which(popLabMod$pop=="lager-saaz")] <- "Lager-Saaz"
popLabMod$pop[which(popLabMod$pop=="lager-frohberg")] <- "Lager-Frohberg"
popLabMod$pop[which(popLabMod$pop=="HOL")] <- "Holarctic"
popLabMod$pop[which(popLabMod$pop=="complexHyb")] <- "Complex Hybrids"

####This is the best###
ggplot(SeubIntroHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabMod, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=1,height=-2), hjust = "inward") + 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), hjust = "inward", vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=7),axis.title=element_text(size=9),text = element_text(size =8))

ggsave("SeubIntroHolPCA.jpg", height=5, width=7)

###Suva Supp###
SuvaSum <- read.table("Suva_PCAsummary.txt", check.names = F, header=T)
popColors <- c("Aust"="#404040", "SA-B"="#6D6D6D", "SA-A"="#C1C1C1", "Holarctic"="#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))
SuvaHol <- read.table("PCAscores_Suva_holarc.txt", check.names = F, sep="\t")
SuvaStrainsHol <-row.names(SuvaHol)
split <- strsplit(SuvaStrainsHol, "_")
SuvaHol$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaHol$syn <- unlist(lapply(split, function(x) x[2]))
SuvaHol$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaHol$pop <- unlist(lapply(split, function(x) x[4]))
SuvaHol$origin <- "Hybrid"
SuvaHol$origin[which(SuvaHol$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")

SuvaHol$popLab <- "Holarctic"
#SuvaHol$popLab[which(SuvaHol$pop=="SuvaXSeub")] <- "SuvaXSeub"
#SuvaHol$popLab[which(SuvaHol$pop=="ScerXSuvaXSeub")] <- "complexHyb"
#SuvaHol$popLab[which(SuvaHol$pop=="ScerXSkudXSuvaXSeub")] <- "complexHyb"
uniPops <- unique(SuvaHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SuvaHol[which(SuvaHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
spLabDF <- data.frame(x=max(SuvaHol$PC1), y=max(SuvaHol$PC2), lab="Holarctic and Hybrids")

####This is the best###
ggplot(SuvaHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=1,height=-2), hjust = "inward") + 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), hjust = "inward", vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=7),axis.title=element_text(size=9),text = element_text(size =6))

ggsave("SuvaHolPCA.jpg", height=5, width=7)

### Introgressed Holarc ###
SuvaHolIntro <- read.table("PCAscores_Suva_IntroHol.txt", check.names = F, sep="\t")
SuvaStrainsHolIntro <-row.names(SuvaHolIntro)
split <- strsplit(SuvaStrainsHolIntro, "_")
SuvaHolIntro$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaHolIntro$syn <- unlist(lapply(split, function(x) x[2]))
SuvaHolIntro$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaHolIntro$pop <- unlist(lapply(split, function(x) x[4]))
SuvaHolIntro$origin <- "Hybrid"
SuvaHolIntro$origin[which(SuvaHolIntro$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")

SuvaHolIntro$popLab <- "Holarctic"
uniPops <- unique(SuvaHolIntro$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SuvaHolIntro[which(SuvaHolIntro$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
spLabDF <- data.frame(x=max(SuvaHolIntro$PC1), y=max(SuvaHolIntro$PC2), lab="Holarctic and Hybrids\nWith minor introgressed hybrids")

####This is the best###
ggplot(SuvaHolIntro, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, alpha=.5, size=3, position=position_jitter(width=1,height=-2), hjust = "inward") + 
  popFill + theme_bw() + labs(x = PC1var, y=PC2var) + 
  geom_text(data=spLabDF, aes(label=lab, x=x, y=y), hjust = "inward", vjust = "inward", size = 4) +
  guides(fill=FALSE)+theme(axis.text=element_text(size=7),axis.title=element_text(size=9),text = element_text(size =6))

ggsave("SuvaIntroHolPCA.jpg", height=5, width=7)

########### Skud ############
SkudSum <- read.table("Skud_PCAsummary.txt", check.names = F, header=T)
SkudAll <- read.table("PCAscores_Skud_All.txt", check.names = F)
SkudStrainsAll <-row.names(SkudAll)
split <- strsplit(SkudStrainsAll, "_")
SkudAll$strainID <- unlist(lapply(split, function(x) x[1]))
SkudAll$syn <- unlist(lapply(split, function(x) x[2]))
SkudAll$spIso <- unlist(lapply(split, function(x) x[3]))
SkudAll$pop <- unlist(lapply(split, function(x) x[4]))
SkudAll$spIso[which(SkudAll$strainID=="yHKS141")] <- "Skud"
SkudAll$pop[which(SkudAll$strainID=="yHKS141")] <- "Europe"
SkudAll$syn[which(SkudAll$strainID=="yHKS141")] <- NA
SkudAll$origin <- "Hybrid"
SkudAll$origin[which(SkudAll$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="All")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="All")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="All")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="All")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="All")]*100, 2), "%", sep="")

pdf("Skud_PCAs.pdf", width=10)


ggplot(SkudAll) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -20, check_overlap = T) + 
  ggtitle("Skud PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SkudAll) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: All strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SkudAll) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: All strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SkudAll) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="\n")), vjust = 0, check_overlap = T) + 
  ggtitle("Skud PCA: All strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)


SkudEUasiaA <- read.table("PCAscores_Skud_EUasiaA.txt", check.names = F)
SkudStrainsEUasiaA <-row.names(SkudEUasiaA)
split <- strsplit(SkudStrainsEUasiaA, "_")
SkudEUasiaA$strainID <- unlist(lapply(split, function(x) x[1]))
SkudEUasiaA$syn <- unlist(lapply(split, function(x) x[2]))
SkudEUasiaA$spIso <- unlist(lapply(split, function(x) x[3]))
SkudEUasiaA$pop <- unlist(lapply(split, function(x) x[4]))
SkudEUasiaA$spIso[which(SkudEUasiaA$strainID=="yHKS141")] <- "Skud"
SkudEUasiaA$pop[which(SkudEUasiaA$strainID=="yHKS141")] <- "EU"
SkudEUasiaA$syn[which(SkudEUasiaA$strainID=="yHKS141")] <- NA
SkudEUasiaA$origin <- "Hybrid"
SkudEUasiaA$origin[which(SkudEUasiaA$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="EUasiaA")]*100, 2), "%", sep="")

SkudEUasiaA$popLab <- SkudEUasiaA$pop
SkudEUasiaA$popLab[which(SkudEUasiaA$pop=="Europe")] <- "EU"
SkudEUasiaA$popLab[which(SkudEUasiaA$spIso=="CIDER")] <- "WineHybs"
SkudEUasiaA$popLab[which(SkudEUasiaA$spIso=="WINE")] <- "WineHybs"
SkudEUasiaA$popLab[which(SkudEUasiaA$spIso=="BEER")] <- "BeerHybs"
uniPops <- unique(SkudEUasiaA$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SkudEUasiaA[which(SkudEUasiaA$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
#popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSkud"| popLabDF$pop=="lager-saaz"),]

####This is the best###
ggplot(SkudEUasiaA, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -5, alpha=.5) + 
  ggtitle("Skud PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SkudEUasiaA_PC1-PC2.pdf", height=8, width=12)


ggplot(SkudEUasiaA) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EUasiaA strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SkudEUasiaA) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EUasiaA strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SkudEUasiaA) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EUasiaA strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SkudEUasiaA) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EUasiaA strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)


##Skud EU##
SkudEU <- read.table("PCAscores_Skud_EU.txt", check.names = F)
SkudStrainsEU <-row.names(SkudEU)
split <- strsplit(SkudStrainsEU, "_")
SkudEU$strainID <- unlist(lapply(split, function(x) x[1]))
SkudEU$syn <- unlist(lapply(split, function(x) x[2]))
SkudEU$spIso <- unlist(lapply(split, function(x) x[3]))
SkudEU$pop <- unlist(lapply(split, function(x) x[4]))
SkudEU$spIso[which(SkudEU$strainID=="yHKS141")] <- "Skud"
SkudEU$pop[which(SkudEU$strainID=="yHKS141")] <- "Europe"
SkudEU$syn[which(SkudEU$strainID=="yHKS141")] <- NA
SkudEU$origin <- "Hybrid"
SkudEU$origin[which(SkudEU$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="EU")]*100, 2), "%", sep="")

SkudEU$popLab <- SkudEU$pop
SkudEU$popLab[which(SkudEU$pop=="Europe")] <- "EU"
SkudEU$popLab[which(SkudEU$spIso=="CIDER")] <- "WineHybs"
SkudEU$popLab[which(SkudEU$spIso=="WINE")] <- "WineHybs"
SkudEU$popLab[which(SkudEU$spIso=="BEER")] <- "BeerHybs"
uniPops <- unique(SkudEU$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SkudEU[which(SkudEU$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}

####This is the best###
ggplot(SkudEU, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -5, alpha=.5) + 
  ggtitle("Seub Holarctic PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SkudEU-PCA.pdf", height=8, width=12)


ggplot(SkudEU) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EU strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SkudEU) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EU strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SkudEU) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EU strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SkudEU) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EU strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

SkudAllIntro <- read.table("PCAscores_Skud_Intro.txt", check.names = F)
SkudStrainsAllIntro <-row.names(SkudAllIntro)
split <- strsplit(SkudStrainsAllIntro, "_")
SkudAllIntro$strainID <- unlist(lapply(split, function(x) x[1]))
SkudAllIntro$syn <- unlist(lapply(split, function(x) x[2]))
SkudAllIntro$spIso <- unlist(lapply(split, function(x) x[3]))
SkudAllIntro$pop <- unlist(lapply(split, function(x) x[4]))
SkudAllIntro$spIso[which(SkudAllIntro$strainID=="yHKS141")] <- "Skud"
SkudAllIntro$pop[which(SkudAllIntro$strainID=="yHKS141")] <- "Europe"
SkudAllIntro$syn[which(SkudAllIntro$strainID=="yHKS141")] <- NA
SkudAllIntro$origin <- "Hybrid"
SkudAllIntro$origin[which(SkudAllIntro$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="Intro")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="Intro")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="Intro")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="Intro")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="Intro")]*100, 2), "%", sep="")

ggplot(SkudAllIntro) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -20, check_overlap = T) + 
  ggtitle("Skud PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SkudAllIntro) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SkudAllIntro) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SkudAllIntro) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

SkudEUIntro <- read.table("PCAscores_Skud_IntroEU.txt", check.names = F)
SkudStrainsEUIntro <-row.names(SkudEUIntro)
split <- strsplit(SkudStrainsEUIntro, "_")
SkudEUIntro$strainID <- unlist(lapply(split, function(x) x[1]))
SkudEUIntro$syn <- unlist(lapply(split, function(x) x[2]))
SkudEUIntro$spIso <- unlist(lapply(split, function(x) x[3]))
SkudEUIntro$pop <- unlist(lapply(split, function(x) x[4]))
SkudEUIntro$spIso[which(SkudEUIntro$strainID=="yHKS141")] <- "Skud"
SkudEUIntro$pop[which(SkudEUIntro$strainID=="yHKS141")] <- "Europe"
SkudEUIntro$syn[which(SkudEUIntro$strainID=="yHKS141")] <- NA
SkudEUIntro$origin <- "Hybrid"
SkudEUIntro$origin[which(SkudEUIntro$spIso=="Skud")] <- "Wild"
PC1var <- paste("PC1 = ", round(SkudSum$PC1[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SkudSum$PC2[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SkudSum$PC3[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SkudSum$PC4[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SkudSum$PC5[which(SkudSum$analysis=="IntroEU")]*100, 2), "%", sep="")

SkudEUIntro$popLab <- SkudEUIntro$pop
SkudEUIntro$popLab[which(SkudEUIntro$pop=="Europe")] <- "EU"
SkudEUIntro$popLab[which(SkudEUIntro$spIso=="CIDER")] <- "WineHybs"
SkudEUIntro$popLab[which(SkudEUIntro$spIso=="WINE")] <- "WineHybs"
SkudEUIntro$popLab[which(SkudEUIntro$spIso=="BEER")] <- "BeerHybs"
uniPops <- unique(SkudEUIntro$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SkudEUIntro[which(SkudEUIntro$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSuvaXSeub" | popLabDF=="ScerXSkudXSeub"),]

####This is the best###
ggplot(SkudEUIntro, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabTrim, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -5, alpha=.5) + 
  ggtitle("Seub Holarctic PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SkudEUIntro-PCA.pdf", height=8, width=12)

ggplot(SkudEUIntro) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EUIntro strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SkudEUIntro) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID, spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EUIntro strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SkudEUIntro) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=24, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Skud PCA: EUIntro strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

dev.off()

###### Suva #########
SuvaSum <- read.table("Suva_PCAsummary.txt", check.names = F, header=T)
popColors <- c("Aust"="#404040", "SA-B"="#6D6D6D", "SA-A"="#C1C1C1", "Holarctic"="#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))

### All ###
SuvaAll <- read.table("PCAscores_Suva_All.txt", check.names = F, sep="\t")
SuvaStrainsAll <-row.names(SuvaAll)
split <- strsplit(SuvaStrainsAll, "_")
SuvaAll$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaAll$syn <- unlist(lapply(split, function(x) x[2]))
SuvaAll$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaAll$pop <- unlist(lapply(split, function(x) x[4]))
SuvaAll$origin <- "Hybrid"
SuvaAll$origin[which(SuvaAll$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="All")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="All")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="All")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="All")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="All")]*100, 2), "%", sep="")

pdf("Suva_PCAs.pdf", width=10)

ggplot(SuvaAll) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SuvaAll) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: All strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SuvaAll) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: All strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SuvaAll) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: All strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### SoAm Hol ###
SuvaSoAmHol <- read.table("PCAscores_Suva_SoAmHolarc.txt", check.names = F, sep="\t")
SuvaStrainsSoAmHol <-row.names(SuvaSoAmHol)
split <- strsplit(SuvaStrainsSoAmHol, "_")
SuvaSoAmHol$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaSoAmHol$syn <- unlist(lapply(split, function(x) x[2]))
SuvaSoAmHol$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaSoAmHol$pop <- unlist(lapply(split, function(x) x[4]))
SuvaSoAmHol$origin <- "Hybrid"
SuvaSoAmHol$origin[which(SuvaSoAmHol$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="SBSAhol")]*100, 2), "%", sep="")

SuvaSoAmHol$popLab <- SuvaSoAmHol$pop
SuvaSoAmHol$popLab[-which(SuvaSoAmHol$pop=="SA-A" | SuvaSoAmHol$pop=="SA-B")] <- "Holarctic"
uniPops <- unique(SuvaSoAmHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SuvaSoAmHol[which(SuvaSoAmHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}

####This is the best###
ggplot(SuvaSoAmHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -10, alpha=.5) + 
  ggtitle("Suva PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SuvaSoAmHol_PC1-PC2.pdf", height=8, width=12)


ggplot(SuvaSoAmHol) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: SoAmHol strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SuvaSoAmHol) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: SoAmHol strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SuvaSoAmHol) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: SoAmHol strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SuvaSoAmHol) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: SoAmHol strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### SA-A Hol ###
SuvaSAaHol <- read.table("PCAscores_Suva_SA-Aholarc.txt", check.names = F, sep="\t")
SuvaStrainsSAaHol <-row.names(SuvaSAaHol)
split <- strsplit(SuvaStrainsSAaHol, "_")
SuvaSAaHol$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaSAaHol$syn <- unlist(lapply(split, function(x) x[2]))
SuvaSAaHol$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaSAaHol$pop <- unlist(lapply(split, function(x) x[4]))
SuvaSAaHol$origin <- "Hybrid"
SuvaSAaHol$origin[which(SuvaSAaHol$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="SAhol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="SAhol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="SAhol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="SAhol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="SAhol")]*100, 2), "%", sep="")

ggplot(SuvaSAaHol) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -5, check_overlap = T) + 
  ggtitle("Suva PCA: SAaHol strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SuvaSAaHol) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: SAaHol strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SuvaSAaHol) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: SAaHol strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SuvaSAaHol) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: SAaHol strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### Hol ###
SuvaHol <- read.table("PCAscores_Suva_holarc.txt", check.names = F, sep="\t")
SuvaStrainsHol <-row.names(SuvaHol)
split <- strsplit(SuvaStrainsHol, "_")
SuvaHol$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaHol$syn <- unlist(lapply(split, function(x) x[2]))
SuvaHol$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaHol$pop <- unlist(lapply(split, function(x) x[4]))
SuvaHol$origin <- "Hybrid"
SuvaHol$origin[which(SuvaHol$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="Holarcitc")]*100, 2), "%", sep="")

SuvaHol$popLab <- "Holarctic"
#SuvaHol$popLab[which(SuvaHol$pop=="SuvaXSeub")] <- "SuvaXSeub"
#SuvaHol$popLab[which(SuvaHol$pop=="ScerXSuvaXSeub")] <- "complexHyb"
#SuvaHol$popLab[which(SuvaHol$pop=="ScerXSkudXSuvaXSeub")] <- "complexHyb"
uniPops <- unique(SuvaHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SuvaHol[which(SuvaHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}

####This is the best###
ggplot(SuvaHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -5, alpha=.5) + 
  ggtitle("Seub Holarctic PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SuvaHolPCA.pdf", height=8, width=12)


ggplot(SuvaHol) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: Hol strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SuvaHol) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: Hol strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SuvaHol) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: Hol strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SuvaHol) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: Hol strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### All Intro ###
SuvaAllIntro <- read.table("PCAscores_Suva_Intro.txt", check.names = F, sep="\t")
SuvaStrainsAllIntro <-row.names(SuvaAllIntro)
split <- strsplit(SuvaStrainsAllIntro, "_")
SuvaAllIntro$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaAllIntro$syn <- unlist(lapply(split, function(x) x[2]))
SuvaAllIntro$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaAllIntro$pop <- unlist(lapply(split, function(x) x[4]))
SuvaAllIntro$origin <- "Hybrid"
SuvaAllIntro$origin[which(SuvaAllIntro$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="Intro")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="Intro")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="Intro")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="Intro")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="Intro")]*100, 2), "%", sep="")

ggplot(SuvaAllIntro) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SuvaAllIntro) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SuvaAllIntro) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SuvaAllIntro) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: AllIntro strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### Introgressed Holarc ###
SuvaHolIntro <- read.table("PCAscores_Suva_IntroHol.txt", check.names = F, sep="\t")
SuvaStrainsHolIntro <-row.names(SuvaHolIntro)
split <- strsplit(SuvaStrainsHolIntro, "_")
SuvaHolIntro$strainID <- unlist(lapply(split, function(x) x[1]))
SuvaHolIntro$syn <- unlist(lapply(split, function(x) x[2]))
SuvaHolIntro$spIso <- unlist(lapply(split, function(x) x[3]))
SuvaHolIntro$pop <- unlist(lapply(split, function(x) x[4]))
SuvaHolIntro$origin <- "Hybrid"
SuvaHolIntro$origin[which(SuvaHolIntro$spIso=="Suva")] <- "Wild"
PC1var <- paste("PC1 = ", round(SuvaSum$PC1[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SuvaSum$PC2[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SuvaSum$PC3[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SuvaSum$PC4[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SuvaSum$PC5[which(SuvaSum$analysis=="IntroHol")]*100, 2), "%", sep="")

SuvaHolIntro$popLab <- "Holarctic"
#SuvaHolIntro$popLab[which(SuvaHolIntro$pop=="SuvaXSeub")] <- "SuvaXSeub"
#SuvaHolIntro$popLab[which(SuvaHolIntro$pop=="ScerXSuvaXSeub")] <- "complexHyb"
#SuvaHolIntro$popLab[which(SuvaHolIntro$pop=="ScerXSkudXSuvaXSeub")] <- "complexHyb"
uniPops <- unique(SuvaHolIntro$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SuvaHolIntro[which(SuvaHolIntro$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}

####This is the best###
ggplot(SuvaHolIntro, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -5, alpha=.5) + 
  ggtitle("Seub Holarctic PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SuvaHolIntroPCA.pdf", height=8, width=12)


ggplot(SuvaHolIntro) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: HolIntro strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SuvaHolIntro) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: HolIntro strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SuvaHolIntro) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: HolIntro strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SuvaHolIntro) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=22, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Suva PCA: HolIntro strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

dev.off()


###### Seub #########
SeubSum <- read.table("Seub_PCAsummary.txt", check.names = F, header=T)
popColors <- c("PA1"="#404040", "PA2"="#6D6D6D", "PB2"="#8A8A8A", "PB1"="#A2A2A2",  "HOL"="#B6B6B6", "ComplexHyb"="#C7C7C7", "complexHyb"="#C7C7C7", "SuvaXSeub"="#C7C7C7", "lager-saaz"="#D7D7D7", "lager-frohberg" = "#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))
### All ###
SeubAll <- read.table("PCAscores_Seub_All.txt", check.names = F, sep="\t")
SeubStrainsAll <-row.names(SeubAll)
split <- strsplit(SeubStrainsAll, "_")
SeubAll$strainID <- unlist(lapply(split, function(x) x[1]))
SeubAll$syn <- unlist(lapply(split, function(x) x[2]))
SeubAll$spIso <- unlist(lapply(split, function(x) x[3]))
SeubAll$pop <- unlist(lapply(split, function(x) x[4]))
SeubAll$extra <- unlist(lapply(split, function(x) x[5]))
SeubAll$extra[which(SeubAll$extra=="1")] <- NA
SeubAll$extra[which(SeubAll$extra=="2")] <- NA
SeubAll$spIso[which(is.na(SeubAll$extra)==F)] <- SeubAll$pop[which(is.na(SeubAll$extra)==F)]
SeubAll$pop[which(is.na(SeubAll$extra)==F)] <- SeubAll$extra[which(is.na(SeubAll$extra)==F)]
SeubAll$origin <- "Hybrid"
SeubAll$origin[which(SeubAll$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")

pdf("Seub_PCAs.pdf", width=10)

SeubAll$popLab <- SeubAll$pop
SeubAll$popLab[which(SeubAll$pop=="lager-")] <- "lager-frohberg"
SeubAll$popLab[which(SeubAll$pop=="SuvaXSeub")] <- "complexHyb"
SeubAll$popLab[which(SeubAll$pop=="SuvaXSkudXSeub")] <- "complexHyb"
uniPops <- unique(SeubAll$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SeubAll[which(SeubAll$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSeub"| popLabDF$pop=="lager-saaz"),]

####This is the best###
ggplot(SeubAll, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabTrim, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -10, alpha=.5) + 
  ggtitle("Seub PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SeubAll_PC1-PC2.pdf", height=8, width=12)


ggplot(SeubAll) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubAll) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubAll) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubAll) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### All LF ###
SeubAll_LF <- read.table("PCAscores_Seub_allPubLF.txt", check.names = F, sep="\t")
SeubStrainsAll_LF <-row.names(SeubAll_LF)
split <- strsplit(SeubStrainsAll_LF, "_")
SeubAll_LF$strainID <- unlist(lapply(split, function(x) x[1]))
SeubAll_LF$syn <- unlist(lapply(split, function(x) x[2]))
SeubAll_LF$spIso <- unlist(lapply(split, function(x) x[3]))
SeubAll_LF$pop <- unlist(lapply(split, function(x) x[4]))
SeubAll_LF$extra <- unlist(lapply(split, function(x) x[5]))
SeubAll_LF$extra[which(SeubAll_LF$extra=="1")] <- NA
SeubAll_LF$extra[which(SeubAll_LF$extra=="2")] <- NA
SeubAll_LF$spIso[which(is.na(SeubAll_LF$extra)==F)] <- SeubAll_LF$pop[which(is.na(SeubAll_LF$extra)==F)]
SeubAll_LF$pop[which(is.na(SeubAll_LF$extra)==F)] <- SeubAll_LF$extra[which(is.na(SeubAll_LF$extra)==F)]
SeubAll_LF$origin <- "Hybrid"
SeubAll_LF$origin[which(SeubAll_LF$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="All_LF")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="All_LF")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="All_LF")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="All_LF")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="All_LF")]*100, 2), "%", sep="")

ggplot(SeubAll_LF) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All_LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubAll_LF) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All_LF strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubAll_LF) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All_LF strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubAll_LF) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: All_LF strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### B and Holarc ###
SeubBholarctic <- read.table("PCAscores_Seub_PopBholarctic.txt", check.names = F, sep="\t")
SeubStrainsBholarctic <-row.names(SeubBholarctic)
split <- strsplit(SeubStrainsBholarctic, "_")
SeubBholarctic$strainID <- unlist(lapply(split, function(x) x[1]))
SeubBholarctic$syn <- unlist(lapply(split, function(x) x[2]))
SeubBholarctic$spIso <- unlist(lapply(split, function(x) x[3]))
SeubBholarctic$pop <- unlist(lapply(split, function(x) x[4]))
SeubBholarctic$extra <- unlist(lapply(split, function(x) x[5]))
SeubBholarctic$extra[which(SeubBholarctic$extra=="1")] <- NA
SeubBholarctic$extra[which(SeubBholarctic$extra=="2")] <- NA
SeubBholarctic$spIso[which(is.na(SeubBholarctic$extra)==F)] <- SeubBholarctic$pop[which(is.na(SeubBholarctic$extra)==F)]
SeubBholarctic$pop[which(is.na(SeubBholarctic$extra)==F)] <- SeubBholarctic$extra[which(is.na(SeubBholarctic$extra)==F)]
SeubBholarctic$origin <- "Hybrid"
SeubBholarctic$origin[which(SeubBholarctic$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="Bholarcitc")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="Bholarcitc")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="Bholarcitc")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="Bholarcitc")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="Bholarcitc")]*100, 2), "%", sep="")

ggplot(SeubBholarctic) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubBholarctic) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubBholarctic) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubBholarctic) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### B and Holarc less Frohberg ###
SeubBholarctic_LF <- read.table("PCAscores_Seub_Bholarctic_lessFrohberg.txt", check.names = F, sep="\t")
SeubStrainsBholarctic_LF <-row.names(SeubBholarctic_LF)
split <- strsplit(SeubStrainsBholarctic_LF, "_")
SeubBholarctic_LF$strainID <- unlist(lapply(split, function(x) x[1]))
SeubBholarctic_LF$syn <- unlist(lapply(split, function(x) x[2]))
SeubBholarctic_LF$spIso <- unlist(lapply(split, function(x) x[3]))
SeubBholarctic_LF$pop <- unlist(lapply(split, function(x) x[4]))
SeubBholarctic_LF$extra <- unlist(lapply(split, function(x) x[5]))
SeubBholarctic_LF$extra[which(SeubBholarctic_LF$extra=="1")] <- NA
SeubBholarctic_LF$extra[which(SeubBholarctic_LF$extra=="2")] <- NA
SeubBholarctic_LF$spIso[which(is.na(SeubBholarctic_LF$extra)==F)] <- SeubBholarctic_LF$pop[which(is.na(SeubBholarctic_LF$extra)==F)]
SeubBholarctic_LF$pop[which(is.na(SeubBholarctic_LF$extra)==F)] <- SeubBholarctic_LF$extra[which(is.na(SeubBholarctic_LF$extra)==F)]
SeubBholarctic_LF$origin <- "Hybrid"
SeubBholarctic_LF$origin[which(SeubBholarctic_LF$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="Bholarcitc_LF")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="Bholarcitc_LF")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="Bholarcitc_LF")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="Bholarcitc_LF")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="Bholarcitc_LF")]*100, 2), "%", sep="")

ggplot(SeubBholarctic_LF) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic_LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubBholarctic_LF) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic_LF strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubBholarctic_LF) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic_LF strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubBholarctic_LF) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Bholarctic_LF strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### Holarc ###
SeubHol <- read.table("PCAscores_Seub_Holarctic.txt", check.names = F, sep="\t")
SeubStrainsHol <-row.names(SeubHol)
split <- strsplit(SeubStrainsHol, "_")
SeubHol$strainID <- unlist(lapply(split, function(x) x[1]))
SeubHol$syn <- unlist(lapply(split, function(x) x[2]))
SeubHol$spIso <- unlist(lapply(split, function(x) x[3]))
SeubHol$pop <- unlist(lapply(split, function(x) x[4]))
SeubHol$extra <- unlist(lapply(split, function(x) x[5]))
SeubHol$extra[which(SeubHol$extra=="1")] <- NA
SeubHol$extra[which(SeubHol$extra=="2")] <- NA
SeubHol$spIso[which(is.na(SeubHol$extra)==F)] <- SeubHol$pop[which(is.na(SeubHol$extra)==F)]
SeubHol$pop[which(is.na(SeubHol$extra)==F)] <- SeubHol$extra[which(is.na(SeubHol$extra)==F)]
SeubHol$origin <- "Hybrid"
SeubHol$origin[which(SeubHol$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="All")]*100, 2), "%", sep="")

SeubHol$popLab <- SeubHol$pop
SeubHol$popLab[which(SeubHol$pop=="lager-")] <- "lager-frohberg"
SeubHol$popLab[which(SeubHol$pop=="SuvaXSeub")] <- "complexHyb"
SeubHol$popLab[which(SeubHol$pop=="ScerXSkudXSeub")] <- "complexHyb"
uniPops <- unique(SeubHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SeubHol[which(SeubHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}
popLabTrim <- popLabDF[-which(popLabDF$pop=="ScerXSkudXSeub"| popLabDF$pop=="lager-saaz"),]

####This is the best###
ggplot(SeubHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -5, alpha=.5) + 
  ggtitle("Seub Holarctic PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SeubHolPCA.pdf", height=8, width=12)

ggplot(SeubHol) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Hol strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubHol) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label =paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Hol strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubHol) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Hol strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

### Holarc less Frohberg ###
SeubHol_LF <- read.table("PCAscores_Seub_Holarctic_lessFrohberg.txt", check.names = F, sep="\t")
SeubStrainsHol_LF <-row.names(SeubHol_LF)
split <- strsplit(SeubStrainsHol_LF, "_")
SeubHol_LF$strainID <- unlist(lapply(split, function(x) x[1]))
SeubHol_LF$syn <- unlist(lapply(split, function(x) x[2]))
SeubHol_LF$spIso <- unlist(lapply(split, function(x) x[3]))
SeubHol_LF$pop <- unlist(lapply(split, function(x) x[4]))
SeubHol_LF$extra <- unlist(lapply(split, function(x) x[5]))
SeubHol_LF$extra[which(SeubHol_LF$extra=="1")] <- NA
SeubHol_LF$extra[which(SeubHol_LF$extra=="2")] <- NA
SeubHol_LF$spIso[which(is.na(SeubHol_LF$extra)==F)] <- SeubHol_LF$pop[which(is.na(SeubHol_LF$extra)==F)]
SeubHol_LF$pop[which(is.na(SeubHol_LF$extra)==F)] <- SeubHol_LF$extra[which(is.na(SeubHol_LF$extra)==F)]
SeubHol_LF$origin <- "Hybrid"
SeubHol_LF$origin[which(SeubHol_LF$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="Holarcitc_LF")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="Holarcitc_LF")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="Holarcitc_LF")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="Holarcitc_LF")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="Holarcitc_LF")]*100, 2), "%", sep="")

ggplot(SeubHol_LF) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Hol_LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubHol_LF) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label =paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Hol_LF strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubHol_LF) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Hol_LF strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubHol_LF) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Hol_LF strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### Intro ###
SeubIntro <- read.table("PCAscores_Seub_Intro.txt", check.names = F, sep="\t")
SeubStrainsIntro <-row.names(SeubIntro)
split <- strsplit(SeubStrainsIntro, "_")
SeubIntro$strainID <- unlist(lapply(split, function(x) x[1]))
SeubIntro$syn <- unlist(lapply(split, function(x) x[2]))
SeubIntro$spIso <- unlist(lapply(split, function(x) x[3]))
SeubIntro$pop <- unlist(lapply(split, function(x) x[4]))
SeubIntro$extra <- unlist(lapply(split, function(x) x[5]))
SeubIntro$extra[which(SeubIntro$extra=="1")] <- NA
SeubIntro$extra[which(SeubIntro$extra=="2")] <- NA
SeubIntro$spIso[which(is.na(SeubIntro$extra)==F)] <- SeubIntro$pop[which(is.na(SeubIntro$extra)==F)]
SeubIntro$pop[which(is.na(SeubIntro$extra)==F)] <- SeubIntro$extra[which(is.na(SeubIntro$extra)==F)]
SeubIntro$origin <- "Hybrid"
SeubIntro$origin[which(SeubIntro$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="Intro")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="Intro")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="Intro")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="Intro")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="Intro")]*100, 2), "%", sep="")

ggplot(SeubIntro) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubIntro) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubIntro) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubIntro) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### Intro_LF ###
SeubIntro_LF <- read.table("PCAscores_Seub_intro_lessFrohberg.txt", check.names = F, sep="\t")
SeubStrainsIntro_LF <-row.names(SeubIntro_LF)
split <- strsplit(SeubStrainsIntro_LF, "_")
SeubIntro_LF$strainID <- unlist(lapply(split, function(x) x[1]))
SeubIntro_LF$syn <- unlist(lapply(split, function(x) x[2]))
SeubIntro_LF$spIso <- unlist(lapply(split, function(x) x[3]))
SeubIntro_LF$pop <- unlist(lapply(split, function(x) x[4]))
SeubIntro_LF$extra <- unlist(lapply(split, function(x) x[5]))
SeubIntro_LF$extra[which(SeubIntro_LF$extra=="1")] <- NA
SeubIntro_LF$extra[which(SeubIntro_LF$extra=="2")] <- NA
SeubIntro_LF$spIso[which(is.na(SeubIntro_LF$extra)==F)] <- SeubIntro_LF$pop[which(is.na(SeubIntro_LF$extra)==F)]
SeubIntro_LF$pop[which(is.na(SeubIntro_LF$extra)==F)] <- SeubIntro_LF$extra[which(is.na(SeubIntro_LF$extra)==F)]
SeubIntro_LF$origin <- "Hybrid"
SeubIntro_LF$origin[which(SeubIntro_LF$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="Intro_LF")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="Intro_LF")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="Intro_LF")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="Intro_LF")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="Intro_LF")]*100, 2), "%", sep="")

ggplot(SeubIntro_LF) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro_LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubIntro_LF) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro_LF strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubIntro_LF) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro_LF strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubIntro_LF) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Intro_LF strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### IntroHol ###
SeubIntroHol <- read.table("PCAscores_Seub_IntroHol.txt", check.names = F, sep="\t")
SeubStrainsIntroHol <-row.names(SeubIntroHol)
split <- strsplit(SeubStrainsIntroHol, "_")
SeubIntroHol$strainID <- unlist(lapply(split, function(x) x[1]))
SeubIntroHol$syn <- unlist(lapply(split, function(x) x[2]))
SeubIntroHol$spIso <- unlist(lapply(split, function(x) x[3]))
SeubIntroHol$pop <- unlist(lapply(split, function(x) x[4]))
SeubIntroHol$extra <- unlist(lapply(split, function(x) x[5]))
SeubIntroHol$extra[which(SeubIntroHol$extra=="1")] <- NA
SeubIntroHol$extra[which(SeubIntroHol$extra=="2")] <- NA
SeubIntroHol$spIso[which(is.na(SeubIntroHol$extra)==F)] <- SeubIntroHol$pop[which(is.na(SeubIntroHol$extra)==F)]
SeubIntroHol$pop[which(is.na(SeubIntroHol$extra)==F)] <- SeubIntroHol$extra[which(is.na(SeubIntroHol$extra)==F)]
SeubIntroHol$origin <- "Hybrid"
SeubIntroHol$origin[which(SeubIntroHol$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="IntroHol")]*100, 2), "%", sep="")

SeubIntroHol$popLab <- SeubIntroHol$pop
SeubIntroHol$popLab[which(SeubIntroHol$pop=="lager-")] <- "lager-frohberg"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="ScerXSkudXSuvaXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="ScerXSkudXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="ScerXSuvaXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(SeubIntroHol$pop=="SuvaXSeub")] <- "complexHyb"
SeubIntroHol$popLab[which(is.na(SeubIntroHol$pop))] <- "complexHyb"
uniPops <- unique(SeubIntroHol$popLab)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- SeubIntroHol[which(SeubIntroHol$popLab==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2)
  popLabDF <- rbind(popLabDF, labTemp)
}

####This is the best###
ggplot(SeubIntroHol, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=popLab, fill=popLab), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = -5, alpha=.5) + 
  ggtitle("Seub Holarctic PCA: Introgressed strains included") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("SeubIntroHolPCA.pdf", height=8, width=12)


ggplot(SeubIntroHol) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubIntroHol) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubIntroHol) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubIntroHol) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### IntroHol_LF ###
SeubIntroHol_LF <- read.table("PCAscores_Seub_introHol_lessFrohberg.txt", check.names = F, sep="\t")
SeubStrainsIntroHol_LF <-row.names(SeubIntroHol_LF)
split <- strsplit(SeubStrainsIntroHol_LF, "_")
SeubIntroHol_LF$strainID <- unlist(lapply(split, function(x) x[1]))
SeubIntroHol_LF$syn <- unlist(lapply(split, function(x) x[2]))
SeubIntroHol_LF$spIso <- unlist(lapply(split, function(x) x[3]))
SeubIntroHol_LF$pop <- unlist(lapply(split, function(x) x[4]))
SeubIntroHol_LF$extra <- unlist(lapply(split, function(x) x[5]))
SeubIntroHol_LF$extra[which(SeubIntroHol_LF$extra=="1")] <- NA
SeubIntroHol_LF$extra[which(SeubIntroHol_LF$extra=="2")] <- NA
SeubIntroHol_LF$spIso[which(is.na(SeubIntroHol_LF$extra)==F)] <- SeubIntroHol_LF$pop[which(is.na(SeubIntroHol_LF$extra)==F)]
SeubIntroHol_LF$pop[which(is.na(SeubIntroHol_LF$extra)==F)] <- SeubIntroHol_LF$extra[which(is.na(SeubIntroHol_LF$extra)==F)]
SeubIntroHol_LF$origin <- "Hybrid"
SeubIntroHol_LF$origin[which(SeubIntroHol_LF$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="IntroHol_LF")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="IntroHol_LF")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="IntroHol_LF")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="IntroHol_LF")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="IntroHol_LF")]*100, 2), "%", sep="")

ggplot(SeubIntroHol_LF) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol_LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubIntroHol_LF) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol_LF strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubIntroHol_LF) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol_LF strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubIntroHol_LF) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Seub PCA: IntroHol_LF strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### FullHybs ###
SeubFullHybs <- read.table("PCAscores_Seub_FullHybs.txt", check.names = F, sep="\t")
SeubStrainsFullHybs <-row.names(SeubFullHybs)
split <- strsplit(SeubStrainsFullHybs, "_")
SeubFullHybs$strainID <- unlist(lapply(split, function(x) x[1]))
SeubFullHybs$syn <- unlist(lapply(split, function(x) x[2]))
SeubFullHybs$spIso <- unlist(lapply(split, function(x) x[3]))
SeubFullHybs$pop <- unlist(lapply(split, function(x) x[4]))
SeubFullHybs$extra <- unlist(lapply(split, function(x) x[5]))
SeubFullHybs$extra[which(SeubFullHybs$extra=="1")] <- NA
SeubFullHybs$extra[which(SeubFullHybs$extra=="2")] <- NA
SeubFullHybs$spIso[which(is.na(SeubFullHybs$extra)==F)] <- SeubFullHybs$pop[which(is.na(SeubFullHybs$extra)==F)]
SeubFullHybs$pop[which(is.na(SeubFullHybs$extra)==F)] <- SeubFullHybs$extra[which(is.na(SeubFullHybs$extra)==F)]
SeubFullHybs$pop[which(is.na(SeubFullHybs$spIso)==T)] <- SeubFullHybs$syn[which(is.na(SeubFullHybs$spIso)==T)]
toChange <- which(is.na(SeubFullHybs$pop)==T)
SeubFullHybs$pop[toChange] <- SeubFullHybs$spIso[toChange]
SeubFullHybs$spIso[toChange] <- SeubFullHybs$syn[toChange]
SeubFullHybs$spIso[which(SeubFullHybs$spIso=="Admix")] <- "Seub"
SeubFullHybs$origin <- "Hybrid"
SeubFullHybs$origin[which(SeubFullHybs$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="FullHybs")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="FullHybs")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="FullHybs")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="FullHybs")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="FullHybs")]*100, 2), "%", sep="")

ggplot(SeubFullHybs) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullHybs strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubFullHybs) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullHybs strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubFullHybs) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullHybs strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubFullHybs) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullHybs strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### FullNoAdHybs ###
SeubFullNoAdHybs <- read.table("PCAscores_Seub_FullNoAdHybs.txt", check.names = F, sep="\t")
SeubStrainsFullNoAdHybs <-row.names(SeubFullNoAdHybs)
split <- strsplit(SeubStrainsFullNoAdHybs, "_")
SeubFullNoAdHybs$strainID <- unlist(lapply(split, function(x) x[1]))
SeubFullNoAdHybs$syn <- unlist(lapply(split, function(x) x[2]))
SeubFullNoAdHybs$spIso <- unlist(lapply(split, function(x) x[3]))
SeubFullNoAdHybs$pop <- unlist(lapply(split, function(x) x[4]))
SeubFullNoAdHybs$extra <- unlist(lapply(split, function(x) x[5]))
SeubFullNoAdHybs$extra[which(SeubFullNoAdHybs$extra=="1")] <- NA
SeubFullNoAdHybs$extra[which(SeubFullNoAdHybs$extra=="2")] <- NA
SeubFullNoAdHybs$spIso[which(is.na(SeubFullNoAdHybs$extra)==F)] <- SeubFullNoAdHybs$pop[which(is.na(SeubFullNoAdHybs$extra)==F)]
SeubFullNoAdHybs$pop[which(is.na(SeubFullNoAdHybs$extra)==F)] <- SeubFullNoAdHybs$extra[which(is.na(SeubFullNoAdHybs$extra)==F)]
SeubFullNoAdHybs$pop[which(is.na(SeubFullNoAdHybs$spIso)==T)] <- SeubFullNoAdHybs$syn[which(is.na(SeubFullNoAdHybs$spIso)==T)]
toChange <- which(is.na(SeubFullNoAdHybs$pop)==T)
SeubFullNoAdHybs$pop[toChange] <- SeubFullNoAdHybs$spIso[toChange]
SeubFullNoAdHybs$spIso[toChange] <- SeubFullNoAdHybs$syn[toChange]
SeubFullNoAdHybs$spIso[which(SeubFullNoAdHybs$spIso=="Admix")] <- "Seub"
SeubFullNoAdHybs$origin <- "Hybrid"
SeubFullNoAdHybs$origin[which(SeubFullNoAdHybs$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="FullNoAdHybs")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="FullNoAdHybs")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="FullNoAdHybs")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="FullNoAdHybs")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="FullNoAdHybs")]*100, 2), "%", sep="")

ggplot(SeubFullNoAdHybs) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAdHybs strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubFullNoAdHybs) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAdHybs strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubFullNoAdHybs) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAdHybs strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubFullNoAdHybs) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAdHybs strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### Full ###
SeubFull <- read.table("PCAscores_Seub_Full.txt", check.names = F, sep="\t")
SeubStrainsFull <-row.names(SeubFull)
split <- strsplit(SeubStrainsFull, "_")
SeubFull$strainID <- unlist(lapply(split, function(x) x[1]))
SeubFull$syn <- unlist(lapply(split, function(x) x[2]))
SeubFull$spIso <- unlist(lapply(split, function(x) x[3]))
SeubFull$pop <- unlist(lapply(split, function(x) x[4]))
SeubFull$extra <- unlist(lapply(split, function(x) x[5]))
SeubFull$extra[which(SeubFull$extra=="1")] <- NA
SeubFull$extra[which(SeubFull$extra=="2")] <- NA
SeubFull$spIso[which(is.na(SeubFull$extra)==F)] <- SeubFull$pop[which(is.na(SeubFull$extra)==F)]
SeubFull$pop[which(is.na(SeubFull$extra)==F)] <- SeubFull$extra[which(is.na(SeubFull$extra)==F)]
SeubFull$pop[which(is.na(SeubFull$spIso)==T)] <- SeubFull$syn[which(is.na(SeubFull$spIso)==T)]
toChange <- which(is.na(SeubFull$pop)==T)
SeubFull$pop[toChange] <- SeubFull$spIso[toChange]
SeubFull$spIso[toChange] <- SeubFull$syn[toChange]
SeubFull$spIso[which(SeubFull$spIso=="Admix")] <- "Seub"
SeubFull$origin <- "Hybrid"
SeubFull$origin[which(SeubFull$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="Full")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="Full")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="Full")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="Full")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="Full")]*100, 2), "%", sep="")

ggplot(SeubFull) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Full strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubFull) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Full strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubFull) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Full strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubFull) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: Full strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

### FullNoAd ###
SeubFullNoAd <- read.table("PCAscores_Seub_FullNoAd.txt", check.names = F, sep="\t")
SeubStrainsFullNoAd <-row.names(SeubFullNoAd)
split <- strsplit(SeubStrainsFullNoAd, "_")
SeubFullNoAd$strainID <- unlist(lapply(split, function(x) x[1]))
SeubFullNoAd$syn <- unlist(lapply(split, function(x) x[2]))
SeubFullNoAd$spIso <- unlist(lapply(split, function(x) x[3]))
SeubFullNoAd$pop <- unlist(lapply(split, function(x) x[4]))
SeubFullNoAd$extra <- unlist(lapply(split, function(x) x[5]))
SeubFullNoAd$extra[which(SeubFullNoAd$extra=="1")] <- NA
SeubFullNoAd$extra[which(SeubFullNoAd$extra=="2")] <- NA
SeubFullNoAd$spIso[which(is.na(SeubFullNoAd$extra)==F)] <- SeubFullNoAd$pop[which(is.na(SeubFullNoAd$extra)==F)]
SeubFullNoAd$pop[which(is.na(SeubFullNoAd$extra)==F)] <- SeubFullNoAd$extra[which(is.na(SeubFullNoAd$extra)==F)]
SeubFullNoAd$pop[which(is.na(SeubFullNoAd$spIso)==T)] <- SeubFullNoAd$syn[which(is.na(SeubFullNoAd$spIso)==T)]
toChange <- which(is.na(SeubFullNoAd$pop)==T)
SeubFullNoAd$pop[toChange] <- SeubFullNoAd$spIso[toChange]
SeubFullNoAd$spIso[toChange] <- SeubFullNoAd$syn[toChange]
SeubFullNoAd$spIso[which(SeubFullNoAd$spIso=="Admix")] <- "Seub"
SeubFullNoAd$origin <- "Hybrid"
SeubFullNoAd$origin[which(SeubFullNoAd$spIso=="Seub")] <- "Wild"
PC1var <- paste("PC1 = ", round(SeubSum$PC1[which(SeubSum$analysis=="FullNoAd")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(SeubSum$PC2[which(SeubSum$analysis=="FullNoAd")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(SeubSum$PC3[which(SeubSum$analysis=="FullNoAd")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(SeubSum$PC4[which(SeubSum$analysis=="FullNoAd")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(SeubSum$PC5[which(SeubSum$analysis=="FullNoAd")]*100, 2), "%", sep="")

ggplot(SeubFullNoAd) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAd strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(SeubFullNoAd) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID, spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAd strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(SeubFullNoAd) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAd strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(SeubFullNoAd) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=23, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID, spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Seub PCA: FullNoAd strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

dev.off()


###### Scer #########
renameOri <- read.table("strainPopOriKey.txt", check.names = F, header = T)
ScerSum <- read.table("Scer_PCAsummary.txt", check.names = F, header=T)
popColors <- c("WildMisc"="#404040", "SakeAsia"="#6D6D6D", "Mosaic"="#8A8A8A", "BreadMixed"="#A2A2A2",  "Beer2"="#B6B6B6", "MedOak"="#C7C7C7", "Wine"="#D7D7D7", "AleBeer1" = "#E6E6E6")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))
### All ###
Scer.galGonAll <- read.table("PCAscores_Scer_galGonAll.txt", check.names = F, sep="\t")
Scer.galGonAll.strains <-row.names(Scer.galGonAll)
split <- strsplit(Scer.galGonAll.strains, "_")
Scer.galGonAll$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.galGonAll$syn <- unlist(lapply(split, function(x) x[2]))
Scer.galGonAll$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.galGonAll$tempPop <- unlist(lapply(split, function(x) x[4]))
Scer.galGonAll$phase <- unlist(lapply(split, function(x) x[5]))
Scer.galGonAll$phase[which(Scer.galGonAll$pop==1)] <- 1
Scer.galGonAll$phase[which(Scer.galGonAll$pop==2)] <- 2
Scer.galGonAll$tempPop[which(Scer.galGonAll$tempPop==1)] <- "Unplaced"
Scer.galGonAll$tempPop[which(Scer.galGonAll$tempPop==2)] <- "Unplaced"
Scer.galGonAll$origin <- NA
Scer.galGonAll$pop <- NA
for (strainID in Scer.galGonAll.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  strainPop <- as.character(renameOri$pop[which(renameOri$rename==strainID)])
  Scer.galGonAll$origin[which(rownames(Scer.galGonAll)==strainID)] <- strainOri
  Scer.galGonAll$pop[which(rownames(Scer.galGonAll)==strainID)] <- strainPop
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="galGonAll")]*100, 2), "%", sep="")

uniPops <- unique(Scer.galGonAll$pop)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- Scer.galGonAll[which(Scer.galGonAll$pop==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF <- rbind(popLabDF, labTemp)
}

pdf("Scer_PCAs.pdf", width=10)

ggplot(Scer.galGonAll) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

####This is the best###
ggplot(Scer.galGonAll, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=pop), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF[-which(popLabDF$pop=="MedOak"),], aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = 0, alpha=.8) + 
  ggtitle("Scer PCA: All strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)
ggsave("ScerAll_PC1-PC2.pdf", height=8, width=12)

#ggplot(Scer.galGonAll) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
#  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
#  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)+ scale_x_continuous(limits=c(40,NA))

ggplot(Scer.galGonAll)+geom_polygon(data=polyDF, aes(x=PC2, y=PC3,group=pop), color="black", fill="grey50", alpha=0.8) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.galGonAll) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.galGonAll) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

Wine <- Scer.galGonAll[which(Scer.galGonAll$PC1>70),]
#ggplot(Wine) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
#  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
#  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)
Beer2 <- Scer.galGonAll[which(Scer.galGonAll$PC1<78 & Scer.galGonAll$PC1>40 & Scer.galGonAll$PC2>(-30)),]
#ggplot(Beer2) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
#  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
#  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)
write.table(Scer.galGonAll, file="Scer_PCAwLabels.txt", quote=F, sep="\t")

### LF ###
Scer.galGonLF <- read.table("PCAscores_Scer_galGonLF.txt", check.names = F, sep="\t")
Scer.galGonLF.strains <-row.names(Scer.galGonLF)
split <- strsplit(Scer.galGonLF.strains, "_")
Scer.galGonLF$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.galGonLF$syn <- unlist(lapply(split, function(x) x[2]))
Scer.galGonLF$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.galGonLF$tempPop <- unlist(lapply(split, function(x) x[4]))
Scer.galGonLF$phase <- unlist(lapply(split, function(x) x[5]))
Scer.galGonLF$phase[which(Scer.galGonLF$pop==1)] <- 1
Scer.galGonLF$phase[which(Scer.galGonLF$pop==2)] <- 2
Scer.galGonLF$tempPop[which(Scer.galGonLF$tempPop==1)] <- "Unplaced"
Scer.galGonLF$tempPop[which(Scer.galGonLF$tempPop==2)] <- "Unplaced"
Scer.galGonLF$origin <- NA
Scer.galGonLF$pop <- NA
for (strainID in Scer.galGonLF.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  strainPop <- as.character(renameOri$pop[which(renameOri$rename==strainID)])
  Scer.galGonLF$origin[which(rownames(Scer.galGonLF)==strainID)] <- strainOri
  Scer.galGonLF$pop[which(rownames(Scer.galGonLF)==strainID)] <- strainPop
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="GalGonLF")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="GalGonLF")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="GalGonLF")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="GalGonLF")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="GalGonLF")]*100, 2), "%", sep="")

uniPops <- unique(Scer.galGonLF$pop)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- Scer.galGonLF[which(Scer.galGonLF$pop==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF <- rbind(popLabDF, labTemp)
}

ggplot(Scer.galGonLF) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.galGonLF, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop), color="black", fill="grey50", alpha=0.7) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = 0, alpha=.7) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

polyDF2 <- data.frame()
popLabDF2 <- data.frame()
for (popName in uniPops) {
  popData <- Scer.galGonLF[which(Scer.galGonLF$pop==popName),]
  needed <- chull(x=popData$PC2, y=popData$PC3)
  wantedPopData <- popData[needed,]
  polyDF2 <- rbind(polyDF2, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF2 <- rbind(popLabDF2, labTemp)
}

ggplot(Scer.galGonLF) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: LF strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.galGonLF, aes(x=PC2, y=PC3)) + 
  geom_polygon(data=polyDF2, aes(x=PC2, y=PC3,group=pop), color="black", fill="grey50", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF2, aes(label = pop, x=PC2, y=PC3), vjust = 0, nudge_y = 0, alpha=.8) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)


polyDF3 <- data.frame()
popLabDF3 <- data.frame()
for (popName in uniPops) {
  popData <- Scer.galGonLF[which(Scer.galGonLF$pop==popName),]
  needed <- chull(x=popData$PC3, y=popData$PC4)
  wantedPopData <- popData[needed,]
  polyDF3 <- rbind(polyDF3, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC3=labPC3, PC4=labPC4, PC4=labPC4, PC5=labPC5)
  popLabDF3 <- rbind(popLabDF3, labTemp)
}

ggplot(Scer.galGonLF) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: LF strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.galGonLF, aes(x=PC3, y=PC4)) + 
  geom_polygon(data=polyDF3, aes(x=PC3, y=PC4,group=pop), color="black", fill="grey50", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF3, aes(label = pop, x=PC3, y=PC4), vjust = 0, nudge_y = 0, alpha=.8) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

polyDF4 <- data.frame()
popLabDF4 <- data.frame()
for (popName in uniPops) {
  popData <- Scer.galGonLF[which(Scer.galGonLF$pop==popName),]
  needed <- chull(x=popData$PC4, y=popData$PC5)
  wantedPopData <- popData[needed,]
  polyDF4 <- rbind(polyDF4, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC4=labPC4, PC5=labPC5, PC5=labPC5, PC5=labPC5)
  popLabDF4 <- rbind(popLabDF4, labTemp)
}

ggplot(Scer.galGonLF) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: LF strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

ggplot(Scer.galGonLF, aes(x=PC4, y=PC5)) + 
  geom_polygon(data=polyDF4, aes(x=PC4, y=PC5,group=pop), color="black", fill="grey50", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF4, aes(label = pop, x=PC4, y=PC5), vjust = 0, nudge_y = 0, alpha=.8) + 
  ggtitle("Scer PCA: All strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

sakeAsian <- Scer.galGonLF[which(Scer.galGonAll$PC1<(-100)),]
ggplot(sakeAsian) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Scer PCA: LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)
write.table(sakeAsian, "sakeAsian_PCA.txt", quote=F, sep="\t")

phiWAetc <- Scer.galGonLF[which(Scer.galGonAll$PC1>(-100) & Scer.galGonAll$PC1<0),]
phiWAetc <- phiWAetc[which(phiWAetc$PC2<(-50)),]
ggplot(phiWAetc) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Scer PCA: LF strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)
write.table(phiWAetc, "phiWA-NA-JPN_PCA.txt", quote=F, sep="\t")


## GalGonNoAsiaWAlF ###
Scer.GalGonNoAsiaWAlF <- read.table("PCAscores_Scer_GalGonNoAsiaWAlF.txt", check.names = F, sep="\t")
Scer.GalGonNoAsiaWAlF.strains <-row.names(Scer.GalGonNoAsiaWAlF)
split <- strsplit(Scer.GalGonNoAsiaWAlF.strains, "_")
Scer.GalGonNoAsiaWAlF$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.GalGonNoAsiaWAlF$syn <- unlist(lapply(split, function(x) x[2]))
Scer.GalGonNoAsiaWAlF$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.GalGonNoAsiaWAlF$pop <- unlist(lapply(split, function(x) x[4]))
Scer.GalGonNoAsiaWAlF$phase <- unlist(lapply(split, function(x) x[5]))
Scer.GalGonNoAsiaWAlF$phase[which(Scer.GalGonNoAsiaWAlF$pop==1)] <- 1
Scer.GalGonNoAsiaWAlF$phase[which(Scer.GalGonNoAsiaWAlF$pop==2)] <- 2
Scer.GalGonNoAsiaWAlF$pop[which(Scer.GalGonNoAsiaWAlF$pop==1)] <- "Unplaced"
Scer.GalGonNoAsiaWAlF$pop[which(Scer.GalGonNoAsiaWAlF$pop==2)] <- "Unplaced"
Scer.GalGonNoAsiaWAlF$origin <- NA
for (strainID in Scer.GalGonNoAsiaWAlF.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.GalGonNoAsiaWAlF$origin[which(rownames(Scer.GalGonNoAsiaWAlF)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="GalGonNoAsiaWAlF")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="GalGonNoAsiaWAlF")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="GalGonNoAsiaWAlF")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="GalGonNoAsiaWAlF")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="GalGonNoAsiaWAlF")]*100, 2), "%", sep="")

ggplot(Scer.GalGonNoAsiaWAlF) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) +
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) +
  ggtitle("Scer PCA: Dom strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.GalGonNoAsiaWAlF) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) +
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) +
  ggtitle("Scer PCA: Dom strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.GalGonNoAsiaWAlF) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) +
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) +
  ggtitle("Scer PCA: Dom strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.GalGonNoAsiaWAlF) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) +
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) +
  ggtitle("Scer PCA: Dom strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)


#### PHI-Sake ###
Scer.PHI_Sake <- read.table("PCAscores_Scer_PHI_Sake.txt", check.names = F, sep="\t")
Scer.PHI_Sake.strains <-row.names(Scer.PHI_Sake)
split <- strsplit(Scer.PHI_Sake.strains, "_")
Scer.PHI_Sake$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.PHI_Sake$syn <- unlist(lapply(split, function(x) x[2]))
Scer.PHI_Sake$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.PHI_Sake$pop <- unlist(lapply(split, function(x) x[4]))
Scer.PHI_Sake$phase <- unlist(lapply(split, function(x) x[5]))
Scer.PHI_Sake$phase[which(Scer.PHI_Sake$pop==1)] <- 1
Scer.PHI_Sake$phase[which(Scer.PHI_Sake$pop==2)] <- 2
Scer.PHI_Sake$pop[which(Scer.PHI_Sake$pop==1)] <- "Unplaced"
Scer.PHI_Sake$pop[which(Scer.PHI_Sake$pop==2)] <- "Unplaced"
Scer.PHI_Sake$origin <- NA
for (strainID in Scer.PHI_Sake.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.PHI_Sake$origin[which(rownames(Scer.PHI_Sake)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="PHI_Sake")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="PHI_Sake")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="PHI_Sake")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="PHI_Sake")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="PHI_Sake")]*100, 2), "%", sep="")

ggplot(Scer.PHI_Sake) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: PHI_Sake strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.PHI_Sake) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: PHI_Sake strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)


#### BreadMixed ###
Scer.BreadMixed <- read.table("PCAscores_Scer_BreadMixed.txt", check.names = F, sep="\t")
Scer.BreadMixed.strains <-row.names(Scer.BreadMixed)
split <- strsplit(Scer.BreadMixed.strains, "_")
Scer.BreadMixed$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.BreadMixed$syn <- unlist(lapply(split, function(x) x[2]))
Scer.BreadMixed$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.BreadMixed$pop <- unlist(lapply(split, function(x) x[4]))
Scer.BreadMixed$phase <- unlist(lapply(split, function(x) x[5]))
Scer.BreadMixed$phase[which(Scer.BreadMixed$pop==1)] <- 1
Scer.BreadMixed$phase[which(Scer.BreadMixed$pop==2)] <- 2
Scer.BreadMixed$pop[which(Scer.BreadMixed$pop==1)] <- "Unplaced"
Scer.BreadMixed$pop[which(Scer.BreadMixed$pop==2)] <- "Unplaced"
Scer.BreadMixed$origin <- NA
for (strainID in Scer.BreadMixed.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.BreadMixed$origin[which(rownames(Scer.BreadMixed)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="BreadMixed")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="BreadMixed")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="BreadMixed")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="BreadMixed")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="BreadMixed")]*100, 2), "%", sep="")

ggplot(Scer.BreadMixed) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: BreadMixed strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.BreadMixed) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: BreadMixed strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.BreadMixed) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: BreadMixed strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

#### Wine ###
Scer.Wine <- read.table("PCAscores_Scer_Wine.txt", check.names = F, sep="\t")
Scer.Wine.strains <-row.names(Scer.Wine)
split <- strsplit(Scer.Wine.strains, "_")
Scer.Wine$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Wine$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Wine$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Wine$pop <- unlist(lapply(split, function(x) x[4]))
Scer.Wine$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Wine$phase[which(Scer.Wine$pop==1)] <- 1
Scer.Wine$phase[which(Scer.Wine$pop==2)] <- 2
Scer.Wine$pop[which(Scer.Wine$pop==1)] <- "Unplaced"
Scer.Wine$pop[which(Scer.Wine$pop==2)] <- "Unplaced"
Scer.Wine$origin <- NA
for (strainID in Scer.Wine.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.Wine$origin[which(rownames(Scer.Wine)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Wine")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Wine")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Wine")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Wine")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Wine")]*100, 2), "%", sep="")

ggplot(Scer.Wine) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Wine strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.Wine) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Wine strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.Wine) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Wine strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

#### Wine_noWA ###
Scer.Wine_noWA <- read.table("PCAscores_Scer_WineNoWA.txt", check.names = F, sep="\t")
Scer.Wine_noWA.strains <-row.names(Scer.Wine_noWA)
split <- strsplit(Scer.Wine_noWA.strains, "_")
Scer.Wine_noWA$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Wine_noWA$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Wine_noWA$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Wine_noWA$pop <- unlist(lapply(split, function(x) x[4]))
Scer.Wine_noWA$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Wine_noWA$phase[which(Scer.Wine_noWA$pop==1)] <- 1
Scer.Wine_noWA$phase[which(Scer.Wine_noWA$pop==2)] <- 2
Scer.Wine_noWA$pop[which(Scer.Wine_noWA$pop==1)] <- "Unplaced"
Scer.Wine_noWA$pop[which(Scer.Wine_noWA$pop==2)] <- "Unplaced"
Scer.Wine_noWA$origin <- NA
for (strainID in Scer.Wine_noWA.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.Wine_noWA$origin[which(rownames(Scer.Wine_noWA)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Wine_noWA")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Wine_noWA")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Wine_noWA")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Wine_noWA")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Wine_noWA")]*100, 2), "%", sep="")

ggplot(Scer.Wine_noWA) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Wine_noWA strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.Wine_noWA) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Wine_noWA strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.Wine_noWA) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Wine_noWA strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.Wine_noWA) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Wine_noWA strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

#### Beer2 ###
Scer.Beer2 <- read.table("PCAscores_Scer_Beer2.txt", check.names = F, sep="\t")
Scer.Beer2.strains <-row.names(Scer.Beer2)
split <- strsplit(Scer.Beer2.strains, "_")
Scer.Beer2$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Beer2$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Beer2$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Beer2$pop <- unlist(lapply(split, function(x) x[4]))
Scer.Beer2$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Beer2$phase[which(Scer.Beer2$pop==1)] <- 1
Scer.Beer2$phase[which(Scer.Beer2$pop==2)] <- 2
Scer.Beer2$pop[which(Scer.Beer2$pop==1)] <- "Unplaced"
Scer.Beer2$pop[which(Scer.Beer2$pop==2)] <- "Unplaced"
Scer.Beer2$origin <- NA
for (strainID in Scer.Beer2.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.Beer2$origin[which(rownames(Scer.Beer2)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Beer2")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Beer2")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Beer2")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Beer2")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Beer2")]*100, 2), "%", sep="")

ggplot(Scer.Beer2) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2 strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.Beer2) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2 strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.Beer2) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2 strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.Beer2) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2 strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

#### Beer2Wine ###
Scer.Beer2Wine <- read.table("PCAscores_Scer_Beer2Wine.txt", check.names = F, sep="\t")
Scer.Beer2Wine.strains <-row.names(Scer.Beer2Wine)
split <- strsplit(Scer.Beer2Wine.strains, "_")
Scer.Beer2Wine$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Beer2Wine$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Beer2Wine$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Beer2Wine$pop <- unlist(lapply(split, function(x) x[4]))
Scer.Beer2Wine$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Beer2Wine$phase[which(Scer.Beer2Wine$pop==1)] <- 1
Scer.Beer2Wine$phase[which(Scer.Beer2Wine$pop==2)] <- 2
Scer.Beer2Wine$pop[which(Scer.Beer2Wine$pop==1)] <- "Unplaced"
Scer.Beer2Wine$pop[which(Scer.Beer2Wine$pop==2)] <- "Unplaced"
Scer.Beer2Wine$origin <- NA
for (strainID in Scer.Beer2Wine.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.Beer2Wine$origin[which(rownames(Scer.Beer2Wine)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Beer2Wine")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Beer2Wine")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Beer2Wine")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Beer2Wine")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Beer2Wine")]*100, 2), "%", sep="")

ggplot(Scer.Beer2Wine) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.Beer2Wine) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.Beer2Wine) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.Beer2Wine) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

#### Beer2Wine.noWA ###
Scer.Beer2Wine.noWA <- read.table("PCAscores_Scer_Beer2Wine.noWA.txt", check.names = F, sep="\t")
Scer.Beer2Wine.noWA.strains <-row.names(Scer.Beer2Wine.noWA)
split <- strsplit(Scer.Beer2Wine.noWA.strains, "_")
Scer.Beer2Wine.noWA$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Beer2Wine.noWA$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Beer2Wine.noWA$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Beer2Wine.noWA$tempPop <- unlist(lapply(split, function(x) x[4]))
Scer.Beer2Wine.noWA$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Beer2Wine.noWA$phase[which(Scer.Beer2Wine.noWA$pop==1)] <- 1
Scer.Beer2Wine.noWA$phase[which(Scer.Beer2Wine.noWA$pop==2)] <- 2
Scer.Beer2Wine.noWA$tempPop[which(Scer.Beer2Wine.noWA$tempPop==1)] <- "Unplaced"
Scer.Beer2Wine.noWA$tempPop[which(Scer.Beer2Wine.noWA$tempPop==2)] <- "Unplaced"
Scer.Beer2Wine.noWA$origin <- NA
Scer.Beer2Wine.noWA$pop <- NA
for (strainID in Scer.Beer2Wine.noWA.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  strainPop <- as.character(renameOri$pop[which(renameOri$rename==strainID)])
  Scer.Beer2Wine.noWA$origin[which(rownames(Scer.Beer2Wine.noWA)==strainID)] <- strainOri
  Scer.Beer2Wine.noWA$pop[which(rownames(Scer.Beer2Wine.noWA)==strainID)] <- strainPop
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")

uniPops <- unique(Scer.Beer2Wine.noWA$pop)
#uniPops <- uniPops[-which(uniPops=="Mosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- Scer.Beer2Wine.noWA[which(Scer.Beer2Wine.noWA$pop==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF <- rbind(popLabDF, labTemp)
}

####This is the best###
ggplot(Scer.Beer2Wine.noWA, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=pop), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_y = 0, alpha=.8) + 
  ggtitle("Scer PCA: Beer2/Wine strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("ScerBeer2Wine_PC1-PC2.pdf", height=8, width=12)


ggplot(Scer.Beer2Wine.noWA) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine.noWA strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.Beer2Wine.noWA) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine.noWA strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.Beer2Wine.noWA) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine.noWA strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.Beer2Wine.noWA) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer2Wine.noWA strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)


#### Beer1 ###
Scer.Beer1 <- read.table("PCAscores_Scer_Beer1.txt", check.names = F, sep="\t")
Scer.Beer1.strains <-row.names(Scer.Beer1)
split <- strsplit(Scer.Beer1.strains, "_")
Scer.Beer1$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Beer1$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Beer1$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Beer1$pop <- unlist(lapply(split, function(x) x[4]))
Scer.Beer1$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Beer1$phase[which(Scer.Beer1$pop==1)] <- 1
Scer.Beer1$phase[which(Scer.Beer1$pop==2)] <- 2
Scer.Beer1$pop[which(Scer.Beer1$pop==1)] <- "Unplaced"
Scer.Beer1$pop[which(Scer.Beer1$pop==2)] <- "Unplaced"
Scer.Beer1$origin <- NA
for (strainID in Scer.Beer1.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.Beer1$origin[which(rownames(Scer.Beer1)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Beer1")]*100, 2), "%", sep="")

ggplot(Scer.Beer1) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1 strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.Beer1) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1 strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.Beer1) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1 strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.Beer1) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1 strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

#### Beer1_lessFro ###
renameOri <- read.table("beer1strainPopKey.txt", check.names = F, header = T)
popColors <- c("aleMosaic"="#666666", "belGerAltKol"="#909090", "US"="#AFAFAF", "britIsles"="#C9C9C9",  "lager-saaz" = "#DFDFDF", "lager-frohberg"="#F2F2F2")
popFill <- scale_fill_manual(values=c(colors, popColors), name="Origin", breaks = c("Wild", "Domesticated", "Hybrid"), labels=c("Wild", "Domesticated", "Hybrid"))

Scer.Beer1_lessFro <- read.table("PCAscores_Scer_Beer1lF.txt", check.names = F, sep="\t")
Scer.Beer1_lessFro.strains <-row.names(Scer.Beer1_lessFro)
split <- strsplit(Scer.Beer1_lessFro.strains, "_")
Scer.Beer1_lessFro$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.Beer1_lessFro$syn <- unlist(lapply(split, function(x) x[2]))
Scer.Beer1_lessFro$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.Beer1_lessFro$tempPop <- unlist(lapply(split, function(x) x[4]))
Scer.Beer1_lessFro$phase <- unlist(lapply(split, function(x) x[5]))
Scer.Beer1_lessFro$phase[which(Scer.Beer1_lessFro$pop==1)] <- 1
Scer.Beer1_lessFro$phase[which(Scer.Beer1_lessFro$pop==2)] <- 2
Scer.Beer1_lessFro$tempPop[which(Scer.Beer1_lessFro$tempPop==1)] <- "Unplaced"
Scer.Beer1_lessFro$tempPop[which(Scer.Beer1_lessFro$tempPop==2)] <- "Unplaced"
Scer.Beer1_lessFro$origin <- NA
Scer.Beer1_lessFro$pop <- NA
for (strainID in Scer.Beer1_lessFro.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  strainPop <- as.character(renameOri$pop[which(renameOri$rename==strainID)])
  Scer.Beer1_lessFro$origin[which(rownames(Scer.Beer1_lessFro)==strainID)] <- strainOri
  Scer.Beer1_lessFro$pop[which(rownames(Scer.Beer1_lessFro)==strainID)] <- strainPop
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="Beer2Wine.noWA")]*100, 2), "%", sep="")

uniPops <- unique(Scer.Beer1_lessFro$pop)
uniPops <- uniPops[-which(uniPops=="aleMosaic")]
polyDF <- data.frame()
popLabDF <- data.frame()
for (popName in uniPops) {
  popData <- Scer.Beer1_lessFro[which(Scer.Beer1_lessFro$pop==popName),]
  needed <- chull(x=popData$PC1, y=popData$PC2)
  wantedPopData <- popData[needed,]
  polyDF <- rbind(polyDF, wantedPopData)
  labPC1 <- min(wantedPopData$PC1)+((max(wantedPopData$PC1)-min(wantedPopData$PC1))/2)
  labPC2 <- min(wantedPopData$PC2)+((max(wantedPopData$PC2)-min(wantedPopData$PC2))/2)
  labPC3 <- min(wantedPopData$PC3)+((max(wantedPopData$PC3)-min(wantedPopData$PC3))/2)
  labPC4 <- min(wantedPopData$PC4)+((max(wantedPopData$PC4)-min(wantedPopData$PC4))/2)
  labPC5 <- min(wantedPopData$PC5)+((max(wantedPopData$PC5)-min(wantedPopData$PC5))/2)
  labTemp <- data.frame(pop = popName, PC1=labPC1, PC2=labPC2, PC3=labPC3, PC4=labPC4, PC5=labPC5)
  popLabDF <- rbind(popLabDF, labTemp)
}

####This is the best###
ggplot(Scer.Beer1_lessFro, aes(x=PC1, y=PC2)) + 
  geom_polygon(data=polyDF, aes(x=PC1, y=PC2,group=pop, fill=pop), color="black", alpha=0.8) + 
  geom_point(aes(fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_label(data=popLabDF, aes(label = pop, x=PC1, y=PC2), vjust = 0, nudge_x = -7, nudge_y = 2, alpha=.8) + 
  ggtitle("Scer PCA: Ale/Beer1 w/ less Frohberg strains") + popFill + theme_bw() + labs(x = PC1var, y=PC2var)

ggsave("ScerBeer1lessFro_PC1-PC2.pdf", height=8, width=12)

ggplot(Scer.Beer1_lessFro) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.Beer1_lessFro) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.Beer1_lessFro) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.Beer1_lessFro) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

lagerSubset <- Scer.Beer1_lessFro[which(Scer.Beer1_lessFro$PC1<29 & Scer.Beer1_lessFro$PC1>5 & Scer.Beer1_lessFro$PC2<(-6) & Scer.Beer1_lessFro$PC2>(-40)),]
ggplot(lagerSubset) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro lager subset") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)
ggplot(lagerSubset) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro lagerSubset") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)
ggplot(lagerSubset) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro lagerSubset") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)
ggplot(lagerSubset) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = 0, check_overlap = T) + 
  ggtitle("Scer PCA: Beer1_lessFro lagerSubset") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)

#### beerUS ###
Scer.beerUS <- read.table("PCAscores_Scer_US.txt", check.names = F, sep="\t")
Scer.beerUS.strains <-row.names(Scer.beerUS)
split <- strsplit(Scer.beerUS.strains, "_")
Scer.beerUS$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerUS$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerUS$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerUS$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerUS$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerUS$phase[which(Scer.beerUS$pop==1)] <- 1
Scer.beerUS$phase[which(Scer.beerUS$pop==2)] <- 2
Scer.beerUS$pop[which(Scer.beerUS$pop==1)] <- "Unplaced"
Scer.beerUS$pop[which(Scer.beerUS$pop==2)] <- "Unplaced"
Scer.beerUS$origin <- NA
for (strainID in Scer.beerUS.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerUS$origin[which(rownames(Scer.beerUS)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="beerUS")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="beerUS")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="beerUS")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="beerUS")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="beerUS")]*100, 2), "%", sep="")

ggplot(Scer.beerUS) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerUS) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerUS) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

ggplot(Scer.beerUS) + geom_point(aes(x=PC4, y=PC5, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC4, y=PC5, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS strains") + fillScale + theme_bw() + labs(x = PC4var, y=PC5var)


#### beerUS.Mosaic ###
Scer.beerUS.Mosaic <- read.table("PCAscores_Scer_US.Mosaic.txt", check.names = F, sep="\t")
Scer.beerUS.Mosaic.strains <-row.names(Scer.beerUS.Mosaic)
split <- strsplit(Scer.beerUS.Mosaic.strains, "_")
Scer.beerUS.Mosaic$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerUS.Mosaic$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerUS.Mosaic$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerUS.Mosaic$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerUS.Mosaic$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerUS.Mosaic$phase[which(Scer.beerUS.Mosaic$pop==1)] <- 1
Scer.beerUS.Mosaic$phase[which(Scer.beerUS.Mosaic$pop==2)] <- 2
Scer.beerUS.Mosaic$pop[which(Scer.beerUS.Mosaic$pop==1)] <- "Unplaced"
Scer.beerUS.Mosaic$pop[which(Scer.beerUS.Mosaic$pop==2)] <- "Unplaced"
Scer.beerUS.Mosaic$origin <- NA
for (strainID in Scer.beerUS.Mosaic.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerUS.Mosaic$origin[which(rownames(Scer.beerUS.Mosaic)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="beerUS.Mosaic")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="beerUS.Mosaic")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="beerUS.Mosaic")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="beerUS.Mosaic")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="beerUS.Mosaic")]*100, 2), "%", sep="")

ggplot(Scer.beerUS.Mosaic) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS.Mosaic strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerUS.Mosaic) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS.Mosaic strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerUS.Mosaic) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS.Mosaic strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

#### beerUS.MosaicLager ###
Scer.beerUS.MosaicLager <- read.table("PCAscores_Scer_US.MosaicLager.txt", check.names = F, sep="\t")
Scer.beerUS.MosaicLager.strains <-row.names(Scer.beerUS.MosaicLager)
split <- strsplit(Scer.beerUS.MosaicLager.strains, "_")
Scer.beerUS.MosaicLager$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerUS.MosaicLager$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerUS.MosaicLager$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerUS.MosaicLager$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerUS.MosaicLager$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerUS.MosaicLager$phase[which(Scer.beerUS.MosaicLager$pop==1)] <- 1
Scer.beerUS.MosaicLager$phase[which(Scer.beerUS.MosaicLager$pop==2)] <- 2
Scer.beerUS.MosaicLager$pop[which(Scer.beerUS.MosaicLager$pop==1)] <- "Unplaced"
Scer.beerUS.MosaicLager$pop[which(Scer.beerUS.MosaicLager$pop==2)] <- "Unplaced"
Scer.beerUS.MosaicLager$origin <- NA
for (strainID in Scer.beerUS.MosaicLager.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerUS.MosaicLager$origin[which(rownames(Scer.beerUS.MosaicLager)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="beerUS.MosaicLager")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="beerUS.MosaicLager")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="beerUS.MosaicLager")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="beerUS.MosaicLager")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="beerUS.MosaicLager")]*100, 2), "%", sep="")

ggplot(Scer.beerUS.MosaicLager) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) +
  ggtitle("Scer PCA: beerUS.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerUS.MosaicLager) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerUS.MosaicLager) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerUS.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)


#### beerBelGerAltKol ###
Scer.beerBelGerAltKol <- read.table("PCAscores_Scer_BelGerALtKol.txt", check.names = F, sep="\t")
Scer.beerBelGerAltKol.strains <-row.names(Scer.beerBelGerAltKol)
split <- strsplit(Scer.beerBelGerAltKol.strains, "_")
Scer.beerBelGerAltKol$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerBelGerAltKol$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerBelGerAltKol$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerBelGerAltKol$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerBelGerAltKol$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerBelGerAltKol$phase[which(Scer.beerBelGerAltKol$pop==1)] <- 1
Scer.beerBelGerAltKol$phase[which(Scer.beerBelGerAltKol$pop==2)] <- 2
Scer.beerBelGerAltKol$pop[which(Scer.beerBelGerAltKol$pop==1)] <- "Unplaced"
Scer.beerBelGerAltKol$pop[which(Scer.beerBelGerAltKol$pop==2)] <- "Unplaced"
Scer.beerBelGerAltKol$origin <- NA
for (strainID in Scer.beerBelGerAltKol.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerBelGerAltKol$origin[which(rownames(Scer.beerBelGerAltKol)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="BelGerAltKol")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="BelGerAltKol")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="BelGerAltKol")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="BelGerAltKol")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="BelGerAltKol")]*100, 2), "%", sep="")

ggplot(Scer.beerBelGerAltKol) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerBelGerAltKol) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerBelGerAltKol) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

#### beerBelGerAltKol.Mosaic ###
Scer.beerBelGerAltKol.Mosaic <- read.table("PCAscores_Scer_BelGerAltKol.Mosaic.txt", check.names = F, sep="\t")
Scer.beerBelGerAltKol.Mosaic.strains <-row.names(Scer.beerBelGerAltKol.Mosaic)
split <- strsplit(Scer.beerBelGerAltKol.Mosaic.strains, "_")
Scer.beerBelGerAltKol.Mosaic$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerBelGerAltKol.Mosaic$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerBelGerAltKol.Mosaic$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerBelGerAltKol.Mosaic$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerBelGerAltKol.Mosaic$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerBelGerAltKol.Mosaic$phase[which(Scer.beerBelGerAltKol.Mosaic$pop==1)] <- 1
Scer.beerBelGerAltKol.Mosaic$phase[which(Scer.beerBelGerAltKol.Mosaic$pop==2)] <- 2
Scer.beerBelGerAltKol.Mosaic$pop[which(Scer.beerBelGerAltKol.Mosaic$pop==1)] <- "Unplaced"
Scer.beerBelGerAltKol.Mosaic$pop[which(Scer.beerBelGerAltKol.Mosaic$pop==2)] <- "Unplaced"
Scer.beerBelGerAltKol.Mosaic$origin <- NA
for (strainID in Scer.beerBelGerAltKol.Mosaic.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerBelGerAltKol.Mosaic$origin[which(rownames(Scer.beerBelGerAltKol.Mosaic)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="BelGerAltKol.Mosaic")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="BelGerAltKol.Mosaic")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="BelGerAltKol.Mosaic")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="BelGerAltKol.Mosaic")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="beerBelGerAltKol.Mosaic")]*100, 2), "%", sep="")

ggplot(Scer.beerBelGerAltKol.Mosaic) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol.Mosaic strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerBelGerAltKol.Mosaic) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol.Mosaic strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerBelGerAltKol.Mosaic) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol.Mosaic strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

#### beerBelGerAltKol.MosaicLager ###
Scer.beerBelGerAltKol.MosaicLager <- read.table("PCAscores_Scer_BelGerAltKol.MosaicLager.txt", check.names = F, sep="\t")
Scer.beerBelGerAltKol.MosaicLager.strains <-row.names(Scer.beerBelGerAltKol.MosaicLager)
split <- strsplit(Scer.beerBelGerAltKol.MosaicLager.strains, "_")
Scer.beerBelGerAltKol.MosaicLager$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerBelGerAltKol.MosaicLager$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerBelGerAltKol.MosaicLager$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerBelGerAltKol.MosaicLager$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerBelGerAltKol.MosaicLager$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerBelGerAltKol.MosaicLager$phase[which(Scer.beerBelGerAltKol.MosaicLager$pop==1)] <- 1
Scer.beerBelGerAltKol.MosaicLager$phase[which(Scer.beerBelGerAltKol.MosaicLager$pop==2)] <- 2
Scer.beerBelGerAltKol.MosaicLager$pop[which(Scer.beerBelGerAltKol.MosaicLager$pop==1)] <- "Unplaced"
Scer.beerBelGerAltKol.MosaicLager$pop[which(Scer.beerBelGerAltKol.MosaicLager$pop==2)] <- "Unplaced"
Scer.beerBelGerAltKol.MosaicLager$origin <- NA
for (strainID in Scer.beerBelGerAltKol.MosaicLager.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerBelGerAltKol.MosaicLager$origin[which(rownames(Scer.beerBelGerAltKol.MosaicLager)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="BelGerAltKol.MosaicLager")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="BelGerAltKol.MosaicLager")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="BelGerAltKol.MosaicLager")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="BelGerAltKol.MosaicLager")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="BelGerAltKol.MosaicLager")]*100, 2), "%", sep="")

ggplot(Scer.beerBelGerAltKol.MosaicLager) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerBelGerAltKol.MosaicLager) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerBelGerAltKol.MosaicLager) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBelGerAltKol.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)


#### beerBritishIsles ###
Scer.beerBritishIsles <- read.table("PCAscores_Scer_BritishIsles.txt", check.names = F, sep="\t")
Scer.beerBritishIsles.strains <-row.names(Scer.beerBritishIsles)
split <- strsplit(Scer.beerBritishIsles.strains, "_")
Scer.beerBritishIsles$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerBritishIsles$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerBritishIsles$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerBritishIsles$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerBritishIsles$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerBritishIsles$phase[which(Scer.beerBritishIsles$pop==1)] <- 1
Scer.beerBritishIsles$phase[which(Scer.beerBritishIsles$pop==2)] <- 2
Scer.beerBritishIsles$pop[which(Scer.beerBritishIsles$pop==1)] <- "Unplaced"
Scer.beerBritishIsles$pop[which(Scer.beerBritishIsles$pop==2)] <- "Unplaced"
Scer.beerBritishIsles$origin <- NA
for (strainID in Scer.beerBritishIsles.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerBritishIsles$origin[which(rownames(Scer.beerBritishIsles)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="BritishIsles")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="BritishIsles")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="BritishIsles")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="BritishIsles")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="BritishIsles")]*100, 2), "%", sep="")

ggplot(Scer.beerBritishIsles) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerBritishIsles) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerBritishIsles) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

#### beerBritishIsles.Mosaic ###
Scer.beerBritishIsles.Mosaic <- read.table("PCAscores_Scer_BritishIsles.Mosaic.txt", check.names = F, sep="\t")
Scer.beerBritishIsles.Mosaic.strains <-row.names(Scer.beerBritishIsles.Mosaic)
split <- strsplit(Scer.beerBritishIsles.Mosaic.strains, "_")
Scer.beerBritishIsles.Mosaic$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerBritishIsles.Mosaic$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerBritishIsles.Mosaic$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerBritishIsles.Mosaic$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerBritishIsles.Mosaic$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerBritishIsles.Mosaic$phase[which(Scer.beerBritishIsles.Mosaic$pop==1)] <- 1
Scer.beerBritishIsles.Mosaic$phase[which(Scer.beerBritishIsles.Mosaic$pop==2)] <- 2
Scer.beerBritishIsles.Mosaic$pop[which(Scer.beerBritishIsles.Mosaic$pop==1)] <- "Unplaced"
Scer.beerBritishIsles.Mosaic$pop[which(Scer.beerBritishIsles.Mosaic$pop==2)] <- "Unplaced"
Scer.beerBritishIsles.Mosaic$origin <- NA
for (strainID in Scer.beerBritishIsles.Mosaic.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerBritishIsles.Mosaic$origin[which(rownames(Scer.beerBritishIsles.Mosaic)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="BritishIsles.Mosaic")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="BritishIsles.Mosaic")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="BritishIsles.Mosaic")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="BritishIsles.Mosaic")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="BritishIsles.Mosaic")]*100, 2), "%", sep="")

ggplot(Scer.beerBritishIsles.Mosaic) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles.Mosaic strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerBritishIsles.Mosaic) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles.Mosaic strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerBritishIsles.Mosaic) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles.Mosaic strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

#### beerBritishIsles.MosaicLager ###
Scer.beerBritishIsles.MosaicLager <- read.table("PCAscores_Scer_BritishIsles.MosaicLager.txt", check.names = F, sep="\t")
Scer.beerBritishIsles.MosaicLager.strains <-row.names(Scer.beerBritishIsles.MosaicLager)
split <- strsplit(Scer.beerBritishIsles.MosaicLager.strains, "_")
Scer.beerBritishIsles.MosaicLager$strainID <- unlist(lapply(split, function(x) x[1]))
Scer.beerBritishIsles.MosaicLager$syn <- unlist(lapply(split, function(x) x[2]))
Scer.beerBritishIsles.MosaicLager$spIso <- unlist(lapply(split, function(x) x[3]))
Scer.beerBritishIsles.MosaicLager$pop <- unlist(lapply(split, function(x) x[4]))
Scer.beerBritishIsles.MosaicLager$phase <- unlist(lapply(split, function(x) x[5]))
Scer.beerBritishIsles.MosaicLager$phase[which(Scer.beerBritishIsles.MosaicLager$pop==1)] <- 1
Scer.beerBritishIsles.MosaicLager$phase[which(Scer.beerBritishIsles.MosaicLager$pop==2)] <- 2
Scer.beerBritishIsles.MosaicLager$pop[which(Scer.beerBritishIsles.MosaicLager$pop==1)] <- "Unplaced"
Scer.beerBritishIsles.MosaicLager$pop[which(Scer.beerBritishIsles.MosaicLager$pop==2)] <- "Unplaced"
Scer.beerBritishIsles.MosaicLager$origin <- NA
for (strainID in Scer.beerBritishIsles.MosaicLager.strains) {
  strainOri <- as.character(renameOri$origin[which(renameOri$rename==strainID)])
  Scer.beerBritishIsles.MosaicLager$origin[which(rownames(Scer.beerBritishIsles.MosaicLager)==strainID)] <- strainOri
}
PC1var <- paste("PC1 = ", round(ScerSum$PC1[which(ScerSum$analysis=="BritishIsles.MosaicLager")]*100, 2), "%", sep="")
PC2var <- paste("PC2 = ", round(ScerSum$PC2[which(ScerSum$analysis=="BritishIsles.MosaicLager")]*100, 2), "%", sep="")
PC3var <- paste("PC3 = ", round(ScerSum$PC3[which(ScerSum$analysis=="BritishIsles.MosaicLager")]*100, 2), "%", sep="")
PC4var <- paste("PC4 = ", round(ScerSum$PC4[which(ScerSum$analysis=="BritishIsles.MosaicLager")]*100, 2), "%", sep="")
PC5var <- paste("PC5 = ", round(ScerSum$PC5[which(ScerSum$analysis=="BritishIsles.MosaicLager")]*100, 2), "%", sep="")

ggplot(Scer.beerBritishIsles.MosaicLager) + geom_point(aes(x=PC1, y=PC2, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC1, y=PC2, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC1var, y=PC2var)

ggplot(Scer.beerBritishIsles.MosaicLager) + geom_point(aes(x=PC2, y=PC3, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC2, y=PC3, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC2var, y=PC3var)

ggplot(Scer.beerBritishIsles.MosaicLager) + geom_point(aes(x=PC3, y=PC4, fill=origin), shape=21, size = 5, alpha=0.8) + 
  geom_text(aes(x=PC3, y=PC4, label = paste(strainID,spIso,pop,sep="-")), vjust = 0, nudge_y = -2, check_overlap = T) + 
  ggtitle("Scer PCA: beerBritishIsles.MosaicLager strains") + fillScale + theme_bw() + labs(x = PC3var, y=PC4var)

dev.off()


