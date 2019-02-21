options(scipen=999)
library("ggplot2")

setwd("~/Documents/Hittinger_Lab/whole_genome/brewingProject/contributionPlots")
#args <- commandArgs(TRUE)
ploidyFileName <- "ScerSkudSuvaSeub_ploidyInfo.txt"
chrLengthsFile <- "ScerSkudSuvaSeub_binnedChrLengths.txt"
mitoParentFile <- "mitoParentKey.txt"
strainOrderFile <- "strainOrder_contributionPlot.txt"

chrLengths <- read.table(chrLengthsFile, header=T)
mitoOrder <- read.table(mitoParentFile, header=T)
strainOrder <- read.table(strainOrderFile, header=T)

speciesColors <- c("Scer"="red", "Spar"="#FFDB00", "Smik"="#49FF00", "Skud"="#00FF92", "SkudZP"="#00FF92", "SkudIFO"="#00FF92", "Sarb"="#0092FF", "Suva"="#4900FF", "Seub"="#FF00DB")

ploidyTable <- read.table(ploidyFileName, header=T)
levels(ploidyTable$species) <- c("Scer", "Skud", "Suva", "Seub")
indicies <- which(ploidyTable$species=="SkudZP")
ploidyTable$species[indicies] <- "Skud"
maxPloidy <- max(ploidyTable$ploidy)
ploidyTable$yStart <- NA
ploidyTable$yEnd <- NA
#uniStrains <- unique(mitoOrder$strainID)
uniStrains <- unique(strainOrder$genomeName)
uniSpecies <- unique(ploidyTable$species)

ploidyTable$ploidyXregLen <- ploidyTable$ploidy * ploidyTable$regionLength
strainSummary <- data.frame()
for (strainName in uniStrains){
  mitoPar <- mitoOrder$mitoSp[which(mitoOrder$strainID==strainName)]
  strainTemp <- data.frame()
  strainData <- ploidyTable[which(ploidyTable$strainID==strainName),]
  total <- 0
  for (speciesName in uniSpecies) {
    spData <- strainData[which(strainData$species==speciesName),]
    spSum <- sum(spData$ploidyXregLen)
    spTemp <- data.frame(strainID=strainName, contVal = spSum, species=speciesName, mitoPa=mitoPar, contPlo = spSum/12000000)
    strainTemp <- rbind(strainTemp, spTemp)
    total <- total + spSum
  }
  strainTemp$contPer <- strainTemp$contVal/total
  strainSummary <- rbind(strainSummary, strainTemp)
}

strainSummary$species[which(strainSummary$species=="SkudZP")] <- "Skud"
levels(strainSummary$species) <- c("Scer", "Skud", "Suva", "Seub")
#levels(strainSummary$mitoPa) <- c("Scer", "Skud", "Suva", "Seub")

uniSpecies <- unique(strainSummary$mitoPa)
levels(uniSpecies) <- c("Scer", "Skud", "Suva", "Seub")
fillScale <- scale_fill_manual(values=speciesColors, name="Species")

ggplot(strainSummary, aes(x=strainID))+geom_bar(aes(fill=species, weight=contVal))+fillScale+facet_grid(. ~ mitoPa, scale='free_x',space="free_x")+labs(x = "Strain", y="Genomic basepair content")+theme_bw()+theme(axis.text.x=element_blank(),strip.text=element_text(face="italic"))

#pdf("mitoInheritancePlots.pdf", width=12)

ggplot(strainSummary, aes(x=strainID))+geom_bar(aes(fill=species, weight=contPer))+fillScale+facet_grid(. ~ factor(mitoPa, levels=c("Scer", "Skud", "Suva", "Seub")), scale='free_x',space="free_x")+labs(x = "Strain", y="Proportion Nuclear Genomic Content")+theme_bw()+theme(axis.text.x=element_blank(),strip.text=element_text(face="italic"), strip.background =element_rect(fill="white"), legend.text=element_text(face="italic"),plot.title = element_text(hjust = 0.5))+scale_y_continuous(limits = c(0,1))+ggtitle("Mitochondrial Genome")
ggsave("MitoInheritanceByGenomeContent_reordered.pdf", height=8, width=12)

ggplot(strainSummary, aes(x=strainID))+geom_bar(aes(fill=species, weight=contPlo))+fillScale+facet_grid( ~ mitoPa, scale='free_x',space="free_x")+ggtitle("Genome Content faceted by mito parent")+labs(x = "Strain", y="Genomic Ploidy")+theme_bw()+theme(axis.text.x=element_blank(),strip.text=element_text(face="italic"))+scale_y_continuous(limits = c(0,NA), breaks=c(0,1,2,3,4,5,6,7,8))

dev.off()

strainSummary$species <- as.character(strainSummary$species)
strainSummary$species[which(strainSummary$species=="SkudZP")] <- "Skud"

options(stringsAsFactors = F)

inherSum <- data.frame()
for (strainName in uniStrains){
  strainData <- strainSummary[which(strainSummary$strainID==strainName),]
  strainMito <- as.character(unique(strainData$mitoPa))
  maxCont <- max(strainData$contPer)
  maxSp <- as.character(strainData$species[which(strainData$contPer==maxCont)])
  spCont <- strainData$contPer[which(strainData$species==strainMito)]
  strainTemp <- data.frame(strainID=strainName, mitoPa=strainMito, maxSp = maxSp, mitoCont = spCont, maxCont = maxCont)
  inherSum <- rbind(inherSum, strainTemp)
}
inherSum$match <- inherSum$mitoPa==inherSum$maxSp
inherSum$match[which(inherSum$match=="FALSE")] <- 0
inherSum$match[which(inherSum$match=="TRUE")] <- 1

fit <- glm(match~maxCont, data=inherSum, family=binomial())
write.csv(capture.output(summary(fit)), "mitoInher_logisticRegression.csv")

ggplot(inherSum, aes(as.character(match), maxCont)) +geom_boxplot(color="grey40", outlier.shape = NA) + geom_jitter(aes(fill=mitoPa, color=mitoPa, pch=maxSp), alpha=.9, size=3)+scale_shape_manual(name = "Max Nulcear \nGenome", values=c("Scer"=21, "Suva"=22, "Seub"=23))+theme_bw()+scale_y_continuous("Largest Proprotion of Nuclear Genomic Content", limits=c(NA,1)) + scale_x_discrete(" Inherits Mitochondria of Species with Largest \nProportion of Nuclear Genomic Content", labels=c("False", "True"))+scale_fill_manual(values=speciesColors, name="Mito Genome", breaks=c("Scer", "Skud", "Suva", "Seub"))+scale_color_manual(values=speciesColors, name="Mito Genome", breaks=c("Scer", "Skud", "Suva", "Seub"))+theme(legend.text=element_text(face="italic"))
ggsave("MitoMatchToGenome.pdf", height=8, width=6)

+ggtitle("Largest Genomic Content Match to Mito")


ggplot(inherSum, aes(as.character(match)))+geom_bar()

library("aod")
wald.test(b=coef(fit), Sigma=vcov(fit), Terms=1)

dev.off()

ggplot(inherSum, aes(as.character(match), maxCont)) +geom_boxplot(color="grey40") + geom_jitter(aes(color=mitoPa), alpha=.9, size=3)+scale_color_manual(values=speciesColors, name="Mito Genome", breaks=uniSpecies)+theme_bw()+scale_y_continuous("Largest Genomic Contribution", limits=c(NA,1)) + scale_x_discrete(" Inherits Mitochondria of Largest Contributing Genome", labels=c("False", "True"))+ggtitle("Largest Genomic Content Match to Mito")+coord_flip()

quantile(inherSum$maxCont[which(inherSum$match==1)])
outliers <- inherSum[which(inherSum$match==0 & inherSum$maxCont>0.6532046),]

### Testing ###

ggplot(inherSum, aes(mitoPa, cont)) + geom_violin(aes(fill=mitoPa)) +fillScale
ggplot(inherSum, aes(match, maxCont)) + geom_point()

ggplot(inherSum) + geom_density(aes(maxCont, group=match, fill=match), alpha=.7, trim=T, adjust=.5)

1-pchisq(169.0 -114.1, 121 - 120)

inherSum$hybrid <- NA
for (i in 1:nrow(inherSum)) {
  strainName <- inherSum$strainID[i]
  splitLen <- length(strsplit(strainName, "X")[[1]])
  inherSum$hybrid[i] <- splitLen
}
inherSum$hybrid[which(inherSum$hybrid==1)] <- 2

ggplot(inherSum, aes(match, maxCont)) + geom_jitter(aes(color=hybrid))

match <- inherSum[which(inherSum$match==1),]
ggplot() + geom_density(data=match, aes(maxCont, fill=mitoPa), alpha=.7, trim=T, adjust=.5)+ geom_density(data=nonMatching, aes(maxCont, fill=mitoPa), alpha=.7, trim=T, adjust=.5)+fillScale

inherSum$resid <- inherSum$maxCont-inherSum$mitoCont
ggplot(inherSum) + geom_density(aes(maxCont, fill=mitoPa), alpha=.7, trim=T, adjust=.5)+facet_grid(match~.)+fillScale
ggplot(inherSum) + geom_density(aes(maxCont, fill=maxSp, color=mitoPa), alpha=.7, trim=T, adjust=.5)+facet_grid(match~.)+fillScale +scale_color_manual(values=speciesColors, name="Species", breaks=uniSpecies)


ggplot(inherSum) + geom_density(aes(mitoCont, fill=mitoPa), alpha=.7, trim=T, adjust=.5)+facet_grid(match~.)+fillScale


ggplot(inherSum) + geom_density(aes(mitoCont, fill=mitoPa), alpha=.7, trim=T, adjust=.5)+facet_grid(match~.)+fillScale
ggplot(inherSum) + geom_density(aes(resid, fill=mitoPa), alpha=.7, trim=T, adjust=.5)+facet_grid(match~.)+fillScale
ggplot(inherSum) + geom_density(aes(x=mitoCont, y=..count.., fill=mitoPa), alpha=.7, trim=T, adjust=.5)+fillScale


ggplot(inherSum) + geom_violin(aes(x=as.character(match), y=mitoCont, fill=mitoPa), alpha=.7)+fillScale
ggplot(inherSum) + geom_violin(aes(x=as.character(match), y=maxCont, fill=mitoPa), alpha=.7)+fillScale


fit2 <- glm(match~maxCont+hybrid, data=inherSum, family=binomial())
fit3 <- glm(match~maxCont+mitoPa, data=inherSum, family=binomial())

fit4 <- glm(match~mitoCont*mitoPa, data=inherSum, family=binomial())


skudEnv <- read.table("../SkudContEnv.txt", header=F, col.names = c("StrainID", "Env", "subPop"))
uniSkud <- unique(skudEnv$StrainID)
SkudCont <- inherSum[grep("Skud",inherSum$strainID),]
SkudCont$subPop<- NA
for (strainName in uniSkud) {
  subPop <- skudEnv$subPop[which(skudEnv$StrainID==strainName)]
  SkudCont$subPop[which(SkudCont$strainID==strainName)] <- subPop
}
ggplot(SkudCont) + geom_density(aes(maxCont, fill=subPop), alpha=.8) + facet_grid(mitoPa~.)
ggplot(SkudCont) + geom_density(aes(mitoCont, fill=subPop), alpha=.8) + facet_grid(mitoPa~.)


pairs <- inherSum[which(inherSum$hybrid==2 & inherSum$match==0),]
pairs$expect <- .5
pairs$ratio <- pairs$mitoCont/pairs$expect
pairsLow <- pairs[which(pairs$ratio<.9),]
summary(pairsLow)
pairs$resid <- pairs$maxCont-pairs$mitoCont

nonMatching <- inherSum[which(inherSum$match==0),]
nonMatching$resid <- nonMatching$maxCont - nonMatching$mitoCont
ggplot(nonMatching) + geom_density(aes(resid)) 

outlie <- nonMatching[which(nonMatching$resid>0.2),]


nonMatching <- inherSum[which(inherSum$match==0),]
nonMatching$hybrid <- NA
for (i in 1:nrow(nonMatching)) {
  strainName <- nonMatching$strainID[i]
  splitLen <- length(strsplit(strainName, "X")[[1]])
  nonMatching$hybrid[i] <- splitLen
}
nonMatching$hybrid[which(nonMatching$hybrid==1)] <- 2

nonMatching$expect <- NA
nonMatching$expect[which(nonMatching$hybrid==2)] <- .5
nonMatching$expect[which(nonMatching$hybrid==3)] <- .33333333
nonMatching$expect[which(nonMatching$hybrid==4)] <- .25

nonMatching$ratio <- nonMatching$mitoCont/nonMatching$expect

ggplot(nonMatching) + geom_histogram(aes(ratio))

nonMatchingLow <- nonMatching[which(nonMatching$ratio<.9),]


mean(inherSum$maxCont[which(inherSum$match==1)])
mean(inherSum$maxCont[which(inherSum$match==0)])

median(inherSum$maxCont[which(inherSum$match==1)])
median(inherSum$maxCont[which(inherSum$match==0)])



predict(fit)
