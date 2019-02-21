options(scipen=999)
options(stringsAsFactors = F)
library("ggplot2")

setwd("~/Documents/Hittinger_Lab/whole_genome/brewingProject/contributionPlots")
hybridEnv <- read.table("strainEnvHyb.txt", header=T)
uniEnv <- unique(hybridEnv$Environment)
uniHyb <- unique(hybridEnv$Hybrid)
uniSpecies <- c("Scer", "Skud", "Suva", "Seub")
speciesColors <- c("Scer"="red", "Spar"="#FFDB00", "Smik"="#49FF00", "Skud"="#00FF92", "SkudZP"="#00FF92", "SkudIFO"="#00FF92", "Sarb"="#0092FF", "Suva"="#4900FF", "Seub"="#FF00DB")
fillScale <- scale_fill_manual(values=speciesColors, name="Species", breaks=uniSpecies)

envDF <- data.frame()
for (isoEnv in uniEnv) {
  envData <- hybridEnv[which(hybridEnv$Environment==isoEnv),]
  envHybs <- unique(envData$Hybrid)
  for (hyb in envHybs) {
    subData <- envData[which(envData$Hybrid==hyb),]
    count <- nrow(subData)
    temp <- data.frame(Env=isoEnv, Hyb=hyb, count=count)
    envDF <- rbind(envDF, temp)
  }
}

envMod <- envDF
envMod$Env[which(envMod$Env=="DISTILLERY")] <- "UNKNOWN"

envColors <- c("BEER"="#FFE0B3", "CIDER"="#FFFF00", "DRINK"="#00FF4D", "FRUIT"="#00E5FF", "UNKNOWN"="#AAAAAA", "WILD"="#0019FF", "WINE"="#4C00FF")
#envFill <- scale_fill_manual(values=gray.colors(length(uniEnv), start = 0.0, end = .9), name="Environment")
envFill <- scale_fill_manual(values=alpha(envColors, 0.5), name="Environment", labels=c("BEER", "CIDER", "DRINK", "FRUIT", "UNKNOWN/\nOTHER", "WILD", "WINE"))
envScale <- scale_color_manual(values=envColors, guide=F)
ggplot(envMod, aes(x=Hyb))+geom_bar(aes(fill=Env, weight=count)) +envScale+envFill+theme_bw()+ labs(x="Hyrbid", y="Isolation Count")

#ggplot(envDF, aes(x=Hyb))+geom_bar(aes(fill=Env, weight=count, color=Env)) +envColor+envFill+theme_bw()+ labs(x="Hyrbid", y="Isolation Count")

cpEnvDF <- envMod
cpEnvDF$Hyb[which(cpEnvDF$Hyb=="lager-Frohberg")] <- "Lager\nFrohberg"
cpEnvDF$Hyb[which(cpEnvDF$Hyb=="lager-Saaz")] <- "Lager\nSaaz"
cpEnvDF$Hyb <- gsub("X", " x\n", cpEnvDF$Hyb)
ggplot(cpEnvDF, aes(x=Hyb))+geom_bar(aes(fill=Env, weight=count, color=Env)) +envScale+envFill+theme_bw()+theme(axis.text.x=element_text(face="italic"))+ labs(x="Hyrbid", y="Isolation Count")
ggsave("hybEnv_newColor_vert.pdf", height=8, width=6)

renameEnv <- envMod
renameEnv$Hyb[which(renameEnv$Hyb=="lager-Frohberg")] <- "Lager\nFrohberg"
renameEnv$Hyb[which(renameEnv$Hyb=="lager-Saaz")] <- "Lager\nSaaz"
renameEnv$Hyb[which(renameEnv$Hyb=="SuvaXSeub")] <- "Seub x\nSuva"
renameEnv$Hyb[which(renameEnv$Hyb=="ScerXSuvaXSeub")] <- "Scer x\nSeub x\nSuva"
renameEnv$Hyb[which(renameEnv$Hyb=="ScerXSkudXSeub")] <- "Scer x\nSeub x\nSkud"
renameEnv$Hyb[which(renameEnv$Hyb=="ScerXSkudXSuvaXSeub")] <- "Scer x\nSeub x\nSkud x\nSuva"
renameEnv$Hyb[which(renameEnv$Hyb=="ScerXSkud")] <- "Scer x\nSkud"
#cpEnvDF$Hyb <- gsub("X", " x\n", cpEnvDF$Hyb)
ggplot(renameEnv, aes(x=Hyb))+geom_bar(aes(fill=Env, weight=count, color=Env)) +envScale+envFill+theme_bw()+theme(axis.text.x=element_text(face="italic"))+ labs(x="Hyrbid", y="Isolation Count")
ggsave("hybEnv_newColor_vert_rename.pdf", height=8, width=6)


#ggplot(envDF, aes(x=Hyb))+geom_bar(aes(weight=count, color=Env), fill= "black") +envColor+envFill+theme_bw()+ labs(x="Hyrbid", y="Isolation Count")
#ggplot(envDF, aes(x=Hyb))+geom_bar(aes(fill=Env, weight=count), position="dodge")+envFill+theme_bw() + labs(x="Hyrbid", y="Isolation Count")
#ggplot(envDF, aes(x=Env))+geom_bar(aes(fill=Hyb, weight=count))

spEnv <- data.frame()
for (spName in uniSpecies) {
  spData <- hybridEnv[grep(spName,hybridEnv$StrainID),]
  if (spName == "Scer" | spName == "Seub") {
    lagerData <- hybridEnv[grep("lager",hybridEnv$StrainID),]
    spData <- rbind(spData, lagerData)
  }
  spUenv <- unique(spData$Environment)
  for (isoEnv  in spUenv) {
    subEnv <- spData[which(spData$Environment==isoEnv),]
    count <- nrow(subEnv)
    temp <- data.frame(Species=spName, Env=isoEnv, count=count)
    spEnv <- rbind(spEnv, temp)
  }
}

ggplot(spEnv, aes(x=Species))+geom_bar(aes(fill=Env, weight=count))
ggplot(spEnv, aes(x=Env))+geom_bar(aes(fill=Species, weight=count))+fillScale


pairEnv <- read.csv("Hybrid_ecology.csv")
pairEnv$pval <- 0
pairEnv$oddsRatio <- 0
for (i in (1:nrow(pairEnv))) {
  tempMatrix <- matrix(unlist(pairEnv[i,c(3,5,4,6)]), nrow=2)
  tempTest <- fisher.test(tempMatrix)
  pairEnv$pval[i] <- tempTest$p.value
  pairEnv$oddsRatio <- tempTest$estimate
}
pairEnv$Correction = p.adjust(pairEnv$pval, method = "bonferroni")
pairEnv$Sig = pairEnv$Correction <= 0.05


fruit <- pairEnv[which(pairEnv$Environment=="Fruit"),]
fruitChi <-  chisq.test(fruit$X1.1)

beer <- pairEnv[which(pairEnv$Environment=="Beer"),]
beerChi <- chisq.test(beer$X1.1)

wine <- pairEnv[which(pairEnv$Environment=="Wine"),]
wineChi <- chisq.test(wine$X1.1)

envChi <- data.frame(Env = c("Fruit", "Wine", "Beer"), p.value=c(fruitChi$p.value, wineChi$p.value, beerChi$p.value), chiSq = c(fruitChi$statistic, wineChi$statistic, beerChi$statistic))
envChi$Correction = p.adjust(envChi$p.value, method = "bonferroni")
envChi$Sig = envChi$Correction <= 0.05
###Fruit assoicated with SuvaXSeub, Beer assoicated with ScerXSeub, and Wine assoicated with ScerXSkud and SuvaXSeub but artifact of 0 for ScerXSeub


lager <- pairEnv[which(pairEnv$Hybrid=="ScerXSeub"),]
lagerChi <- chisq.test(lager$X1.1)

ScerXSkud <- pairEnv[which(pairEnv$Hybrid=="ScerXSkud"),]
ScerXSkudChi <- chisq.test(ScerXSkud$X1.1)

SuvaXSeub <- pairEnv[which(pairEnv$Hybrid=="SuvaXSeub"),]
SuvaXSeubChi <- chisq.test(SuvaXSeub$X1.1)

hybChi <- data.frame(Hyb = c("ScerXSeub", "ScerXSkud", "SuvaXSeub"), p.value=c(lagerChi$p.value, ScerXSkudChi$p.value, SuvaXSeubChi$p.value), chiSq = c(lagerChi$statistic, ScerXSkudChi$statistic, SuvaXSeubChi$statistic))
hybChi$Correction = p.adjust(hybChi$p.value, method = "bonferroni")
hybChi$Sig = hybChi$Correction <= 0.05
###ScerXSeub very assoicated with beer, SuvaXSeub no assoication, ScerXSkud assoicated with Wine and Beer but maybe an artifiact from a 0 for Fruit





countMatrix <- read.csv("Saccharomyces_ecology.csv")
countMatrix <- countMatrix[,1:6]

countMatrix$pval <- 0
countMatrix$oddsRatio <- 0
for (i in (1:nrow(countMatrix))) {
  tempMatrix <- matrix(unlist(countMatrix[i,c(3,5,4,6)]), nrow=2)
  tempTest <- fisher.test(tempMatrix)
  countMatrix$pval[i] <- tempTest$p.value
  countMatrix$oddsRatio <- tempTest$estimate
}

countMatrix$Correction = p.adjust(countMatrix$pval, method = "bonferroni")
countMatrix$Sig = countMatrix$Correction <= 0.05

##Scer and Seub positively assoicated with Beer. Driven by lagers. Try without lagers

Suva <- countMatrix[which(countMatrix$Hybrid=="Suva"),]
Suva$pval <- 0
Suva$oddsRatio <- 0
for (i in (1:nrow(Suva))) {
  tempMatrix <- matrix(unlist(Suva[i,c(3,5,4,6)]), nrow=2)
  tempTest <- fisher.test(tempMatrix)
  Suva$pval[i] <- tempTest$p.value
  Suva$oddsRatio <- tempTest$estimate
}
Suva$Correction = p.adjust(Suva$pval, method = "bonferroni")
Suva$Sig = Suva$Correction <= 0.05

Skud <- countMatrix[which(countMatrix$Hybrid=="Skud"),]
Skud$pval <- 0
Skud$oddsRatio <- 0
for (i in (1:nrow(Skud))) {
  tempMatrix <- matrix(unlist(Skud[i,c(3,5,4,6)]), nrow=2)
  tempTest <- fisher.test(tempMatrix)
  Skud$pval[i] <- tempTest$p.value
  Skud$oddsRatio <- tempTest$estimate
}
Skud$Correction = p.adjust(Skud$pval, method = "bonferroni")
Skud$Sig = Skud$Correction <= 0.05


SeubBeer <- matrix(unlist(countMatrix[1,c(3,5,4,6)]), nrow=2)
SeubBeerTest <- fisher.test(SeubBeer)



#Contains Galactose
Galactose = matrix(0,nrow =2, ncol = 2)
Galactose[1,1] = 5
Galactose[1,2] = 2
Galactose[2,1] = 0
Galactose[2,2] = 25
Galactose_fisher = fisher.test(Galactose)



+fillScale+facet_grid(. ~ mitoPa, scale='free_x',space="free_x")+ggtitle("Genome Content faceted by mito parent")+labs(x = "Strain", y="Percent Genomic Content")+theme_bw()+theme(axis.text.x=element_blank(),strip.text=element_text(face="italic"))+scale_y_continuous(limits = c(0,1))