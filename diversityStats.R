options(scipen=999)
options(stringsAsFactors = F)
library(PopGenome)

reduced <- readData("reduced", include.unknown=F)
reducedUN <- readData("reduced", include.unknown=T)
reduced <- diversity.stats(reduced, pi=TRUE)
reducedUN <- diversity.stats(reducedUN, pi=TRUE)
nSitesRed <- reduced@n.sites
nSitesRedUN <- reducedUN@n.sites
validSitesRed <- reduced@n.valid.sites
validSitesRedUN <- reducedUN@n.valid.sites
piValuesRed <- reduced@Pi
piValuesRedUN <- reducedUN@Pi
genomCovRed <- validSitesRed/nSitesRed
genomCovRedUN <- validSitesRedUN/nSitesRedUN
adjustPiRed <- 100*(piValuesRed/validSitesRed)
adjustPiRedUN <- 100*(piValuesRedUN/validSitesRedUN)
redSummary <- data.frame(adjPi = adjustPiRed, rawPi = piValuesRed, validGenomePer=genomCovRed, validSites=validSitesRed, nSites=nSitesRed)
redUnSummary <- data.frame(adjPi = adjustPiRedUN, rawPi = piValuesRedUN, validGenomePer=genomCovRedUN, validSites=validSitesRedUN, nSites=nSitesRedUN)
totalSum <- rbind(redSummary, redUnSummary)
colnames(totalSum) <- c("adjPi", "rawPi", "validGenomePer", "validSites", "nSites")

popWGunknown <- readData("popWG", include.unknown=T)
popWGunknown <- diversity.stats(popWGunknown, pi=TRUE)
nSitesPopsUN <- popWGunknown@n.sites
validSitesPopsUN <- popWGunknown@n.valid.sites
piValuesPopsUN <- popWGunknown@Pi
genomCovUN <- popWGunknown@n.valid.sites/popWGunknown@n.sites
adjustPiUN <- 100*(popWGunknown@Pi/popWGunknown@n.valid.sites)
popSummaryUN <- data.frame(adjPi = adjustPiUN, rawPi = piValuesPopsUN, validGenomePer=genomCovUN, validSites=validSitesPopsUN, nSites=nSitesPopsUN)



britIsles <- read.big.fasta("britIsles.fasta", include.unknown=F)
allAle <- readData("aleWG", include.unknown=F, big.data=T, FAST=T)

####Scer####
wineHyb <- readData("wineHyb", include.unknown=FALSE)
beer2hyb <- readData("beer2hyb", include.unknown=FALSE)
beer1hyb <- readData("beer1hyb", include.unknown=FALSE)

chrNoHyb <- readData("chrNoHybs", include.unknown=FALSE)
chrWhyb <- readData("chrWhybs", include.unknown=FALSE)
WGnoHyb <- readData("WGnoHybs", include.unknown=FALSE)
WGwHyb <- readData("WGwHybs", include.unknown=FALSE)
ScerNoHyb <- readData("noHybs", include.unknown=FALSE)
ScerWhyb <- readData("wHybs", include.unknown=FALSE)
save.image(file="Scer_popGenome.RData")
noHybPops <- c("all", "PHI_Sake", "BreadMixed", "Wine", "Beer2", "Beer1", "beerUS", "beerBelGerAltKol", "beerBritishIsles")
hybPops <- c("allHyb", "WineHyb", "Beer2hyb", "Beer1hyb", "lager", "lagerSaaz", "lagerFrohberg")
all <- list(scan("Scer_All.txt", what="", sep="\n"))
allHyb <- list(scan("Scer_AllHyb.txt", what="", sep="\n"))
PHI.Sake <- list(scan("Scer_PHI_Sake.txt", what="", sep="\n"))
bread <- list(scan("Scer_BreadMixed.txt", what="", sep="\n"))
wine <- list(scan("Scer_Wine.txt", what="", sep="\n"))
wineHyb <- list(scan("Scer_WineHyb.txt", what="", sep="\n"))
beer2 <- list(scan("Scer_Beer2.txt", what="", sep="\n"))
beer2hyb <- list(scan("Scer_Beer2hyb.txt", what="", sep="\n"))
beer1 <- list(scan("Scer_Beer1.txt", what="", sep="\n"))
beer1hyb <- list(scan("Scer_Beer1hyb.txt", what="", sep="\n"))
beerUS <- list(scan("Scer_beerUS.txt", what="", sep="\n"))
beerBel <- list(scan("Scer_beerBelGerAltKol.txt", what="", sep="\n"))
beerBri <- list(scan("Scer_beerBritishIsles.txt", what="", sep="\n"))
lager <- list(scan("Scer_lager.txt", what="", sep="\n"))
lagerS <- list(scan("Scer_lagerSaaz.txt", what="", sep="\n"))
lagerF <- list(scan("Scer_lagerFrohberg.txt", what="", sep="\n"))
ScerNoHyb <- set.populations(ScerNoHyb, c(all, PHI.Sake, bread, wine, beer2, beer1, beerUS, beerBel, beerBri))
ScerWhyb <- set.populations(ScerWhyb, c(allHyb, wineHyb, beer2hyb, beer1hyb, lager, lagerS, lagerF))
chrNoHyb <- set.populations(chrNoHyb, c(all, PHI.Sake, bread, wine, beer2, beer1, beerUS, beerBel, beerBri))
chrWhyb <- set.populations(chrWhyb, c(allHyb, wineHyb, beer2hyb, beer1hyb, lager, lagerS, lagerF))
WGnoHyb <- set.populations(WGnoHyb, c(all, PHI.Sake, bread, wine, beer2, beer1, beerUS, beerBel, beerBri))
WGwHyb <- set.populations(WGwHyb, c(allHyb, wineHyb, beer2hyb, beer1hyb, lager, lagerS, lagerF))
ScerNoHyb <- diversity.stats(ScerNoHyb, pi=TRUE)
ScerWhyb <- diversity.stats(ScerWhyb, pi=TRUE)
nSitesNH <- ScerNoHyb@n.sites
validSitesNH <- ScerNoHyb@n.valid.sites
piValuesNH <- ScerNoHyb@Pi
nSitesHyb <- ScerWhyb@n.sites
validSitesHyb <- ScerWhyb@n.valid.sites
piValuesHyb <- ScerWhyb@Pi
outDF <- data.frame()
for (i in 1:length(noHybPops)) {
  temp <- data.frame(pop = noHybPops[i], pi = 100*(piValuesNH[i]/validSitesNH), validGenomePer=100*(validSitesNH/nSitesNH), validSites = validSitesNH)
  outDF <- rbind(outDF, temp)
}
for (j in 1:length(hybPops)){
  temp <- data.frame(pop = hybPops[j], pi = 100*(piValuesHyb[j]/validSitesHyb), validGenomePer=100*(validSitesHyb/nSitesHyb), validSites = validSitesHyb)
  outDF <- rbind(outDF, temp)
}
write.table(outDF, "Scer_piValues.txt", quote=F, sep="\t")

wineHyb <- diversity.stats(wineHyb, pi=TRUE)
beer2hyb <- diversity.stats(beer2hyb, pi=TRUE)
nSitesWine <- wineHyb@n.sites
validSitesWine <- wineHyb@n.valid.sites
piValuesWine <- wineHyb@Pi
nSitesBeer2 <- beer2hyb@n.sites
validSitesBeer2 <- beer2hyb@n.valid.sites
piValuesBeer2 <- beer2hyb@Pi
temp <- data.frame(pop = "wineHybs", pi = 100*(piValuesWine/validSitesWine), validGenomePer=100*(validSitesWine/nSitesWine), validSites = validSitesWine)
temp <- rbind(temp, data.frame(pop = "Beer2ybs", pi = 100*(piValuesBeer2/validSitesBeer2), validGenomePer=100*(validSitesBeer2/nSitesBeer2), validSites = validSitesBeer2))


chrNoHyb <- diversity.stats(chrNoHyb, pi=TRUE)
chrWhyb <- diversity.stats(chrWhyb, pi=TRUE)
WGnoHyb <- diversity.stats(WGnoHyb, pi=TRUE)
WGwHyb <- diversity.stats(WGwHyb, pi=TRUE)
nSitesNHc <- chrNoHyb@n.sites
validSitesNHc <- chrNoHyb@n.valid.sites
piValuesNHc<- chrNoHyb@Pi
nSitesHybc <- chrWhyb@n.sites
validSitesHybc <- chrWhyb@n.valid.sites
piValuesHybc <- chrWhyb@Pi
nSitesNHwg <- WGnoHyb@n.sites
validSitesNHwg <- WGnoHyb@n.valid.sites
piValuesNHwg <- WGnoHyb@Pi
nSitesHybwg <- WGwHyb@n.sites
validSitesHybwg <- WGwHyb@n.valid.sites
piValuesHybwg <- WGwHyb@Pi
outDF <- data.frame()
for (i in 1:length(noHybPops)) {
  temp <- data.frame(pop = noHybPops[i], pi = 100*(piValuesNHwg[i]/validSitesNHwg), validGenomePer=100*(validSitesNHwg/nSitesNHwg), validSites = validSitesNHwg)
  temp <- rbind(temp, data.frame(pop = noHybPops[i], pi = 100*(piValuesNHc[i]/validSitesNHc), validGenomePer=100*(validSitesNHc/nSitesNHc), validSites = validSitesNHc))
  outDF <- rbind(outDF, temp)
}
for (j in 1:length(hybPops)){
  temp <- data.frame(pop = hybPops[j], pi = 100*(piValuesHybwg[j]/validSitesHybwg), validGenomePer=100*(validSitesHybwg/nSitesHybwg), validSites = validSitesHybwg)
  temp <- rbind(temp, data.frame(pop = hybPops[j], pi = 100*(piValuesHybc[j]/validSitesHybc), validGenomePer=100*(validSitesHybc/nSitesHybc), validSites = validSitesHybc))
  outDF <- rbind(outDF, temp)
}
write.table(outDF, "Scer_chrWG_piValues.txt", quote=F, sep="\t")


####Seub###
popWG <- readData("popWG", include.unknown=FALSE)
popWG <- diversity.stats(popWG, pi=TRUE)
nSitesPops <- popWG@n.sites
validSitesPops <- popWG@n.valid.sites
piValuesPops <- popWG@Pi
genomCov <- popWG@n.valid.sites/popWG@n.sites
adjustPi <- 100*(popWG@Pi/popWG@n.valid.sites)
popSummary <- data.frame(rawPi = piValuesPops, adjPi = adjustPi, validGenomePer=genomCov, validSites=validSitesPops, nSites=nSitesPops)

popWGunknown <- readData("popWG", include.unknown=T)
popWGunknown <- diversity.stats(popWGunknown, pi=TRUE)
nSitesPopsUN <- popWGunknown@n.sites
validSitesPopsUN <- popWGunknown@n.valid.sites
piValuesPopsUN <- popWGunknown@Pi
genomCovUN <- popWGunknown@n.valid.sites/popWGunknown@n.sites
adjustPiUN <- 100*(popWGunknown@Pi/popWGunknown@n.valid.sites)
popSummaryUN <- data.frame(adjPi = adjustPiUN, rawPi = piValuesPopsUN, validGenomePer=genomCovUN, validSites=validSitesPopsUN, nSites=nSitesPopsUN)



WGnoHyb <- readData("WGnoHybs", include.unknown=FALSE)
WGwHyb <- readData("WGwHybs", include.unknown=FALSE)
WGlessFro <- readData("WGlessFro", include.unknown=FALSE)
chrNoHyb <- readData("chrNoHybs", include.unknown=FALSE)
chrWhyb <- readData("chrWhybs", include.unknown=FALSE)
SeubNoHyb <- readData("noHybs", include.unknown=FALSE)
SeubWhyb <- readData("wHybs", include.unknown=FALSE)
save.image(file="Seub_popGenome.RData")
noHybPops <- c("all", "PopA", "PopB", "Hol")
hybPops <- c("PopBhyb", "HolHyb", "Hybs", "complexHybs", "SuvaXSeub", "lager", "lagerSaaz", "lagerFrohberg")
all <- list(scan("Seub_all.txt", what="", sep="\n"))
allHyb <- list(scan("Seub_allHyb.txt", what="", sep="\n"))
popA <- list(scan("Seub_popA.txt", what="", sep="\n"))
popB <- list(scan("Seub_popB.txt", what="", sep="\n"))
popBhyb <- list(scan("Seub_popBhyb.txt", what="", sep="\n"))
hol <- list(scan("Seub_Hol.txt", what="", sep="\n"))
holHyb <- list(scan("Seub_HolHyb.txt", what="", sep="\n"))
comHyb <- list(scan("Seub_complexHybs.txt", what="", sep="\n"))
SuvaXSeub <- list(scan("Seub_SuvaXSeub.txt", what="", sep="\n"))
lager <- list(scan("Seub_lager.txt", what="", sep="\n"))
lagerS <- list(scan("Seub_lagerSaaz.txt", what="", sep="\n"))
lagerF <- list(scan("Seub_lagerFrohberg.txt", what="", sep="\n"))
WGnoHyb <- set.populations(WGnoHyb, c(all, popA, popB, hol))
WGnoA <- set.populations(WGnoA, c(popBhyb, holHyb, comHyb, SuvaXSeub, lager, lagerS, lagerF))
WGnoHyb <- diversity.stats(WGnoHyb, pi=TRUE)
WGnoA <- diversity.stats(WGnoA, pi=TRUE)
SeubNoHyb <- set.populations(SeubNoHyb, c(all, popA, popB, hol))
SeubNoHyb <- diversity.stats(SeubNoHyb, pi=TRUE)
nSitesNH <- SeubNoHyb@n.sites
validSitesNH <- SeubNoHyb@n.valid.sites
piValuesNH <- SeubNoHyb@Pi
nSitesHyb <- WGnoA@n.sites
validSitesHyb <- WGnoA@n.valid.sites
piValuesHyb <- WGnoA@Pi
outDF <- data.frame()
for (i in 1:length(noHybPops)) {
  temp <- data.frame(pop = noHybPops[i], pi = 100*(piValuesNH[i]/validSitesNH), validGenomePer=100*(validSitesNH/nSitesNH), validSites = validSitesNH)
  outDF <- rbind(outDF, temp)
}
for (j in 1:length(hybPops)){
  temp <- data.frame(pop = hybPops[j], pi = 100*(piValuesHyb[j]/validSitesHyb), validGenomePer=100*(validSitesHyb/nSitesHyb), validSites = validSitesHyb)
  outDF <- rbind(outDF, temp)
}
write.table(outDF, "Seub_piValues.txt", quote=F, sep="\t")


###Skud###
chrNoHyb <- readData("chrNoHybs", include.unknown=FALSE)
chrWhyb <- readData("chrWhybs", include.unknown=FALSE)
WGnoHyb <- readData("WGnoHybs", include.unknown=FALSE)
WGwHyb <- readData("WGwHybs", include.unknown=FALSE)
SkudNoHyb <- readData("noHybs", include.unknown=FALSE)
SkudWhyb <- readData("wHybs", include.unknown=FALSE)
save.image(file="Skud_popGenome.RData")
noHybPops <- c("all", "EU")
hybPops <- c("allHyb", "EUhyb", "hyb")
all <- list(scan("Skud_noHybAll.txt", what="", sep="\n"))
EU <- list(scan("Skud_noHybEU.txt", what="", sep="\n"))
allHyb <- list(scan("Skud_hybAll.txt", what="", sep="\n"))
EUhyb <- list(scan("Skud_hybEU.txt", what="", sep="\n"))
hyb <- list(scan("Skud_hybOnly.txt", what="", sep="\n"))
SkudNoHyb <- set.populations(SkudNoHyb, c(all, EU))
SkudWhyb <- set.populations(SkudWhyb, c(allHyb, EUhyb, hyb))
chrNoHyb <- set.populations(chrNoHyb, c(all, EU))
chrWhyb <- set.populations(chrWhyb, c(allHyb, EUhyb, hyb))
WGnoHyb <- set.populations(WGnoHyb, c(all, EU))
WGwHyb <- set.populations(WGwHyb, c(allHyb, EUhyb, hyb))
SkudNoHyb <- diversity.stats(SkudNoHyb, pi=TRUE)
SkudWhyb <- diversity.stats(SkudWhyb, pi=TRUE)
nSitesNH <- SkudNoHyb@n.sites
validSitesNH <- SkudNoHyb@n.valid.sites
piValuesNH <- SkudNoHyb@Pi
nSitesHyb <- SkudWhyb@n.sites
validSitesHyb <- SkudWhyb@n.valid.sites
piValuesHyb <- SkudWhyb@Pi
outDFcom <- data.frame()
for (i in 1:length(noHybPops)) {
  temp <- data.frame(pop = noHybPops[i], pi = 100*(piValuesNH[i]/validSitesNH), rawPi=piValuesNH[i], validGenomePer=100*(validSitesNH/nSitesNH), validSites = validSitesNH, nSites=nSitesNH)
  outDFcom <- rbind(outDF, temp)
}
for (j in 1:length(hybPops)){
  temp <- data.frame(pop = hybPops[j], pi = 100*(piValuesHyb[j]/validSitesHyb), rawPi=piValuesHyb[j], validGenomePer=100*(validSitesHyb/nSitesHyb), validSites = validSitesHyb, nSites=nSitesHyb)
  outDFcom <- rbind(outDF, temp)
}
#write.table(outDFcom, "Skud_piValues.txt", quote=F, sep="\t")

chrNoHyb <- diversity.stats(chrNoHyb, pi=TRUE)
chrWhyb <- diversity.stats(chrWhyb, pi=TRUE)
WGnoHyb <- diversity.stats(WGnoHyb, pi=TRUE)
WGwHyb <- diversity.stats(WGwHyb, pi=TRUE)
nSitesNHc <- chrNoHyb@n.sites
validSitesNHc <- chrNoHyb@n.valid.sites
piValuesNHc<- chrNoHyb@Pi
nSitesHybc <- chrWhyb@n.sites
validSitesHybc <- chrWhyb@n.valid.sites
piValuesHybc <- chrWhyb@Pi
nSitesNHwg <- WGnoHyb@n.sites
validSitesNHwg <- WGnoHyb@n.valid.sites
piValuesNHwg <- WGnoHyb@Pi
nSitesHybwg <- WGwHyb@n.sites
validSitesHybwg <- WGwHyb@n.valid.sites
piValuesHybwg <- WGwHyb@Pi
outDF <- data.frame()
for (i in 1:length(noHybPops)) {
  temp <- data.frame(pop = noHybPops[i], pi = 100*(piValuesNHwg[i]/validSitesNHwg), validGenomePer=100*(validSitesNHwg/nSitesNHwg), validSites = validSitesNHwg)
  temp <- rbind(temp, data.frame(pop = noHybPops[i], pi = 100*(piValuesNHc[i]/validSitesNHc), validGenomePer=100*(validSitesNHc/nSitesNHc), validSites = validSitesNHc))
  outDF <- rbind(outDF, temp)
}
for (j in 1:length(hybPops)){
  temp <- data.frame(pop = hybPops[j], pi = 100*(piValuesHybwg[j]/validSitesHybwg), validGenomePer=100*(validSitesHybwg/nSitesHybwg), validSites = validSitesHybwg)
  temp <- rbind(temp, data.frame(pop = hybPops[j], pi = 100*(piValuesHybc[j]/validSitesHybc), validGenomePer=100*(validSitesHybc/nSitesHybc), validSites = validSitesHybc))
  outDF <- rbind(outDF, temp)
}
write.table(outDF, "Skud_chrWG_piValues.txt", quote=F, sep="\t")


###Suva###
chrNoHyb <- readData("chrNoHybs", include.unknown=FALSE)
chrWhyb <- readData("chrWhybs", include.unknown=FALSE)
WGnoHyb <- readData("WGnoHybs", include.unknown=FALSE)
WGwHyb <- readData("WGwHybs", include.unknown=FALSE)
SuvaNoHyb <- readData("noHybs", include.unknown=FALSE)
SuvaWhyb <- readData("wHybs", include.unknown=FALSE)
save.image(file="Suva_popGenome.RData")
noHybPops <- c("all", "Aust", "noAust", "PopB", "PopA", "PopAhol", "hol")
hybPops <- c("allHybs", "PopAholHybs", "holHybs", "complexHybs")
all <- list(scan("Suva_all.txt", what="", sep="\n"))
allHyb <- list(scan("Suva_allHyb.txt", what="", sep="\n"))
aust <- list(scan("Suva_Aust.txt", what="", sep="\n"))
noAust <- list(scan("Suva_noAust.txt", what="", sep="\n"))
popB <- list(scan("Suva_popB.txt", what="", sep="\n"))
popA <- list(scan("Suva_popA.txt", what="", sep="\n"))
popAhol <- list(scan("Suva_popAhol.txt", what="", sep="\n"))
hol <- list(scan("Suva_hol.txt", what="", sep="\n"))
popAholHybs <- list(scan("Suva_popAholHybs.txt", what="", sep="\n"))
holHybs <- list(scan("Suva_holHybs.txt", what="", sep="\n"))
hybs <- list(scan("Suva_hybs.txt", what="", sep="\n"))
WGnoHyb <- set.populations(WGnoHyb, c(all, aust, noAust, popB, popA, popAhol, hol))
WGwHyb <- set.populations(WGwHyb, c(allHyb, popAholHybs, holHybs, hybs))
WGnoHyb <- diversity.stats(WGnoHyb, pi=TRUE)
WGwHyb <- diversity.stats(WGwHyb, pi=TRUE)
nSitesNH <- WGnoHyb@n.sites
validSitesNH <- WGnoHyb@n.valid.sites
piValuesNH <- WGnoHyb@Pi
nSitesHyb <- WGwHyb@n.sites
validSitesHyb <- WGwHyb@n.valid.sites
piValuesHyb <- WGwHyb@Pi
outDF <- data.frame()
for (i in 1:length(noHybPops)) {
  temp <- data.frame(pop = noHybPops[i], pi = 100*(piValuesNH[i]/validSitesNH), validGenomePer=100*(validSitesNH/nSitesNH), validSites = validSitesNH)
  outDF <- rbind(outDF, temp)
}
for (j in 1:length(hybPops)){
  temp <- data.frame(pop = hybPops[j], pi = 100*(piValuesHyb[j]/validSitesHyb), validGenomePer=100*(validSitesHyb/nSitesHyb), validSites = validSitesHyb)
  outDF <- rbind(outDF, temp)
}
write.table(outDF, "Suva_piValues.txt", quote=F, sep="\t")


SuvaNoHyb <- set.populations(SuvaNoHyb, c(all, aust, noAust, popB, popA, popAhol, hol))
SuvaWhyb <- set.populations(SuvaWhyb, c(allHyb, popAholHybs, holHybs, hybs))
SuvaNoHyb <- diversity.stats(SuvaNoHyb, pi=TRUE)
SuvaWhyb <- diversity.stats(SuvaWhyb, pi=TRUE)
nSitesNH <- SuvaNoHyb@n.sites
validSitesNH <- SuvaNoHyb@n.valid.sites
piValuesNH <- SuvaNoHyb@Pi
nSitesHyb <- SuvaWhyb@n.sites
validSitesHyb <- SuvaWhyb@n.valid.sites
piValuesHyb <- SuvaWhyb@Pi
outDF <- data.frame()
for (i in 1:length(noHybPops)) {
  temp <- data.frame(pop = noHybPops[i], pi = 100*(piValuesNH[i]/validSitesNH), validGenomePer=100*(validSitesNH/nSitesNH), validSites = validSitesNH)
  outDF <- rbind(outDF, temp)
}
for (j in 1:length(hybPops)){
  temp <- data.frame(pop = hybPops[j], pi = 100*(piValuesHyb[j]/validSitesHyb), validGenomePer=100*(validSitesHyb/nSitesHyb), validSites = validSitesHyb)
  outDF <- rbind(outDF, temp)
}
write.table(outDF, "Suva_chrWG_piValues.txt", quote=F, sep="\t")
