library(pvclust)

hybPloidyDF <- read.csv("allHyb_winPloidy10kb_ScerSkudSuvaSeub.csv", check.names=F, row.names=1)

hybPloidyM <- as.matrix(hybPloidyDF)

fit = pvclust(hybPloidyM, method.hclust = "ward.D", method.dist = "euclidean", nboot = 1000)

pdf("StrainCluster.pdf", width=16)
plot(fit, print.num = FALSE)
pvrect(fit, alpha = 0.95, pv = "au", max.only = FALSE)
dev.off()

pdf("StrainCluster_MaxOnly.pdf", width=16)
plot(fit, print.num = FALSE, print.pv=F)
pvrect(fit, alpha = 0.95, pv = "au", max.only = T)
dev.off()


options(stringAsFactors=FALSE)
ploidyMatrix <- read.csv("allHyb_winPloidy10kb_ScerSkudSuvaSeub.csv", check.name=F, row.names=1)
ploidyT <- t(ploidyMatrix)
ord_All <- hclust( dist(ploidyT, method = "euclidean"), method = "ward.D" )
strainPreSort <- row.names(ploidyT)
ordNums <- ord_All$order
strainOrder <- c()
for (num in ordNums) {
  strainN <- strainPreSort[[num]]
  strainOrder <- c(strainOrder, strainN)
}
revOrder <- rev(strainOrder)

decodingKey <- read.table("ScerSkudSuvaSeub_strainOrder.txt", header=F)
orderDF <- data.frame(genomeName=as.character(), strainKey=as.character())
for (strain in revOrder) {
  strainDF <- data.frame(genomeName=strain, strainKey=decodingKey[which(decodingKey[,1]==strain),2])
  orderDF <- rbind(orderDF, strainDF)
}
write.table(orderDF, "strainOrder_contributionPlot.txt", quote = F, sep="\t")
