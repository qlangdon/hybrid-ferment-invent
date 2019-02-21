args <- commandArgs(TRUE)
strains <- read.table(args[1])
output <- data.frame()
for (i in 1:length(strains[,1])) {
  strainName = toString(strains[i,1])
  strainFile <- read.table(paste(strainName, ".bedgraph", sep=""), header=FALSE, col.names = c("chrom", "chromstart", "value"))
  strainOut <- data.frame(strain=NA, q10=NA, q99=NA)
  strainOut$strain <- strainName
  strainOut$q10<- quantile(strainFile$value, 0.1)
  strainOut$q99<- quantile(strainFile$value, 0.99)
  strainOut$mean<- mean(strainFile$value)
  strainOut$median <- median(strainFile$value)
  output <- rbind(output, strainOut)
}
write.table(output, file="coverageQuant_forMask", row.names=F, sep = "\t")