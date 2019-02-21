options(stringsAsFactors = FALSE)
library("seqinr")
library("igraph")
library("reshape")

PAD1 = read.alignment("PAD1_moreSuvaXSeub_AA_psiTrimmed_noSingles.fas", format = "FASTA", forceToLower = F)
HaploKey = read.delim("PAD1_haplotypeKey (1).txt", sep = "\t", header = F)
colnames(HaploKey) = c("Strain", "Origin", "Species", "Population", "Pop2")

PAD1_m = as.matrix.alignment(PAD1)
colnames(PAD1_m) = c(1:ncol(PAD1_m))

SameCount = matrix(0, nrow = 1, ncol = ncol(PAD1_m))
colnames(SameCount) = c(1:ncol(SameCount))
for(i in 1:ncol(PAD1_m)){
  SameCount[1,i] = length(which(unlist(PAD1_m[,i] == PAD1_m[1,i]) == TRUE))
}

VARIABLESITES = as.numeric(colnames(SameCount)[which(SameCount[1,] != max(SameCount[1,]))])

PAD1_var = PAD1_m[,VARIABLESITES]

#write.csv(PAD1_var, file = "PAD1_var.csv")

################## Find Haplotypes ######################
STNAMES = rownames(PAD1_var)
SITES = colnames(PAD1_var)

Haplotype_INFO = matrix(0, nrow = length(STNAMES), ncol = length(STNAMES))
rownames(Haplotype_INFO) = STNAMES
colnames(Haplotype_INFO) = STNAMES
for(i in 1:length(STNAMES)){
  for(j in 1:length(STNAMES)){
    TEMP = PAD1_var[c(i,j),]
    COMPARISON = TEMP[2,] == TEMP[1,]
    Haplotype_INFO[i,j] = length(which(COMPARISON == FALSE))
  }
}

HaploTEMP = melt(Haplotype_INFO)
colnames(HaploTEMP) = c("Haplotype_A", "Haplotype_B", "Differences")

Haplotype_INFO_lim = HaploTEMP[which(HaploTEMP$Differences == 0),]

Haplo_list =list()
for(i in 1:length(STNAMES)){
  TEMP1 = Haplotype_INFO_lim[which(Haplotype_INFO_lim$Haplotype_A == STNAMES[[i]]),2]
  TEMP2 = Haplotype_INFO_lim[which(Haplotype_INFO_lim$Haplotype_A == STNAMES[[i]]),1]
  TEMP3 = Haplotype_INFO_lim[which(Haplotype_INFO_lim$Haplotype_B == STNAMES[[i]]),2]
  TEMP4 = Haplotype_INFO_lim[which(Haplotype_INFO_lim$Haplotype_B == STNAMES[[i]]),1]
  Haplo_list[[i]] = list()
  Haplo_list[[i]] = unique(c(TEMP1, TEMP2, TEMP3, TEMP4))
}
names(Haplo_list) = STNAMES

SetDiff_m = matrix(0, nrow = length(names(Haplo_list)), ncol = length(names(Haplo_list)))
rownames(SetDiff_m) = names(Haplo_list)
colnames(SetDiff_m) = names(Haplo_list)
for(i in 1:length(Haplo_list)){
  for(j in 1:length(Haplo_list)){
    if(length(Haplo_list[[i]]) < (ncol(SetDiff_m) - 1) | length(Haplo_list[[j]]) < (ncol(SetDiff_m)-1)){
      TEMP1 =setdiff(Haplo_list[[i]], Haplo_list[[j]])  
      TEMP2 =setdiff(Haplo_list[[j]], Haplo_list[[i]]) 
      TEMP =unique(c(TEMP1, TEMP2))
      SetDiff_m[i,j] = length(TEMP)
    }else{
      SetDiff_m[i,j] = ncol(SetDiff_m)-1
    }
  }
}

Haplotypes_df = data.frame(HaploGroup = character(), Count = numeric())
HaploGroup_list = list()
TEMPNAMES = STNAMES

k = 1
while(length(TEMPNAMES) > 0){
  HaploGroup_list[[k]] = list()
  TEMP = SetDiff_m[which(rownames(SetDiff_m) == TEMPNAMES[1]),]
  DROPNAMES = names(which(TEMP == 0))
  TEMPNAMES = TEMPNAMES[-which(TEMPNAMES %in% DROPNAMES)]
  HaploGroup_list[[k]] = DROPNAMES
  Haplotypes_df[k,1] = k
  Haplotypes_df[k,2] = length(DROPNAMES)
  k = k+1
}

Representatives = data.frame(Group = numeric(), Strain = character())
for(i in 1:nrow(Haplotypes_df)){
  Representatives[i,1] = i
  Representatives[i,2] = HaploGroup_list[[i]][1]
}

Rep_m = matrix(0, ncol = ncol(PAD1_var), nrow = 1)
#rownames(Rep_m) = Representatives$Group
colnames(Rep_m) = colnames(PAD1_var)
Rep_m = as.data.frame(Rep_m)

for(i in 1:nrow(Representatives)){
  TEMP = t(as.matrix(PAD1_var[which(rownames(PAD1_var) == Representatives[i,2]),]))
  Rep_m = rbind(Rep_m, TEMP)
}
Rep_m = Rep_m[-1,]
rownames(Rep_m) = Representatives$Group

Haplotype_df = data.frame(HaploGroup_A = character(), HaploGroup_B = character(), Differences = numeric())
k = 1
for(i in 1:(nrow(Rep_m)-1)){
  for(j in (i+1):nrow(Rep_m)){
    TEMP = Rep_m[c(i,j),]
    COMPARISON = TEMP[2,] == TEMP[1,]
    Haplotype_df[k,1] = rownames(Rep_m)[i]
    Haplotype_df[k,2] = rownames(Rep_m)[j]
    Haplotype_df[k,3] = length(which(COMPARISON == FALSE))
    k = k+1
  }
}
GROUPS = Representatives$Group
Haplotype_df_lim = data.frame(HaploGroup_A = character(), HaploGroup_B = character(), Differences = numeric())
for(i in 1:(length(GROUPS)-1)){
  TEMP = Haplotype_df[which(Haplotype_df$HaploGroup_A == GROUPS[[i]]),]
  MINTEMP = TEMP[which(TEMP$Differences == min(TEMP$Differences, na.rm = FALSE)),]
  Haplotype_df_lim = rbind(Haplotype_df_lim, MINTEMP)
}

Haplotype_df_lim$lty = rep(1, nrow(Haplotype_df_lim))
Haplotype_df_lim[which(Haplotype_df_lim$Differences >= 100 ), "lty"] = 2

Haplogroup_df = data.frame(Group = character(), Count = numeric())
for(i in 1:length(HaploGroup_list)){
  Haplogroup_df[i,1] = i
  Haplogroup_df[i,2] = length(HaploGroup_list[[i]])
}

Haplogroup_df$Ratio = Haplogroup_df$Count/sum(Haplogroup_df$Count)

Haplogroup_df$Scaled = Haplogroup_df$Ratio
for(i in 1:length(unique(Haplogroup_df$Ratio))){
  Haplogroup_df[which(Haplogroup_df$Ratio == min(Haplogroup_df$Scaled)),"Scaled"] = i 
}

HaploStrain = data.frame(Strain = character(), Haplogroup = character())
k = 1
for(i in 1:length(HaploGroup_list)){
  for(j in 1:length(HaploGroup_list[[i]])){
    HaploStrain[k,1] = HaploGroup_list[[i]][j]
    HaploStrain[k,2] = i
    k = k+1
    
  }
}
write.csv(HaploStrain, file = "HaploStrain_PAD1_AA_condensed_more.csv")

#ORIGIN = unique(HaploKey$Origin)
#ORIGIN = ORIGIN[-which(ORIGIN == 0)]
#pieValues_Origin = list()
# for(i in 1:length(HaploGroup_list)){
#   pieValues_Origin[[i]] = list()
#   TEMP=HaploGroup_list[[i]]
#   TEMPVECTOR = vector()
#   for(j in 1:length(ORIGIN)){
#     TEMPGROUP = HaploKey[which(HaploKey$Strain %in% TEMP),]
#     TEMPVALUE = length(which(TEMPGROUP$Origin == ORIGIN[[j]]))
#     TEMPVECTOR = c(TEMPVECTOR, TEMPVALUE) 
#   }
#   pieValues_Origin[[i]] = TEMPVECTOR
# 
# }


# SPECIES = unique(HaploKey$Species)
# pieValues_Species = list()
# for(i in 1:length(HaploGroup_list)){
#   pieValues_Species[[i]] = list()
#   TEMP=HaploGroup_list[[i]]
#   TEMPVECTOR = vector()
#   for(j in 1:length(SPECIES)){
#     TEMPGROUP = HaploKey[which(HaploKey$Strain %in% TEMP),]
#     TEMPVALUE = length(which(TEMPGROUP$Species == SPECIES[[j]]))
#     TEMPVECTOR = c(TEMPVECTOR, TEMPVALUE) 
#   }
#   pieValues_Species[[i]] = TEMPVECTOR
#   
# }
# 

POPULATION = unique(HaploKey$Pop2)
pieValues_Population = list()
for(i in 1:length(HaploGroup_list)){
  pieValues_Population[[i]] = list()
  TEMP=HaploGroup_list[[i]]
  TEMPVECTOR = vector()
  for(j in 1:length(POPULATION)){
    TEMPGROUP = HaploKey[which(HaploKey$Strain %in% TEMP),]
    TEMPVALUE = length(which(TEMPGROUP$Pop2 == POPULATION[[j]]))
    TEMPVECTOR = c(TEMPVECTOR, TEMPVALUE) 
  }
  pieValues_Population[[i]] = TEMPVECTOR
  
}

Edges_Graph = Haplotype_df_lim[,c(1,2)]
g = graph.data.frame(Edges_Graph, directed = FALSE, vertices = Haplogroup_df)

NodeSize = (V(g)$Scaled*2)
NodeLabel = names(V(g))

minC <- rep(-Inf, vcount(g)*100)
maxC <- rep(Inf, vcount(g))
minC[1] <- maxC[1] <- 0
fc = fastgreedy.community(g, weights = Haplotype_df_lim$Differences)

c1 = cluster_fast_greedy(g, weights = Haplotype_df_lim$Differences)

coords <- layout_with_dh(g)
coords <- norm_coords(coords, ymin=-10, ymax=10, xmin=-10, xmax=10)

colors <- c("#0092FF", "#FF0000", "#00FF92", "#FFDB00", "#707070")
#colors <- rainbow(max(membership(fc)))

V(g)$pie.color = list(c("#CC79A7",
                        "#999999",
                        "#868686",
                        "#707070",
                        "#E8E8E8",
                        "#D2D2D2",
                        "#B8B8B8",
                        "#DEDEDE",
                        "#C6C6C6",
                        "#525252",
                        "#E1B0CB",
                        "#0092FF",
                        "#000000",
                        "#AAAAAA",
                        "#F3F3F3"))
pdf("HaplotypeNetwork_Population_AA_PAD1_community_condensed_more_labeled.pdf", height = 28, width = 25)
plot(c1,g,mark.col= adjustcolor(colors, alpha = 0.5), mark.border = colors,
     vertex.shape = "pie",vertex.pie = pieValues_Population,
#     vertex.label = NA,
    vertex.label = NodeLabel,
     vertex.label.cex = 2,
     vertex.label.color = "cyan",
     vertex.size = NodeSize*10,
     edge.label = Haplotype_df_lim$Differences,
     edge.label.color = "black",
     edge.label.cex = 2,
     edge.width = 2,
     edge.color="black", 
     edge.lty = Haplotype_df_lim$lty,
     layout=coords, xlim = range(coords[,1]), ylim = range(coords[,2]), rescale=F)
dev.off()

save(list = ls(),file = "PAD1_AA_Haplotype.RData")

plot(1,1)

POPCOLS = c("#0092FF",  
            "#000000",
            "#525252",
            "#707070",
            "#868686",
            "#999999",
            "#AAAAAA",
            "#B8B8B8",
            "#C6C6C6",
            "#D2D2D2",
            "#DEDEDE",
            "#E8E8E8",  
            "#F3F3F3", 
            "#E1B0CB",  
            "#CC79A7")
NAMES = c("Non-Scer", 
          "Wild Misc.", 
          "Sake/Asian", 
          "Bread/Mixed", 
          "Beer2", 
          "Wine",
          "MedOak",
          "Mosaic", 
          "Belgian/German/Alt-Kolsch", 
          "British Isles",
          "US", 
          "Ale/Beer1 Mosaic",
          "Wheat", "Lager", 
          "Non-Lager Hybrids")


pdf("Legend.pdf", height = 25, width = 20)
plot(0)
legend("topright", 
       legend=  NAMES , 
      pt.bg = POPCOLS, pch=21 , col = "black", 
       pt.cex = 3.5, cex = 1.5, text.col="black" , horiz = F, 
       inset = c(0.1, 0.1))
dev.off()

# V(g)$pie.color = list(c("#C9C9C9",
#                         "#ff0000",
#                         "#F2F2F2",
#                         "#000000",
#                         "#FFBF00",
#                         "#FF00DB",
#                         "#00FF92",
#                         "#00FF40",
#                         "#B1B1B1",
#                         "#939393",
#                         "#4900FF"))
# 
# 
# pdf("HaplotypeNetwork_Species_AA_PAD1_community.pdf", height = 25, width = 20)
# plot(c1,g,vertex.color=colors[membership(fc)], 
#      vertex.shape = "pie",vertex.pie = pieValues_Species,
#      vertex.label = NodeLabel,
#      vertex.label.cex = 2,
#      vertex.label.color = "black",
#      vertex.size = NodeSize*10,
#      edge.label = Haplotype_df_lim$Differences,
#      edge.label.color = "black",
#      edge.label.cex = 2,
#      edge.width = 2,
#      edge.color="black", 
#      #edge.lty = Haplotype_df_lim$lty,
#      layout=coords, xlim = range(coords[,1]), ylim = range(coords[,2]), rescale=F)
# dev.off()
# 
# 
# V(g)$pie.color = list(c("#cc79a7", "#d52900", "#0072b2"))
# pdf("HaplotypeNetwork_Origin_AA_PAD1_community.pdf", height = 25, width = 20)
# plot(c1,g,vertex.color=colors[membership(fc)], 
#      vertex.shape = "pie",vertex.pie = pieValues_Origin,
#      vertex.label = NodeLabel,
#      vertex.label.cex = 2,
#      vertex.label.color = "black",
#      vertex.size = NodeSize*10,
#      edge.label = Haplotype_df_lim$Differences,
#      edge.label.color = "black",
#      edge.label.cex = 2,
#      edge.width = 2,
#      edge.color="black", 
#      #edge.lty = Haplotype_df_lim$lty,
#      layout=coords, xlim = range(coords[,1]), ylim = range(coords[,2]), rescale=F)
# dev.off()
# 
# 
