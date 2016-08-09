immune <- read.delim('immune_GSE18865.txt', sep = '\t')


cell_type <- as.character(immune$Description)
ID <- colnames(immune)[2:90]

ES <- t(immune[,2:90])
ES <- scale(ES, center = TRUE, scale = TRUE)
#ES <- sweep(ES, 1, apply(ES, 1, median, na.rm = TRUE))
colnames(ES) <- cell_type
## heatmap(ES)

library(ConsensusClusterPlus)

results_col = ConsensusClusterPlus(ES,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                               title='consensus',
                               clusterAlg="hc",
                               innerLinkage = "ward.D2",
                               finalLinkage = "ward.D2",
                               distance="euclidean",
                               plot="pdf")

results_row = ConsensusClusterPlus(t(ES),maxK=10,reps=5000,pItem=0.8,pFeature=1,
                               title='consensus_by_row',
                               clusterAlg="hc",
                               innerLinkage = "ward.D2",
                               finalLinkage = "ward.D2",
                               distance="euclidean",
                               plot="pdf")

library(NMF)
library(Cairo)
CairoPDF(file = 'immune_TCIA_Lung.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES, 
              hclustfun=function(d) hclust(d, method="ward.D2"),
              #Colv = colv,
              #annRow = ann,
              #annCol = ann_col,
              #annColors = ann_colors_e,
              #cex = 2,
              #labRow = rep('',dim(data_hc_e)[1]),
              #labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              Colv = results_col[[3]]$consensusTree,
              Rowv = results_row[[4]]$consensusTree
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()
