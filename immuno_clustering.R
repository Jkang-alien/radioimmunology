immune <- read.delim('ssGSEA_GSE18865.txt')

ssGSEA <- read.delim('TCGA_UM.gct', sep = '\t')
cell_type <- as.character(ssGSEA$Description)
ID <- colnames(ssGSEA)[3:82]

ES <- t(ssGSEA[,3:82])
#ES <- scale(ES, center = TRUE, scale = TRUE)
#ES <- sweep(ES, 1, apply(ES, 1, median, na.rm = TRUE))
colnames(ES) <- cell_type
## heatmap(ES)
Chr3 <- rownames(ES) %in% ID_3qNor
Chr3 <- factor(Chr3, 
               levels = c(TRUE, FALSE),
               labels = c('Disomy', 'Monosomy'))

hc = hclust(dist(ES), 'ward.D2')

group_immune <- factor(cutree(hc, k = 5), levels = 1:5,
                       labels = c('C', 'A', 'B', 'D', 'E'))

ann <- data.frame(Chr3 = Chr3, group = group_immune)

ann$group <- as.character(ann$group)

ann$group[(ann$Chr3 == 'Monosomy') & (ann$group == 'C')] <- 'D'
ann$group <- factor(ann$group)
df_ann <- data.frame(ID = rownames(ann), ann)
df_ES <- data.frame(ID = rownames(ES), ES)

data_ES <- merge(df_ann, df_ES, by = 'ID')

data_ssGEEA_cli <- merge(data, data_ES, by = 'ID')

boxplot(Th2.cells ~ group, data_ES)

par (mar=c(5.1, 5, 1, 1.1))
bp <- boxplot(Th2.cells ~ Chr3, data_ES,
              ylab = 'Enrichment Score',
              lab = '',
              #border = "white",
              frame = FALSE,
              names = rep('',2),
              xaxt='n',
              outpch = NA
              #ylim = c(3,8.5)
              
)

axis(mgp=c(3.5, 2, 0),
     side = 1,
     at = 1:2,
     pos = -8000,
     tck = -0.01,
     labels=c('Disomy',
              "Monosomy"))

stripchart(Th2.cells ~ Chr3, data_ES,
           vertical=T,pch=1, method="jitter",cex=1,add=T)

x <- c(1, 1, 3, 3)
y <- c(7.9, 8.0, 8.0, 7.9)-0.9
lines (x,y)
x <- c(2, 2, 4, 4)
y <- c(7.9, 8.0, 8.0, 7.9)-0.1
lines (x,y)


p_value <- c()
for (i in 4:28){
  p <- summary(aov(data_ES[,i] ~ data_ES$group))[[1]][['Pr(>F)']][1]
  p_value <- append(p_value, p)
}

hist(log(p_value))

colnames(ES)[log(p_value) < -30]

aov <- aov(Th17.cells ~ group, data_ES)
summary(aov)
library(ConsensusClusterPlus)

results = ConsensusClusterPlus(ES,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                               title='consensus',
                               clusterAlg="hc",
                               innerLinkage = "ward.D2",
                               finalLinkage = "ward.D2",
                               distance="euclidean",
                               plot="pdf")

library(NMF)
library(Cairo)
CairoPDF(file = 'ssGSEA_TCGA.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES, 
              hclustfun=function(d) hclust(d, method="ward.D2"),
              #Colv = colv,
              annRow = ann,
              #annCol = ann_col,
              #annColors = ann_colors_e,
              #cex = 2,
              #labRow = rep('',dim(data_hc_e)[1]),
              #labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              Colv = results[[4]]$consensusTree
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()
