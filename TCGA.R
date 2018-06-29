library('RTCGAToolbox')


getFirehoseDatasets()
getFirehoseRunningDates(last = NULL)
readData = getFirehoseData (dataset="LUAD", runDate="20160128",forceDownload = TRUE,
                            Clinic=TRUE, Mutation=FALSE, Methylation=FALSE, RNAseq2_Gene_Norm=TRUE)

clin = getData(readData, "Clinical")
clin[1:5,]
mRNA_adeno <- getData(readData, 'RNASeq2GeneNorm')

mRNA_adeno[1:10, 1:10]
dim(mRNA_adeno)
colnames(mRNA_adeno)[duplicated(gsub('-...-...-....-..', '', colnames(mRNA_adeno)), fromLast = TRUE)]
colnames(mRNA_adeno)
sum(duplicated(gsub('-...-...-....-..', '', colnames(mRNA_adeno))))

#colnames(mRNA_adeno) <- gsub('-...-...-....-..', '', colnames(mRNA_adeno))

write.table(mRNA_adeno, file = 'TCGA_mixture_LUSC.txt', append = FALSE, 
            quote = FALSE, sep = '\t',
            row.names = TRUE,
            col.names = TRUE 
)
dim(mRNA_adeno)

######################################################
############## ssGSEA data analysis ##################
immune <- read.delim('TCGA_LUAD_mixture.PROJ.gct', sep = '\t')


cell_type <- as.character(immune$Description)
ID <- colnames(immune)[3:dim(immune)[2]]

ES <- t(immune[,3:dim(immune)[2]])
ES_TCGA <- scale(ES_TCGA, center = TRUE, scale = TRUE)
#ES <- sweep(ES, 1, apply(ES, 1, median, na.rm = TRUE))
colnames(ES) <- cell_type

rownames(ES)[duplicated(gsub('\\....\\....\\.....\\...', '', rownames(ES)), fromLast = FALSE)]
recur_T <- c("TCGA.50.5946.02A.11R.2090.07", "TCGA.50.5066.02A.11R.2090.07")

library(ConsensusClusterPlus)

results_LUAD_col = ConsensusClusterPlus(ES_TCGA,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_TCGA_luad_col',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

results_LUAD_row = ConsensusClusterPlus(t(ES_TCGA),maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_TCGA_luad_row',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

library(NMF)
library(Cairo)

rownames(ES_TCGA)[1:10]
ann <- data.frame(case_class = factor(results_LUAD_row[[4]]$consensusClass))
#ann$TN[grep('TCGA.........11*', rownames(ES))] <- 'Normal'
#ann$TN[is.na(ann$TN)] <- 'Tumor'

ann_col <- data.frame(immune_class = factor(results_LUAD_col[[3]]$consensusClass))
ann_colors <- list(case_class = c('yellow', 'blue', 'orange', 'purple')) #, 'green', 'darkslategray1'),
                   #TN = c('darkolivegreen','deepskyblue'))
CairoPDF(file = 'immune_TCGA_LUAD.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES_TCGA, 
              hclustfun=function(d) hclust(d, method="ward.D2"),
              #Colv = colv,
              annRow = ann,
              annCol = ann_col,
              annColors = ann_colors,
              #cex = 2,
              #labRow = rep('',dim(data_hc_e)[1]),
              #labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              Colv = results_LUAD_col[[3]]$consensusTree,
              Rowv = results_LUAD_row[[4]]$consensusTree
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()

write.table(results_LUAD_row[[4]]$consensusClass, 'class_TCGA.txt',
            row.names = TRUE,
            col.names = FALSE)

data_clin <- data.frame(ID = gsub('\\....\\....\\.....\\...', '', rownames(ann)),
                        ann)
data_clin <- data_clin[(rownames(data_clin) %in% recur_T == FALSE) & (data_clin$TN == "Tumor"),]

clin$ID <- toupper(rownames(clin))

data <- merge(data_clin, clin, by = 'ID')
class(data$days_to_death)
data[, c(5, 7, 8)] <- sapply(data[, c(5, 7, 8)], as.character)

data[, c(5, 7, 8)] <- sapply(data[, c(5, 7, 8)], as.numeric)
surv_months <- pmax(data$days_to_death, 
                    data$days_to_last_followup,
                    na.rm = TRUE)/30.4

library(survival)
library(rms)


diff = survdiff(Surv(surv_months, vital_status == 1)~ case_class, 
                data = data)
diff

svg(file = "Figure3.svg", pointsize = 10,
    width = 7.5 , height = 4,)
layout(matrix(c(1,2), ncol = 2, byrow = TRUE))
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))

fit = npsurv(Surv(surv_months, vital_status == 1)~ case_class, 
             data = data)

strata = levels(data$case_class)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:6),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data$case_class,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:6), col=data$case_class, cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.052', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')


cox <- coxph(Surv(surv_months, vital_status == 1)~ case_class, 
             data = data)

summary(cox)

######################## LUSC ########################

immune <- read.delim('TCGA_mixture_LUSC.PROJ.gct', sep = '\t', skip = 2)


cell_type <- as.character(immune$Description)
ID <- colnames(immune)[3:dim(immune)[2]]

ES <- t(immune[,3:dim(immune)[2]])
ES_TCGA <- scale(ES_TCGA, center = TRUE, scale = TRUE)
#ES <- sweep(ES, 1, apply(ES, 1, median, na.rm = TRUE))
colnames(ES) <- cell_type


library(ConsensusClusterPlus)

results_LUSC_col = ConsensusClusterPlus(ES,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_TCGA_lusc_col',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

results_LUSC_row = ConsensusClusterPlus(t(ES),maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_TCGA_lusc_row',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")


library(NMF)
library(Cairo)

rownames(ES)[1:10]
ann <- data.frame(case_class = factor(results_LUSC_row[[4]]$consensusClass))
#ann$TN[grep('TCGA.........11*', rownames(ES))] <- 'Normal'
#ann$TN[is.na(ann$TN)] <- 'Tumor'

ann_col <- data.frame(immune_class = factor(results_LUSC_col[[3]]$consensusClass))

CairoPDF(file = 'immune_TCGA_LUSC.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES_TCGA, 
              hclustfun=function(d) hclust(d, method="ward.D2"),
              #Colv = colv,
              annRow = ann,
              annCol = ann_col,
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

