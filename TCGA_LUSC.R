library(ssGSEA4TCGA)
library(ConsensusClusterPlus)
library(NMF)
library(Cairo)

###### TCGA_LUSC ######
gs <- gs_gmt ('bindea_gmt')

rdate <- getFirehoseRunningDates(last = NULL)

dset <- getFirehoseDatasets()
dset <- 'LUSC'
ES_TCGA_LUSC <- ssGSEA(dset, rdate[1], gs)
ES_TCGA_LUSC <- t(ES_TCGA_LUSC)


###### TCGA_LUAD ######


dset <- 'LUAD'
ES_TCGA_LUAD <- ssGSEA(dset, rdate[1], gs)
colnames(ES_TCGA_LUAD)
ES_TCGA_LUAD <- t(ES_TCGA_LUAD)

###### Combind ad and sq ###############
dataset <- c('LUAD', 'LUSC')
ES_TCGA_combind <- pancancer_ssGSEA(dataset, rdate[1], gs)
ES_TCGA_combind <- t(ES_TCGA_combind)
###### Clustering ########

###### LUAD ##############
ES_TCGA_LUAD <- scale(ES_TCGA_LUAD, center = TRUE, scale = TRUE)

results_LUAD_col = ConsensusClusterPlus(ES_TCGA_LUAD,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                        title='consensus_TCGA_luad_col',
                                        clusterAlg="hc",
                                        innerLinkage = "ward.D2",
                                        finalLinkage = "ward.D2",
                                        distance="euclidean",
                                        plot="pdf")

results_LUAD_row = ConsensusClusterPlus(t(ES_TCGA_LUAD),maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                        title='consensus_TCGA_luad_row',
                                        clusterAlg="hc",
                                        innerLinkage = "ward.D2",
                                        finalLinkage = "ward.D2",
                                        distance="euclidean",
                                        plot="pdf")


rownames(ES_TCGA_LUAD)[1:10]
ann <- data.frame(case_class = factor(results_LUAD_row[[4]]$consensusClass))
#ann$TN[grep('TCGA.........11*', rownames(ES))] <- 'Normal'
#ann$TN[is.na(ann$TN)] <- 'Tumor'

ann_col <- data.frame(immune_class = factor(results_LUAD_col[[3]]$consensusClass))
ann_colors <- list(case_class = c('yellow', 'blue', 'orange', 'purple')) #, 'green', 'darkslategray1'),
#TN = c('darkolivegreen','deepskyblue'))
CairoPDF(file = 'immune_TCGA_LUAD.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES_TCGA_LUAD, 
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

write.table(ES_TCGA_LUAD, file = "immune_TCGA_LUAD.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE)

###### LUSC ##############
ES_TCGA_LUSC <- scale(ES_TCGA_LUSC, center = TRUE, scale = TRUE)

results_LUSC_col = ConsensusClusterPlus(ES_TCGA_LUSC,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                        title='consensus_TCGA_lusq_col',
                                        clusterAlg="hc",
                                        innerLinkage = "ward.D2",
                                        finalLinkage = "ward.D2",
                                        distance="euclidean",
                                        plot="pdf")

results_LUSC_row = ConsensusClusterPlus(t(ES_TCGA_LUSC),maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                        title='consensus_TCGA_lusq_row',
                                        clusterAlg="hc",
                                        innerLinkage = "ward.D2",
                                        finalLinkage = "ward.D2",
                                        distance="euclidean",
                                        plot="pdf")


rownames(ES_TCGA_LUSC)[1:10]
ann <- data.frame(case_class = factor(results_LUSC_row[[4]]$consensusClass))
#ann$TN[grep('TCGA.........11*', rownames(ES))] <- 'Normal'
#ann$TN[is.na(ann$TN)] <- 'Tumor'

ann_col <- data.frame(immune_class = factor(results_LUSC_col[[3]]$consensusClass))
ann_colors <- list(case_class = c('yellow', 'blue', 'orange', 'purple')) #, 'green', 'darkslategray1'),
#TN = c('darkolivegreen','deepskyblue'))
CairoPDF(file = 'immune_TCGA_LUSC.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES_TCGA_LUSC, 
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
              Colv = results_LUSC_col[[3]]$consensusTree,
              Rowv = results_LUSC_row[[4]]$consensusTree
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()

write.table(results_LUSC_row[[4]]$consensusClass, 'class_TCGA_LUSC.txt',
            row.names = TRUE,
            col.names = FALSE)

write.table(ES_TCGA_LUSC, file = "immune_TCGA_LUSC.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE)

########################################
############# Combined #################
########################################

ES_TCGA_combind <- scale(ES_TCGA_combind, center = TRUE, scale = TRUE)

results_combind_col = ConsensusClusterPlus(ES_TCGA_combind,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                        title='consensus_TCGA_combind_col',
                                        clusterAlg="hc",
                                        innerLinkage = "ward.D2",
                                        finalLinkage = "ward.D2",
                                        distance="euclidean",
                                        plot="pdf")

results_combind_row = ConsensusClusterPlus(t(ES_TCGA_combind),maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                        title='consensus_TCGA_combind_row',
                                        clusterAlg="hc",
                                        innerLinkage = "ward.D2",
                                        finalLinkage = "ward.D2",
                                        distance="euclidean",
                                        plot="pdf")


rownames(ES_TCGA_combind)[1:10]
ann <- data.frame(case_class = factor(results_combind_row[[4]]$consensusClass))
#ann$TN[grep('TCGA.........11*', rownames(ES))] <- 'Normal'
#ann$TN[is.na(ann$TN)] <- 'Tumor'

ann_col <- data.frame(immune_class = factor(results_combind_col[[3]]$consensusClass))
ann_colors <- list(case_class = c('yellow', 'blue', 'orange', 'purple')) #, 'green', 'darkslategray1'),
#TN = c('darkolivegreen','deepskyblue'))
CairoPDF(file = 'immune_TCGA_combind.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES_TCGA_combind, 
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
              Colv = results_combind_col[[3]]$consensusTree,
              Rowv = results_combind_row[[4]]$consensusTree
              #reorderfun = function(d, w) reorder(d, 10)
)

dev.off()



write.table(ES_TCGA_combind, file = "immune_TCGA_combine.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE)


###### Clnical data #####

getFirehoseData('LUAD', runDate = rdate[1], gistic2_Date = NULL,
                                 RNAseq_Gene = FALSE, Clinic = TRUE, miRNASeq_Gene = FALSE,
                                 RNAseq2_Gene_Norm = FALSE, CNA_SNP = FALSE, CNV_SNP = FALSE,
                                 CNA_Seq = FALSE, CNA_CGH = FALSE, Methylation = FALSE,
                                 Mutation = FALSE, mRNA_Array = FALSE, miRNA_Array = FALSE,
                                 RPPA = FALSE, RNAseqNorm = "raw_counts",
                                 RNAseq2Norm = "normalized_count", forceDownload = FALSE, destdir = ".",
                                 fileSizeLimit = 500, getUUIDs = FALSE)

Clinical_LUAD_df <- getData(clinical_LUAD,"Clinical")
Clinical_LUAD_df <- survivalTCGA(Clinical_LUAD_df)
Clinical_LUAD_df$ID <- gsub('\\.', '-', Clinical_LUAD_df$ID)
Clinical_LUAD_df$ID %in% rownames(ES_TCGA_LUAD)
data_LUAD <- data.frame(ID = rownames(ES_TCGA_LUAD),
           class = results_LUAD_row[[4]]$consensusClass,
           ES_TCGA_LUAD)

data_LUAD <- merge(Clinical_LUAD_df, data_LUAD, by = 'ID', all.x = FALSE)

write.table(data_LUAD, file = "data_TCGA_LUAD.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

############ LUSC #######################
clinical_LUSC <- getFirehoseData('LUSC', runDate = rdate[1], gistic2_Date = NULL,
                RNAseq_Gene = FALSE, Clinic = TRUE, miRNASeq_Gene = FALSE,
                RNAseq2_Gene_Norm = FALSE, CNA_SNP = FALSE, CNV_SNP = FALSE,
                CNA_Seq = FALSE, CNA_CGH = FALSE, Methylation = FALSE,
                Mutation = FALSE, mRNA_Array = FALSE, miRNA_Array = FALSE,
                RPPA = FALSE, RNAseqNorm = "raw_counts",
                RNAseq2Norm = "normalized_count", forceDownload = FALSE, destdir = ".",
                fileSizeLimit = 500, getUUIDs = FALSE)

Clinical_LUSC_df <- getData(clinical_LUSC,"Clinical")
Clinical_LUSC_df <- survivalTCGA(Clinical_LUSC_df)
Clinical_LUSC_df$ID <- gsub('\\.', '-', Clinical_LUSC_df$ID)
Clinical_LUSC_df$ID %in% rownames(ES_TCGA_LUSC)
data_LUSC <- data.frame(ID = rownames(ES_TCGA_LUSC),
                        class = results_LUSC_row[[4]]$consensusClass,
                        ES_TCGA_LUSC)

data_LUSC <- merge(Clinical_LUSC_df, data_LUSC, by = 'ID', all.x = FALSE)

write.table(data_LUSC, file = "data_TCGA_LUSC.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)


image_id <- read.delim('TCIA case no_TCGA LUSC.txt', header = FALSE)

write.table(data_LUSC [data_LUSC$ID %in% image_id$V1,], file = "data_TCGA_LUSC_with_image.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)


#################################################
###################### CLASS transform ##########
library(reshape2)
data_LUAD$class == 3
combine_factor(data_LUAD$class)
class_seperate <- data.frame(ID = append(data_LUAD$ID, data_LUSC$ID), class = append(data_LUAD$class == 3, data_LUSC$class == 4))
class_seperate$class <- factor(class_seperate$class, levels = c(TRUE, FALSE), labels = c('Immune suppressive', 'immune active'))

class_combine <- data.frame(ID = rownames(ES_TCGA_combind), class = results_combind_row[[4]]$consensusClass)
class_seperate$ID == class_combine$ID 
class_combine$class <- factor(class_combine$class == 3, levels = c(TRUE, FALSE), labels = c('Immune suppressive', 'immune active'))

table(class_seperate$class[1:515], class_combine$class[1:515])
table(class_seperate$class[516:1016], class_combine$class[516:1016])
chisq.test(class_seperate$class, class_combine$class)

write.table(class_combine, 'class_TCGA_combined.txt',
            row.names = FALSE,
            col.names = TRUE)
write.table(class_seperate, 'class_TCGA_seperate.txt',
            row.names = FALSE,
            col.names = TRUE)
