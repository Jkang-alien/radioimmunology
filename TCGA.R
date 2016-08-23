library('RTCGAToolbox')


getFirehoseDatasets()
getFirehoseRunningDates(last = NULL)
readData = getFirehoseData (dataset="LUAD", runDate="20151101",forceDownload = TRUE,
                            Clinic=TRUE, Mutation=FALSE, Methylation=FALSE, RNAseq2_Gene_Norm=TRUE)

clin = getData(readData, "Clinical")
clin[1:5,]
mRNA <- getData(readData, 'RNASeq2GeneNorm')
mRNA[1:10, 1:10]
dim(mRNA)
colnames(mRNA)[duplicated(gsub('-...-...-....-..', '', colnames(mRNA)), fromLast = TRUE)]
colnames(mRNA)
sum(duplicated(gsub('-...-...-....-..', '', colnames(mRNA))))

#colnames(mRNA) <- gsub('-...-...-....-..', '', colnames(mRNA))

write.table(mRNA, file = 'TCGA_mixture.txt', append = FALSE, 
            quote = FALSE, sep = '\t',
            row.names = TRUE,
            col.names = TRUE 
)
dim(mRNA)
