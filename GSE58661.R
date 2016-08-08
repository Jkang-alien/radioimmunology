#GSE58661

library(Biobase)
library(GEOquery)

GSEDATA <- getGEO ("GSE58661", GSEMatrix = T, AnnotGPL = FALSE)
print(summary(exprs(GSEDATA[[1]])[, 1:2]))
eset <- GSEDATA[[1]]

###############################################################################
########https://www.r-bloggers.com/creating-annotated-data-frames-from-geo-with-the-geoquery-package/
###############################################################################

getGEOdataObjects <- function(x, getGSEobject=FALSE){
  # Make sure the GEOquery package is installed
  require("GEOquery")
  # Use the getGEO() function to download the GEO data for the id stored in x
  GSEDATA <- getGEO(x, GSEMatrix=T, AnnotGPL=FALSE)
  # Inspect the object by printing a summary of the expression values for the first 2 columns
  print(summary(exprs(GSEDATA[[1]])[, 1:2]))
  
  # Get the eset object
  eset <- GSEDATA[[1]]
  # Save the objects generated for future use in the current working directory
  save(GSEDATA, eset, file=paste(x, ".RData", sep=""))
  
  # check whether we want to return the list object we downloaded on GEO or
  # just the eset object with the getGSEobject argument
  if(getGSEobject) return(GSEDATA) else return(eset)
}

# Store the dataset ids in a vector GEO_DATASETS just in case you want to loop through several GEO ids
GEO_DATASETS <- c("GSE58661")

# Use the function we created to return the eset object
eset <- getGEOdataObjects(GEO_DATASETS[1])
# Inspect the eset object to get the annotation GPL id
eset 

# Get the annotation GPL id (see Annotation: GPL10558)
gpl <- getGEO('GPL15048', destdir=".")
Meta(gpl)$title

# Inspect the table of the gpl annotation object
colnames(Table(gpl))

# Get the gene symbol and entrez ids to be used for annotations
Table(gpl)[1:10, c(1, 6, 9, 12)]
dim(Table(gpl))

# Get the gene expression data for all the probes with a gene symbol
geneProbes <- which(!is.na(Table(gpl)$GeneSymbol))
probeids <- as.character(Table(gpl)$ID[geneProbes])

probes <- intersect(probeids, rownames(exprs(eset)))
length(probes)

geneMatrix <- exprs(eset)[probes, ]

inds <- which(Table(gpl)$ID %in% probes)
# Check you get the same probes
head(probes)
head(as.character(Table(gpl)$ID[inds]))

# Get gene symbol
Genesymbol <- Table(gpl)[inds, 4]

# Create the expression matrix with gene ids
geneMatTable <- cbind(as.character(Genesymbol), geneMatrix)
head(geneMatTable)

# Save a copy of the expression matrix as a csv file
write.csv(geneMatTable, paste(GEO_DATASETS[1], "_DataMatrix.csv", sep=""), row.names=T)


######## https://www.bioconductor.org/packages/3.3/bioc/vignettes/GEOquery/inst/doc/GEOquery.html#series
gse<- getGEO('GSE58661')

gse[[1]]
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform})
head(gsmplatforms)

gsmlist = Filter(function(gsm) {Meta(gsm)$platform=='GPL6098'},GSMList(gse))
length(gsmlist)

Table(gsmlist[[1]])[1:5,]
Columns(gsmlist[[1]])[1:5,]
Table(GPLList(gse)[[1]])[2]

probesets <- Table(GPLList(gse)[[1]])[2]
# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
{tab <- Table(x)
mymatch <- match(probesets,tab$ID_REF)
return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})

require(Biobase)
# go through the necessary steps to make a compliant ExpressionSet
rownames(data.matrix) <- probesets$Search_key
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
eset2

#data.matrix <- log2(data.matrix)
data.matrix[1:5,]
dim(data.matrix)
data.matrix <- data.matrix[complete.cases(data.matrix),]
data_GSEA <- c()
data_GSEA <- cbind(data.matrix[, colnames(data.matrix) %in% ID_Chr3_monosomy],
                   data.matrix[, colnames(data.matrix) %in% ID_Chr3_disomy])

sum(colnames(data.matrix) %in% ID_Chr3_monosomy)
sum(colnames(data.matrix) %in% ID_Chr3_disomy)
GSEA <- data.frame(NAME = rownames(data.matrix), 
                   DESCRIPTION = rep('na', dim(data_GSEA)[1]), data_GSEA)

write.table(GSEA, file = 'GSEA_GSE39717.txt', 
            quote = FALSE, row.names = FALSE, sep = '\t')

cls <- c(rep(0, sum(colnames(data.matrix) %in% ID_Chr3_monosomy)),
         rep(1, sum(colnames(data.matrix) %in% ID_Chr3_disomy)))

write(cls, file = 'cls_GSE39717.cls', append = FALSE, 
      #quote = FALSE, 
      sep = '\t',
      #row.names = FALSE,
      #col.names = TRUE, 
      ncolumns = length(cls)
)


data_P1 <- read.csv('Supplementary_Table_S1_Worley_et_al.csv', header = TRUE,
                    na.strings = c('NA','ND'))

## ND = not done, NA = not applicable

summary(data_P1)

sample <- read.csv('sample_GSE39717.csv', header = TRUE)
sample <- subset(sample, Collection == 'Fresh frozen sample from enucleated eye')
sample$ID <- gsub('MM ', '', sample$Title)
sample$ID

clin <- merge(data_P1, sample, by.x = 'MM', by.y = 'ID')
summary(clin)
summary(sample)

ID_Chr3_disomy <- clin$Accession[clin$Chr3_status_aCGH == 'Disomy']
ID_Chr3_monosomy <- clin$Accession[clin$Chr3_status_aCGH == 'Monosomy']
