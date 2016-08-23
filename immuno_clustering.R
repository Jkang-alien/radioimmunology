setwd('/home/jun/radioimmunology/')

immune <- read.delim('ssGSEA_GSE18865.txt', sep = '\t')


cell_type <- as.character(immune$Description)
ID <- colnames(immune)[2:90]

ES <- t(immune[,2:90])
ES <- scale(ES, center = TRUE, scale = TRUE)
#ES <- sweep(ES, 1, apply(ES, 1, median, na.rm = TRUE))
colnames(ES) <- cell_type

GEO_ID <- read.delim('GEO_ID.txt')
rownames(ES) == GEO_ID$ID_chip
rownames(ES) <- GEO_ID$ID_lung

#################### Patient info ##############################
clin <- read.delim('patient.txt')
summary(clin)
clin$Contrast..0.non..1.cont <- factor(clin$Contrast..0.non..1.cont,
                                       levels = c(0, 1),
                                       labels = c('pre', 'post'))

rownames(ES)
clin$PatientID <- sub('3-[0]*', '_', tolower (clin$PatientID))

rownames(ES) == clin$PatientID

clin$Histology_simple <- rep(NA, 89)
clin$Histology_simple [grep('Adeno', clin$Histology)] <- 'Adenocarcinoma'
clin$Histology_simple [grep('Solid', clin$Histology)] <- 'Adenocarcinoma'
clin$Histology_simple [grep('Squamous', clin$Histology)] <- 'Squamous Cell Carcinoma'
clin$Histology_simple [grep('Non-Small Cell', clin$Histology)] <- 'Non-Small Cell Carcinoma'
clin$Histology_simple [grep('Large Cell', clin$Histology)] <- 'Large Cell Neuroendocrine Carcinoma'

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

ann <- data.frame(case_class = factor(results_row[[4]]$consensusClass),
                  Histologic_type = clin$Histology_simple
                  )
ann_col <- data.frame(immune_class = factor(results_col[[3]]$consensusClass))

CairoPDF(file = 'immune_TCIA_Lung.pdf',
         width =7.5, height = 7.5, pointsize = 16)
ah<- aheatmap(ES, 
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

################## Radiomics ########################
radio <- read.delim('lung_3_radio.txt')

summary(radio)
mat_radio <- data.matrix (radio[,4:164])

############ ML ########################

library(mlbench)
library(caret)
correlationMatrix <- cor(mat_radio,use="complete.obs")
# summarize the correlation matrix
heatmap(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
# print indexes of highly correlated attributes
print(highlyCorrelated)

dim(mat_radio[, -highlyCorrelated])

control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(diabetes~., data=PimaIndiansDiabetes, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)

library(rpart)

data_rpart_post <- data.frame(class = factor(results_row[[4]]$consensusClass),
                                        mat_radio[, -highlyCorrelated])[clin$Contrast..0.non..1.cont == 'post',]

rpart_post <- rpart(class~., data = data_rpart_post,
                    control = rpart.control(minsplit = 10))

attributes(rpart_post)
plot(rpart_post)
text(rpart_post, use.n=T,
     cex = 0.5)
rpart_post$ordered

library(party)
ctree_post <- ctree(class~., data = data_rpart_post,
                    controls = ctree_control(maxsurrogate = 4))

table(predict(ctree_post), data_rpart_post$class)
plot(ctree_post)

library(randomForest)
rf <- randomForest(class~., data = data_rpart_post, ntree=100, proximity=TRUE)
table(predict(rf), data_rpart_post$class)

sink('histology_enhance.txt')
table(clin$Contrast..0.non..1.cont, clin$Histology_simple)
sink()
  
############################################################
################# cor p-value ##############################
###https://www.r-bloggers.com/more-on-exploring-correlations-in-r/

cor.prob <- function (X, Y, vp = c()) {
  for (i in 1:dim(X)[2]){
    for (j in 1:dim(Y)[2]){
  a <- cor.test(X[,i],Y[,j])
  p <- a$p.value
  vp <- append(vp, p)
    }
  }
matrix(vp, nrow = dim(X)[2], byrow = TRUE,
       dimnames = list(colnames(X), colnames(Y)))
}

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}

index_adeno_post <- clin$Contrast..0.non..1.cont == 'post' & clin$Histology_simple == 'Adenocarcinoma'

cor_radio_immue <- cor(mat_radio[index_adeno_post,],
                       ES[index_adeno_post,],
                       use="complete.obs")
ind_sig <- which(abs(cor_radio_immue) > 0.5, arr.ind = TRUE)
colnames(cor_radio_immue)[ind_sig[,2]]

aheatmap(cor_radio_immue)

index_adeno_post <- clin$Contrast..0.non..1.cont == 'post' & clin$Histology_simple == 'Adenocarcinoma'

cor_radio_immue_sig <- cor.prob(mat_radio[index_adeno_post,],
                       ES[index_adeno_post,])

aheatmap(log(cor_radio_immue_sig))
ind_sig <- which(log(cor_radio_immue_sig) < -6, arr.ind = TRUE)
colnames(cor_radio_immue_sig)[ind_sig[,2]]

plot(mat_radio[index_adeno_post,60], ES[index_adeno_post,3])
aheatmap(cor_radio_immue)