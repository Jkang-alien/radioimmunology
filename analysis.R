########### LASSO ####################
library(glmnet)
radio_remove_na <- mat_radio[,complete.cases(t(mat_radio))] 
radio_remve_inf <- radio_remove_na[, apply(radio_remove_na, 2, function(x) all(is.finite(x)))] 


################

f_lasso <- function (cell){
  data <- data.frame(cell = ES[,colnames(ES)==cell], radio_remve_inf)
  set.seed(1)
  train = sample(c(TRUE, FALSE), nrow(data), rep = TRUE)
  test = (!train)
  x = model.matrix(cell~., data)[,-1]
  y = data$cell
  grid=10^seq(10, -2, length = 100)
  lasso.mod = glmnet(x, y, alpha=1, lambda = grid)
  #plot(lasso.mod)
  
  cv.out<-cv.glmnet(x[train,],y[train],alpha=1)
  #plot(cv.out)
  bestlam<-cv.out$lambda.min
  
  lasso.pred<-predict(lasso.mod,s=bestlam,newx=x[test,])
  RSS <- mean((lasso.pred-y[test])^2)
  return(RSS)
  
}

RSS <- c()
for (i in colnames(ES)){
  a <- f_lasso(i)
  RSS <- append(RSS, a)
}

f_lasso(colnames(ES)[9])

colnames(ES)[23]


######################################################

data <- data.frame(Th2 = ES[,23], radio_remve_inf)
set.seed(1)
train = sample(c(TRUE, FALSE), nrow(data), rep = TRUE)
test = (!train)
x = model.matrix(Th2~., data)[,-1]
y = data$Th2
grid=10^seq(10, -2, length = 100)
lasso.mod = glmnet(x, y, alpha=1, lambda = grid)
plot(lasso.mod)

cv.out<-cv.glmnet(x[train,],y[train],alpha=1)
plot(cv.out)
bestlam<-cv.out$lambda.min

lasso.pred<-predict(lasso.mod,s=bestlam,newx=x[test,])
mean((lasso.pred-y[test])^2)

#using the full data
model<-glmnet(x,y,alpha=1,lambda=grid)
lasso.coef<-predict(model,type="coefficients",s=bestlam)
lasso.coef[lasso.coef!=0]

plot(data$Th2[test], lasso.pred)

#########################################################
############ Best sybset ################################
library(leaps)
regfit.best <- regsubsets(Th2~., data = data[train,], nvmax = 10)
test.mat <- model.matrix(Th2~., data[test,])
val.errors <- rep(NA,3)
for (i in 1:10){
  coefi = coef(regfit.best, id = i)
  pred = test.mat[, names(coefi)]%*%coefi
  val.errors[i] = mean((data$Th2[test]-pred)^2)
}

val.errors
which.min(val.errors)
coef(regfit.best, 2)
