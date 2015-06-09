#!/c5/shared/R/3.0.2/bin/Rscript


#install.packages("devtools")
library(devtools)
library(caret)
library(pbapply)
#install_github('caretEnsemble', 'zachmayer')
#install.packages("/Bis/home/icortes/Desktop/PCM/DHFR/ensemble/caretEmsemble",repos=NULL)
library(caretEnsemble)
library(gbm)
source("setup.R")


#dataa < readRDS("dataa.rds")
typeof(dataa)

all.models <- list()

models <- as.vector(read.table("models_ensemble_best")$V1)

for (i in 1:length(models)){

model_load = paste("readRDS('",models[i],"')",sep="")
assign(paste("model_",i,sep=""), eval(parse(text=model_load)))

all.models[[length(all.models)+1]] <- eval(parse(text=paste("model_",i,sep="")))

}

names(all.models) <- sapply(all.models, function(x) x$method)
sort(sapply(all.models, function(x) min(as.vector(na.omit(x$results$RMSE)))))

# make a greedy ensemble - currently can only use RMSE
greedy <- caretEnsemble(all.models, iter=1000L)
sort(greedy$weights, decreasing=TRUE)

# make a linear regression ensemble
linear <- caretStack(all.models, method='glm')#, trControl=dataa$trainControl(method='cv'))
summary(linear$ens_model$finalModel)
linear$error

# make a rf regression ensemble
#tune.grid <- expand.grid(.mtry = seq(1,length(all.models),1))
#nonlinear <- caretStack(all.models, method='rf', trControl=trainControl(method='cv'), tune.grid=tune.grid)
#summary(nonlinear$ens_model$finalModel)
#nonlinear$error

#Predict for test set:
preds <- data.frame(sapply(all.models, predict, newdata=dataa$x.test))
preds$ENS_greedy <- predict(greedy, newdata=dataa$x.test)
preds$ENS_linear <- predict(linear, newdata=dataa$x.test)
#preds$ENS_nonlinear <- predict(nonlinear, newdata=x.test)
sort(sqrt(colMeans((preds - y.test) ^ 2)))

# Calculate metrics
source("../functions.R")
print("Q2")
Q2s = apply(preds,2, function(x) Qsquared(x,dataa$y.test))
print("R2")
R2s <- apply(preds,2, function(x) Rsquared(x,dataa$y.test))
print("R20")
R20s <- apply(preds,2, function(x) Rsquared0(x,dataa$y.test))
print("RMSE")
RMSEs <- apply(preds,2, function(x) RMSE(x,dataa$y.test))

RMSEs
Q2s

save(greedy, file="greedy_best.RData")
save(linear, file="linear_best.RData")

save(list=ls(), file="ensemble_COX_best.RData")
