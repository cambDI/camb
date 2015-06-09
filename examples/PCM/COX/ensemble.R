library(devtools) || stop("package devtools need to be installed")
library(pbapply) || stop("package pbapply need to be installed")


source("setup_ensemble.R")

all.models <- list()

models <- as.vector(read.table("models_ensemble")$V1)

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
linear <- caretStack(all.models, method='glm')
summary(linear$ens_model$finalModel)

# make Elastic Net ensemble
enet_ens <- caretStack(all.models, method='enet')
coefs_enet_ens <- enet_ens$ens_model$finalModel$beta.pure[ncol(enet_ens$ens_model$finalModel$beta.pure)+1,]

# make SVM linear ensemble
trControl <- trainControl(method = "cv",  number=5)
tune.grid <- expand.grid(.C=expGrid(power.from=-14,power.to=10,power.by=1,base=2))
linear_svm <- caretStack(all.models, method='svmLinear',trControl=trControl,tuneGrid=tune.grid)

# make SVM radial ensemble
trControl <- trainControl(method = "cv",  number=5)
tune.grid <- expand.grid(.sigma=expGrid(power.from=-14,power.to=10,power.by=1,base=2),
                         .C=expGrid(power.from=-14,power.to=10,power.by=2,base=2))
radial_svm <- caretStack(all.models, method='svmRadial',trControl=trControl,tuneGrid=tune.grid)


# Predict the bioactivities for test set:
preds <- data.frame(sapply(all.models, predict, newdata=dataa$x.test))
preds$ENS_greedy <- predict(greedy, newdata=dataa$x.test)
preds$ENS_linear <- predict(linear, newdata=dataa$x.test)
preds$ENS_enet <- predict(enet_ens, newdata=x.test)
preds$ENS_SVMrad <- predict(radial_svm, newdata=x.test)c
preds$ENS_SVMlin <- predict(linear_svm, newdata=x.test)


# Calculate metrics (We could also apply Validation instead.)
Q2s = apply(preds,2, function(x) Qsquared(x,dataa$y.test))
R2s <- apply(preds,2, function(x) Rsquared(x,dataa$y.test))
R20s <- apply(preds,2, function(x) Rsquared0(x,dataa$y.test))
RMSEs <- apply(preds,2, function(x) RMSE(x,dataa$y.test))

# We sort the models by RMSE
sort(RMSEs)
