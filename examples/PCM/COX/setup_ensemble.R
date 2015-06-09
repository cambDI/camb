#!/c5/shared/R/3.0.2/bin/Rscript
library(caret)

load("../GBM.RData")

dataa <- list()
dataa$transformation <- transformation
dataa$x.train <- x.train
dataa$y.train <- y.train
dataa$x.test <- x.test
dataa$y.test <- y.test

set.seed(1)
folds=5
repeats=1
trControl <- trainControl(method='cv', number=folds, repeats=repeats, returnResamp='none',
                          returnData=FALSE, savePredictions=TRUE,
                          verboseIter=TRUE, allowParallel=TRUE,
                          index=createMultiFolds(y.train, k=folds, times=repeats))

dataa$trControl <- trControl
saveRDS(dataa, file="data_ensemble.rds")
