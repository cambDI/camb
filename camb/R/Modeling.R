#################################################################################
## Modeling 
#################################################################################
asNumeric <- function(x) as.numeric(as.character(x))
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)],   
                                                   asNumeric))

ReplaceInfinitesWithNA <- function(d) {
  do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x), NA)))
}
                         
ImputeFeatures <- function(d, k=10,...) {
 if (!("impute" %in% rownames(installed.packages()))){
 cat("Installing package impute...\n")
 source("http://bioconductor.org/biocLite.R")
 biocLite("impute")
 library(impute)
 }
  as.data.frame(impute.knn(as.matrix(factorsNumeric(d)), k = k, ...)$data)
} 

SplitSet <- function(ids, x, y, percentage = 20, seed = 1) {
  holdout.size <- round(nrow(x) * (percentage/100))
  set.seed(seed)
  holdout.indexes <- sample(1:nrow(x), holdout.size, replace=FALSE)
  train.indexes <- (1:length(y))[-holdout.indexes]
  x.train <- x[train.indexes, ]
  x.holdout <- x[holdout.indexes, ]
  y.train <- y[train.indexes]
  y.holdout <- y[holdout.indexes]
  l <- list()
  l$ids <- ids
  l$holdout.indexes <- holdout.indexes
  l$train.indexes <- train.indexes
  l$x.train <- x.train
  l$x.holdout <- x.holdout
  l$y.train <- y.train
  l$y.holdout <- y.holdout
  l
}

RemoveNearZeroVarianceFeatures <- function(dataset, frequencyCutoff = 30/1,...) {
  nzv.columns <- nearZeroVar(dataset$x.train, freqCut = frequencyCutoff,...)
  if (length(nzv.columns) != 0) {
    message(paste(length(nzv.columns), "features removed with variance below cutoff"))
    dataset$x.train <- dataset$x.train[, -nzv.columns]
    dataset$x.holdout <- dataset$x.holdout[, -nzv.columns]
  }
  else {
    message("no features removed")
  }
  dataset
}

RemoveHighlyCorrelatedFeatures <- function(dataset, correlationCutoff = 0.95,...) {
  hc.columns <- findCorrelation(cor(dataset$x.train), correlationCutoff,...)
  if (length(hc.columns) != 0) {
    message(paste(length(hc.columns), "features removed with correlation above cutoff"))
    dataset$x.train <- dataset$x.train[, -hc.columns]
    dataset$x.holdout <- dataset$x.holdout[, -hc.columns]
  }
  else {
    message("no features removed")
  }
  dataset
}

PreProcess <- function(dataset, steps = c("center", "scale"),...) {
  transformation <- preProcess(dataset$x.train, method = steps,...)
  dataset$x.train <- predict(transformation, dataset$x.train)
  dataset$x.holdout <- predict(transformation, dataset$x.holdout)
  dataset$transformation <- transformation
  dataset
}

GetCVTrainControl <- function(dataset, seed = 1, folds = 5, repeats = 1, method='cv', returnResamp='none', returnData=TRUE, savePredictions=TRUE, verboseIter=TRUE, allowParallel=TRUE,...) {
  set.seed(seed)
  dataset$trControl <- trainControl(method='cv', number=folds, repeats=repeats, returnResamp='none',
                                    returnData=TRUE, savePredictions=TRUE,
                                    verboseIter=TRUE, allowParallel=TRUE,
                                    index=createMultiFolds(dataset$y.train, k=folds, times=repeats),...)
  dataset
}

expGrid <- function(power.from, power.to, power.by, base){
  grid <- c()
  for (i in seq(power.from, power.to, power.by)){
    grid <- append(grid,base^i)
  }
  return(grid)
}

YScrambling <- function(y,percent){
if (percent < 0 || percent > 1){stop("The percent value needs to be between 0 and 1")}
inds_to_resamp <- sample.int(length(y), length(y)*percent)
resamp_vals <- y[inds_to_resamp]
resamp_vals2 <- sample(resamp_vals,length(resamp_vals),replace = FALSE, prob = NULL)
y[inds_to_resamp] <- resamp_vals2
return(y)
}


