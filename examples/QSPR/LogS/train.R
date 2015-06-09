library(caret)
library(kernlab)
library(doMC)

expGrid <- function(power.from, power.to, power.by, base){
  grid <- c()
  for (i in seq(power.from, power.to, power.by)){
    grid <- append(grid,base^i)
  }
  return(grid)
}

registerDoMC(cores=10)
dataset <- readRDS("dataset.rds")
method <- "rf"
tune.grid <- expand.grid(.mtry = seq(5,100,5))
model <- train(dataset$x.train, dataset$y.train, method, tuneGrid=tune.grid, trControl=dataset$trControl)
saveRDS(model, file=paste(method,".rds",sep=""))

method <- "gbm"
tune.grid <- expand.grid(.n.trees=c(500,1000), .interaction.depth=c(25), .shrinkage = c(0.04, 0.08, 0.16))
model <- train(dataset$x.train, dataset$y.train, method, tuneGrid=tune.grid, trControl=dataset$trControl)
saveRDS(model, file=paste(method,".rds",sep=""))