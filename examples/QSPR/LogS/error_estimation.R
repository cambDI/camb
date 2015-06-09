#########################################
# Error estimation
#########################################
dataset <- readRDS("dataset.rds")
model <- readRDS("rf.rds")
library(errorestimatoR)
library(kernlab)
estimator4 = BuildCaretErrorEstimator(dataset$x.train, model, Nmax=20, cores=1)

sigpred4 <- PredictSigmas(x=dataset$x.holdout, estimator4)

preds <- predict(model, newdata = dataset$x.holdout)
errors <- preds - dataset$y.holdout

plot(errors, sigpred4$sigmas)
plot(errors, sigpred$sigma.matrix[,2])
x <- ss$x.test
CEC(sigpred4$sigmas, abs(errors))

CEC(sigpred4$sigma.matrix[,2], abs(errors))

plot(density(apply(sigpred4$sigma.matrix, 2, CEC, abs(errors)), bw=0.01))

plot(apply(sigpred4$sigma.matrix, 2, CEC, abs(errors)))

k <- rbfdot(sigma = 0.2)
k

?train

trainsig <- PredictSigmas(x=ss$x.train, estimator)
trainpreds <- predict(model, newdata = ss$x.train)
trainerrors <- trainpreds - ss$y.train

plot(trainerrors, trainsig$sigmas)

CEC(abs(trainerrors), trainsig$sigmas)

CEC(abs(trainerrors), trainsig$sigma.matrix[,2])

length(abs(trainerrors))
length(trainsig$sigma.matrix[,2])

dim(estimator$x)
dim(x)

dm <- CreateNewDistanceMatrix(estimator$x, x)
nm <- CreateNewEstimatorMatrix(estimator$Nmax, dm, estimator$errors, 
                               estimator$obs, estimator$preds)
sigma.matrix <- GetNewSigmaMatrix(nm, estimator$sl)
sigmas <- GetNewSigmas(sigma.matrix, estimator$weights)
prediction <- list()
prediction$dm <- dm
prediction$nm <- nm
prediction$sigma.matrix <- sigma.matrix
prediction$sigmas <- sigmas
prediction