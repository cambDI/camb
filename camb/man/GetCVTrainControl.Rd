\name{GetCVTrainControl}
\alias{GetCVTrainControl}
\title{Sets up the control parameters of caret's \code{\link{train}} function.}
\usage{
GetCVTrainControl(dataset, seed = 1, folds = 5, repeats = 1,
  method = "cv", returnResamp = "none", returnData = TRUE,
  savePredictions = TRUE, verboseIter = TRUE, allowParallel = TRUE, ...)
}
\arguments{
\item{dataset}{The training dataset returned by \code{\link{SplitSet}}}

\item{seed}{The seed for randomization so that the fold selection can be done in a repeatable way if desired}

\item{folds}{The number of folds to use during cross-validation}

\item{repeats}{For repeated k-fold cross-validation only: the number of complete sets of folds to compute}

\item{returnResamp}{A character string indicating how much of the resampled summary metrics should be saved.
Values can be final, all or none}

\item{returnData}{A logical for saving the data}

\item{savePredictions}{A logical to save the hold-out predictions for each resample}

\item{verboseIter}{A logical for printing a training log.}

\item{allowParallel}{If a parallel backend is loaded and available, should the function use it?}
}
\value{
a dataset with the traincontrol saved within for future training
}
\description{
Calls caret's \code{\link{trainControl}} function to set up the parameters of the \code{\link{train}} function. This control variable is saved
into the \code{dataset} list as \code{dataset$trControl}.
}
\author{
Daniel Murrell <dsmurrell@gmail.com> and Isidro Cortes <isidrolauscher@gmail.com>
}

