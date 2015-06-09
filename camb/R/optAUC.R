greedOptAUC <- function(X, Y, iter = 100L){ 
  #require('caTools') || stop("Pacakge 'caTools' is required")

  if(is.character(Y)){
    Y <- factor(Y)
  }
  stopifnot(is.factor(Y))
  
  N           <- ncol(X)
  weights     <- rep(0L, N)
  pred        <- 0 * X
  sum.weights <- 0L
  
  while(sum.weights < iter) {
    
    sum.weights   <- sum.weights + 1L
    pred          <-(pred + X) * (1L / sum.weights)
    errors        <- colAUC(pred, Y)
    best          <- which.max(errors)
    weights[best] <- weights[best] + 1L
    pred          <- pred[, best] * sum.weights
  }
  return(weights)
}
