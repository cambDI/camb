
caretEnsemble <- function(all.models, optFUN=NULL, ...){
  
  #TODO: Add progressbar argument and move optFUN to an all.models control argument
  
  #Libraries
  #require('caret') || stop("Pacakge 'caret' is required")
  #require('pbapply') || stop("Pacakge 'pbapply' is required")

  #Check the models, and make a matrix of obs and preds
  predobs <- makePredObsMatrix(all.models)

  #If the optimization function is NULL, choose default
  if (is.null(optFUN)){
    if (predobs$type=='Classification') {
      optFUN <- greedOptAUC
    } else { optFUN <- greedOptRMSE }
  }
  
  #Determine weights
  weights <- optFUN(predobs$preds, predobs$obs, ...)
  weights[! is.finite(weights)] <- 0
  
  #Normalize and name weights
  weights <- weights/sum(weights)
  names(weights) <- sapply(all.models, function(x) x$method)
  
  #Remove 0-weighted models
  keep <- which(weights != 0)
  
  #Determine RMSE
  if (predobs$type == "Regression"){
    error <- RMSE(predobs$preds %*% weights, predobs$obs)
    names(error) <- 'RMSE'
  } else {
    metric <- 'AUC'
    error <- colAUC(predobs$preds %*% weights, predobs$obs)
    names(error) <- 'AUC'
  }

  #Return final model
  out <- list(models=all.models[keep], weights=weights[keep], error=error)
  class(out) <- 'caretEnsemble'
  return(out) 
}

predict.caretEnsemble <- function(ensemble, ...){
  type <- checkModels_extractTypes(ensemble$models)
  preds <- multiPredict(ensemble$models, type, ...)
  out <- as.numeric(preds %*% ensemble$weights)
  return(out)
}

