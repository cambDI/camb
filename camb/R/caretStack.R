
caretStack <- function(all.models, ...){
  #Libraries
  #require('caret')
  
  #Check the models, and make a matrix of obs and preds
  predobs <- makePredObsMatrix(all.models)
  
  #Build a caret model
  model <- train(predobs$preds, predobs$obs, ...)
  
  #Return final model
  out <- list(models=all.models, ens_model=model, error=model$results)
  class(out) <- 'caretStack'
  return(out) 
}

predict.caretStack <- function(ensemble, newdata=NULL, ...){
  type <- checkModels_extractTypes(ensemble$models)
  preds <- multiPredict(ensemble$models, newdata=newdata, type)
  out <- predict(ensemble$ens_model, newdata=preds, ...)
  return(out)
}

