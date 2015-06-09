
PredictExternal <- function(structures.file, standardisation.options, descriptor.types, dataset, model) {
  model = model # this is done to make sure model gets evaluated successfully first
  
  standardised.file <- tempfile("standardised", fileext=".sdf")
  removed.file <- tempfile("removed", fileext=".sdf")
  properties.file <- tempfile("properties", fileext=".csv")
  
  StandardiseMolecules(structures.file = structures.file, 
                       standardised.file = standardised.file, 
                       removed.file = removed.file,
                       properties.file = properties.file,
                       standardisation.options)
  
  descriptors <- GeneratePadelDescriptors(standardised.file = standardised.file, types=c("2D"), threads = 1)
  descriptors <- RemoveStandardisedPrefix(descriptors)
  
  ids <- descriptors[,1]
  x <- descriptors[,2:ncol(descriptors)]
  
  x <- x[, names(dataset$x.train)]
  x.finite <- ReplaceInfinitesWithNA(x)
  x.imputed <- ImputeFeatures(x.finite)
  
  x.preprocessed <- predict(dataset$transformation, x.imputed)
  predictions <- predict(model, newdata=x.preprocessed)
  return(data.frame(id=ids, prediction=predictions))
}
