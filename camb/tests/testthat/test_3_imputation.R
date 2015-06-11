context("Testing Imputation")

test_that("ReplaceInfinitesWithNA and ImputeFeatures provides outputs consistent with reference outputs", {
  #properties <- read.table("properties.csv", header=TRUE, sep="\t")
  #properties <- properties[properties$Kept==1, ]
  #head(properties)
  #targets <- data.frame(Name = properties$NAME, target = properties$EXPT)
  
  #all <- merge(x=targets, y=descriptors, by="Name")
  #ids <- all$Name
  #x <- all[3:ncol(all)]
  #y <- all$target
  
  #x.finite <- ReplaceInfinitesWithNA(x)
  #x.imputed <- ImputeFeatures(x.finite)
})