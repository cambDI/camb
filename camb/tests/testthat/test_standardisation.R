context("Testing Compound Standardisation")

library(camb)
library(tools) 

parseToSame <- function(file1, file2) {
  a <- parse(file = file1)
  b <- parse(file = file2)
  attributes(a) <- NULL
  attributes(b) <- NULL
  identical(a,b)
}

test_that("StandardiseMolecules provides outputs consistent with reference outputs", {
  std.options <- StandardiseMolecules(structures.file="solubility_2007_ref2.sdf", 
                                      standardised.file="standardised.sdf", 
                                      removed.file="removed.sdf",
                                      properties.file = "properties.csv",
                                      remove.inorganic=TRUE, 
                                      fluorine.limit=3, 
                                      chlorine.limit=3, 
                                      bromine.limit=3, 
                                      iodine.limit=3, 
                                      min.mass.limit=20, 
                                      max.mass.limit=900)
  
  expect_equal(readRDS("reference_standardisation_options.rds"), std.options)
  expect_equal(readChar("reference_standardisation.log", file.info("reference_standardisation.log")$size), 
               readChar("standardisation.log", file.info("standardisation.log")$size))
  expect_equal(read.csv("reference_properties.csv"), read.csv("properties.csv"))
  #expect_equal(readChar("reference_removed.sdf", file.info("reference_removed.sdf")$size), # need to fix these
  #             readChar("removed.sdf", file.info("removed.sdf")$size))
  #expect_equal(readChar("reference_standardised.sdf", file.info("reference_standardised.sdf")$size), 
  #             readChar("standardised.sdf", file.info("standardised.sdf")$size))
})



