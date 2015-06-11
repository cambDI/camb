context("Testing Compound Standardisation")

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
                                      max.mass.limit=900,
                                      number.processed=200)
  
  expect_equal(readRDS("reference_standardisation_options.rds"), std.options)
  expect_equal(readChar("reference_standardisation.log", file.info("reference_standardisation.log")$size), 
               readChar("standardisation.log", file.info("standardisation.log")$size))
  expect_equal(read.csv("reference_properties.csv"), read.csv("properties.csv"))

  reference_lines <- readLines("reference_removed.sdf")
  reference_lines <- reference_lines[-grep("INDIGO", reference_lines)]
  lines <- readLines("removed.sdf")
  lines <- lines[-grep("INDIGO", lines)]                                
  expect_equal(reference_lines, lines)
  
  reference_lines <- readLines("reference_standardised.sdf")
  reference_lines <- reference_lines[-grep("INDIGO", reference_lines)]
  lines <- readLines("standardised.sdf")
  lines <- lines[-grep("INDIGO", lines)]                                
  expect_equal(reference_lines, lines)
})

?StandardiseMolecules


