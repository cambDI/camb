context("Testing Descriptor Generation")

test_that("GeneratePadelDescriptors provides outputs consistent with reference outputs", {
  descriptor.types <- c("2D") 
  descriptors <- GeneratePadelDescriptors(standardised.file = "standardised.sdf", types=descriptor.types, threads = 1)
  descriptors <- RemoveStandardisedPrefix(descriptors)
  expect_equal(readRDS("reference_descriptors.rds"), descriptors)
})