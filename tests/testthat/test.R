test_that("errors if bad parameters", {
  library(HiBED)
  library(FlowSorted.DLPFC.450k)
  library(minfi)
  data(FlowSorted.DLPFC.450k)
  expect_error(HiBED_deconvolution(FlowSorted.DLPFC.450k,
                                   h=2
  ))
})
