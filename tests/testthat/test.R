test_that("errors if bad parameters", {
  library(HiBED)
  library(FlowSorted.DLPFC.450k)
  library(minfi)
  data(FlowSorted.DLPFC.450k)
  expect_error(estimateCellCounts2(FlowSorted.DLPFC.450k,
                                   h=2
  ), "object is of class 'RGChannelSet', but needs to be a processed Beta matrix")
})
