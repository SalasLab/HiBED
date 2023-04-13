test_that("errors if bad parameters", {
  library(HiBED)
  library(FlowSorted.DLPFC.450k)
  data(FlowSorted.DLPFC.450k)
  expect_error(HiBED_deconvolution(FlowSorted.DLPFC.450k,
                                   h=2
  ))
  Mset<-minfi::preprocessRaw(FlowSorted.DLPFC.450k)
  Examples_Betas<-minfi::getBeta(Mset)
  expect_error(HiBED_deconvolution(Examples_Betas,
                                   h=2,"brain"
  ))
  expect_error(HiBED_deconvolution(Mset,
                                   h=1
  ))
})
