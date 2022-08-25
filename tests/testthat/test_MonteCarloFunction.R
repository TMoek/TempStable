test_that("TemperedEstim_Simulation_with_Classic_ML_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(1.5,1,1,1,1,0),c(0.5,1,1,1,1,0)),
      SampleSizes = 2,
      MCparam = 2,
      TemperedType = "Classic",
      Estimfct = "ML",
      saveOutput = FALSE)


    expect_equal(length(TestObject$outputMat),64)
    expect_equal(TestObject$outputMat[1],1.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "delta+T",
                                                  "delta-T", "lambda+T",
                                                  "lambda-T", "muT",
                                                  "data size", "seed",
                                                  "alphaE", "delta+E",
                                                  "delta-E", "lambda+E",
                                                  "lambda-E", "muE", "failure",
                                                  "time"))

  })

  #expect_error()

})



