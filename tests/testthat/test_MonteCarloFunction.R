
### CTS ###
# For ML, warnings will occur as NaNs will be produced --> suppressWarning()



test_that("TemperedEstim_Simulation_with_CTS_ML_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(1.5,1,1,1,1,0),c(0.5,1,1,1,1,0)),
      SampleSizes = 2,
      MCparam = 2,
      TemperedType = "CTS",
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
})


test_that("TemperedEstim_Simulation_with_CTS_GMM_gives_correct_return", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
    SampleSizes = 4, MCparam = 2, TemperedType = "CTS",
    Estimfct = "GMM", saveOutput = FALSE, algo = "2SGMM",
    regularization = "cut-off", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


  expect_equal(length(TestObject$outputMat),32)
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

test_that("TemperedEstim_Simulation_with_CTS_Cgmm_gives_correct_return", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)), SampleSizes = 4,
    MCparam = 1, TemperedType = "CTS", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "2SCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)


  expect_equal(length(TestObject$outputMat), 16)
  expect_equal(TestObject$outputMat[1],1.45)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "delta+T",
                                                 "delta-T", "lambda+T",
                                                 "lambda-T", "muT",
                                                 "data size", "seed",
                                                 "alphaE", "delta+E",
                                                 "delta-E", "lambda+E",
                                                 "lambda-E", "muE", "failure",
                                                 "time"))
})

test_that("TemperedEstim_Simulation_with_CTS_GMC_gives_correct_return", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)), SampleSizes = 4,
    MCparam = 2, TemperedType = "CTS", Estimfct = "GMC",
    saveOutput = FALSE, algo = "2SGMC", alphaReg = 0.01,
    WeightingMatrix = "OptAsym", regularization = "cut-off", ncond = 8)


  expect_equal(length(TestObject$outputMat), 32)
  expect_equal(TestObject$outputMat[1],1.45)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "delta+T",
                                                 "delta-T", "lambda+T",
                                                 "lambda-T", "muT",
                                                 "data size", "seed",
                                                 "alphaE", "delta+E",
                                                 "delta-E", "lambda+E",
                                                 "lambda-E", "muE", "failure",
                                                 "time"))
})


### TSS ###

test_that("TemperedEstim_Simulation_with_TSS_ML_gives_correct_return",
          {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1)), SampleSizes = 4, MCparam = 1,
      TemperedType = "TSS", Estimfct = "ML", saveOutput = FALSE)


    expect_equal(length(TestObject$outputMat), 10)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                   "lambdaT", "data size",
                                                   "seed", "alphaE", "deltaE",
                                                   "lambdaE", "failure","time"))

  })
})

test_that("TemperedEstim_Simulation_with_TSS_GMM_gives_correct_return",
  {

    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1)), SampleSizes = 4, MCparam = 1,
      TemperedType = "TSS", Estimfct = "GMM", saveOutput = FALSE,
      algo = "2SGMM", regularization = "cut-off", WeightingMatrix = "OptAsym",
      t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


    expect_equal(length(TestObject$outputMat), 10)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                   "lambdaT", "data size",
                                                   "seed", "alphaE", "deltaE",
                                                   "lambdaE", "failure","time"))
})


test_that("TemperedEstim_Simulation_with_TSS_Cgmm_gives_correct_re", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.45,0.55,1)), SampleSizes = 4,
    MCparam = 1, TemperedType = "TSS", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "2SCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)


  expect_equal(length(TestObject$outputMat), 10)
  expect_equal(TestObject$outputMat[1],0.45)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                 "lambdaT", "data size",
                                                 "seed", "alphaE", "deltaE",
                                                 "lambdaE", "failure","time"))
})

test_that("TemperedEstim_Simulation_with_TSS_GMC_gives_correct_re",
          {

            TestObject <- TemperedEstim_Simulation(
              ParameterMatrix = rbind(c(0.45,0.55,1)), SampleSizes = 4,
              MCparam = 2, TemperedType = "TSS", Estimfct = "GMC",
              saveOutput = FALSE, algo = "2SGMC", alphaReg = 0.01,
              WeightingMatrix = "OptAsym", regularization = "cut-off",
              ncond = 8)


            expect_equal(length(TestObject$outputMat), 20)
            expect_equal(TestObject$outputMat[1],0.45)
            expect_equal(colnames(TestObject$outputMat),
                         c("alphaT", "deltaT", "lambdaT", "data size", "seed",
                           "alphaE", "deltaE", "lambdaE", "failure","time"))
          })


### NTS ###

test_that("TemperedEstim_Simulation_with_NTS_ML_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1,1,1)), SampleSizes = 4, MCparam = 1,
      TemperedType = "NTS", Estimfct = "ML", saveOutput = FALSE)


    expect_equal(length(TestObject$outputMat), 14)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT", "deltaT",
                                                   "lambdaT", "muT",
                                                   "data size", "seed",
                                                   "alphaE", "betaE", "deltaE",
                                                   "lambdaE", "muE", "failure",
                                                   "time"))

  })
})

test_that("TemperedEstim_Simulation_with_NTS_GMM_gives_correct_return", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.5,1,1,1,1)), SampleSizes = 4, MCparam = 2,
    TemperedType = "NTS", Estimfct = "GMM", saveOutput = FALSE,
    algo = "2SGMM", regularization = "cut-off", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


  expect_equal(length(TestObject$outputMat), 28)
  expect_equal(TestObject$outputMat[1],0.5)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT", "deltaT",
                                                 "lambdaT", "muT",
                                                 "data size", "seed",
                                                 "alphaE", "betaE", "deltaE",
                                                 "lambdaE", "muE", "failure",
                                                 "time"))
})

test_that("TemperedEstim_Simulation_with_NTS_Cgmm_gives_correct_return", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.55,0.55,1,1,1)), SampleSizes = 4,
    MCparam = 1, TemperedType = "NTS", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "2SCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)


  expect_equal(length(TestObject$outputMat), 14)
  expect_equal(TestObject$outputMat[1],0.55)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT", "deltaT",
                                                 "lambdaT", "muT",
                                                 "data size", "seed",
                                                 "alphaE", "betaE", "deltaE",
                                                 "lambdaE", "muE", "failure",
                                                 "time"))
})









