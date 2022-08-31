
### Classic ###

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
})


test_that("TemperedEstim_Simulation_with_Classic_GMM_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
      SampleSizes = 40, MCparam = 40, TemperedType = "Classic",
      Estimfct = "GMM", saveOutput = FALSE, algo = "2SGMM",
      regularization = "cut-off", WeightingMatrix = "OptAsym",
      t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


    expect_equal(length(TestObject$outputMat),640)
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

test_that("TemperedEstim_Simulation_with_Classic_Cgmm_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)), SampleSizes = 4,
      MCparam = 4, TemperedType = "Classic", Estimfct = "Cgmm",
      saveOutput = FALSE, algo = "2SCgmm", alphaReg = 0.01, subdivisions = 20,
      IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
      s_max= 1)


    expect_equal(length(TestObject$outputMat), 64)
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
})

test_that("TemperedEstim_Simulation_with_Classic_GMC_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)), SampleSizes = 4,
      MCparam = 4, TemperedType = "Classic", Estimfct = "GMC",
      saveOutput = FALSE, algo = "2SGMC", alphaReg = 0.01,
      WeightingMatrix = "OptAsym", regularization = "cut-off", ncond = 8)


    expect_equal(length(TestObject$outputMat), 64)
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
})


### Subordinator ###

test_that("TemperedEstim_Simulation_with_Subordinator_ML_gives_correct_return",
          {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1)), SampleSizes = 10, MCparam = 10,
      TemperedType = "Subordinator", Estimfct = "ML", saveOutput = FALSE)


    expect_equal(length(TestObject$outputMat), 100)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                   "lambdaT", "data size",
                                                   "seed", "alphaE", "deltaE",
                                                   "lambdaE", "failure","time"))

  })
})

test_that("TemperedEstim_Simulation_with_Subordinator_GMM_gives_correct_return",
  {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1)), SampleSizes = 40, MCparam = 40,
      TemperedType = "Subordinator", Estimfct = "GMM", saveOutput = FALSE,
      algo = "2SGMM", regularization = "cut-off", WeightingMatrix = "OptAsym",
      t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


    expect_equal(length(TestObject$outputMat), 400)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                   "lambdaT", "data size",
                                                   "seed", "alphaE", "deltaE",
                                                   "lambdaE", "failure","time"))

  })
})


test_that("TemperedEstim_Simulation_with_Subordinator_Cgmm_gives_correct_re", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.45,0.55,1)), SampleSizes = 4,
      MCparam = 4, TemperedType = "Subordinator", Estimfct = "Cgmm",
      saveOutput = FALSE, algo = "2SCgmm", alphaReg = 0.01, subdivisions = 20,
      IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
      s_max= 1)


    expect_equal(length(TestObject$outputMat), 40)
    expect_equal(TestObject$outputMat[1],0.45)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                   "lambdaT", "data size",
                                                   "seed", "alphaE", "deltaE",
                                                   "lambdaE", "failure","time"))

  })
})

test_that("TemperedEstim_Simulation_with_Subordinator_GMC_gives_correct_re",
          {

            suppressWarnings({
              TestObject <- TemperedEstim_Simulation(
                ParameterMatrix = rbind(c(0.45,0.55,1)), SampleSizes = 4,
                MCparam = 4, TemperedType = "Subordinator", Estimfct = "GMC",
                saveOutput = FALSE, algo = "2SGMC", alphaReg = 0.01,
                WeightingMatrix = "OptAsym", regularization = "cut-off",
                ncond = 8)


              expect_equal(length(TestObject$outputMat), 40)
              expect_equal(TestObject$outputMat[1],0.45)
              expect_equal(colnames(TestObject$outputMat),
                           c("alphaT", "deltaT", "lambdaT", "data size", "seed",
                             "alphaE", "deltaE", "lambdaE", "failure","time"))

            })
          })


### Normal ###

test_that("TemperedEstim_Simulation_with_Normal_ML_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1,1,1)), SampleSizes = 10, MCparam = 10,
      TemperedType = "Normal", Estimfct = "ML", saveOutput = FALSE)


    expect_equal(length(TestObject$outputMat), 140)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT", "deltaT",
                                                   "lambdaT", "muT",
                                                   "data size", "seed",
                                                   "alphaE", "betaE", "deltaE",
                                                   "lambdaE", "muE", "failure",
                                                   "time"))

  })
})

test_that("TemperedEstim_Simulation_with_Normal_GMM_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1,1,1)), SampleSizes = 40, MCparam = 4,
      TemperedType = "Normal", Estimfct = "GMM", saveOutput = FALSE,
      algo = "2SGMM", regularization = "cut-off", WeightingMatrix = "OptAsym",
      t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


    expect_equal(length(TestObject$outputMat), 56)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT", "deltaT",
                                                   "lambdaT", "muT",
                                                   "data size", "seed",
                                                   "alphaE", "betaE", "deltaE",
                                                   "lambdaE", "muE", "failure",
                                                   "time"))

  })
})

test_that("TemperedEstim_Simulation_with_Normal_Cgmm_gives_correct_return", {

  suppressWarnings({
    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.55,0.55,1,1,1)), SampleSizes = 4,
      MCparam = 4, TemperedType = "Normal", Estimfct = "Cgmm",
      saveOutput = FALSE, algo = "2SCgmm", alphaReg = 0.01, subdivisions = 20,
      IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
      s_max= 1)


    expect_equal(length(TestObject$outputMat), 56)
    expect_equal(TestObject$outputMat[1],0.55)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT", "deltaT",
                                                   "lambdaT", "muT",
                                                   "data size", "seed",
                                                   "alphaE", "betaE", "deltaE",
                                                   "lambdaE", "muE", "failure",
                                                   "time"))

  })
})









