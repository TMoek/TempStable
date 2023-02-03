
### Classic ###


test_that("TemperedEstim_Simulation_with_Classic_GMM_CueGMM_cutoff", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
    SampleSizes = 4, MCparam = 1, TemperedType = "Classic",
    Estimfct = "GMM", saveOutput = FALSE, algo = "CueGMM",
    regularization = "cut-off", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))

  expect_equal(length(TestObject$outputMat),16)
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

test_that("TemperedEstim_Simulation_with_Classic_GMM_CueGMM_Tikhonov", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
    SampleSizes = 4, MCparam = 1, TemperedType = "Classic",
    Estimfct = "GMM", saveOutput = FALSE, algo = "CueGMM",
    regularization = "Tikhonov", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


  expect_equal(length(TestObject$outputMat),16)
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

test_that("TemperedEstim_Simulation_with_Classic_GMM_CueGMM_LF", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
    SampleSizes = 4, MCparam = 1, TemperedType = "Classic",
    Estimfct = "GMM", saveOutput = FALSE, algo = "CueGMM",
    regularization = "LF", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


  expect_equal(length(TestObject$outputMat),16)
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

test_that("TemperedEstim_Simulation_with_Classic_GMM_ITGMM_cutoff", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
    SampleSizes = 4, MCparam = 1, TemperedType = "Classic",
    Estimfct = "GMM", saveOutput = FALSE, algo = "ITGMM",
    regularization = "cut-off", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))


  expect_equal(length(TestObject$outputMat),16)
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

test_that("TemperedEstim_Simulation_with_Classic_GMM_CueGMM_cutoff_IC", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
    SampleSizes = 4, MCparam = 1, TemperedType = "Classic",
    Estimfct = "GMM", saveOutput = FALSE, algo = "CueGMM",
    regularization = "cut-off", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12),
    IterationControl = list(NbIter = 8, PrintIterlogical = TRUE,
                            RelativeErrMax = 1e-10))


  expect_equal(length(TestObject$outputMat),16)
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

test_that("TemperedEstim_Simulation_with_Classic_Cgmm_CueCgmm", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)), SampleSizes = 4,
    MCparam = 1, TemperedType = "Classic", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "CueCgmm", alphaReg = 0.01, subdivisions = 20,
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

test_that("TemperedEstim_Simulation_with_Classic_Cgmm_ITCgmm", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)), SampleSizes = 2,
    MCparam = 2, TemperedType = "Classic", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "ITCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)


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

test_that("TemperedEstim_Simulation_with_Normal_Cgmm_ITCgmm", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.5,1,1,1,0)), SampleSizes = 4,
    MCparam = 1, TemperedType = "Normal", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "ITCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)

  expect_equal(length(TestObject$outputMat), 14)
  expect_equal(TestObject$outputMat[1],0.5)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT",
                                                 "deltaT", "lambdaT", "muT",
                                                 "data size", "seed",
                                                 "alphaE", "betaE",
                                                 "deltaE", "lambdaE",
                                                 "muE", "failure",
                                                 "time"))
})

test_that("TemperedEstim_Simulation_with_Normal_Cgmm_CueCgmm", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.5,1,1,1,0)), SampleSizes = 4,
    MCparam = 1, TemperedType = "Normal", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "CueCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)

  expect_equal(length(TestObject$outputMat), 14)
  expect_equal(TestObject$outputMat[1],0.5)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT",
                                                 "deltaT", "lambdaT", "muT",
                                                 "data size", "seed",
                                                 "alphaE", "betaE",
                                                 "deltaE", "lambdaE",
                                                 "muE", "failure",
                                                 "time"))
})


test_that("TemperedEstim_Simulation_with_Classic_GMC_IT", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)), SampleSizes = 4,
    MCparam = 1, TemperedType = "Classic", Estimfct = "GMC",
    saveOutput = FALSE, algo = "ITGMC", alphaReg = 0.01,
    WeightingMatrix = "OptAsym", regularization = "cut-off", ncond = 8)

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


### Subordinator ###

test_that("TemperedEstim_Simulation_with_Subordinator_GMM_IT",
  {

    TestObject <- TemperedEstim_Simulation(
      ParameterMatrix = rbind(c(0.5,1,1)), SampleSizes = 4, MCparam = 1,
      TemperedType = "Subordinator", Estimfct = "GMM", saveOutput = FALSE,
      algo = "ITGMM", regularization = "cut-off", WeightingMatrix = "OptAsym",
      t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))

    expect_equal(length(TestObject$outputMat), 10)
    expect_equal(TestObject$outputMat[1],0.5)
    expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                   "lambdaT", "data size",
                                                   "seed", "alphaE", "deltaE",
                                                   "lambdaE", "failure","time"))
})

test_that("TemperedEstim_Simulation_with_Subordinator_GMM_Cue",
    {

      TestObject <- TemperedEstim_Simulation(
        ParameterMatrix = rbind(c(0.5,1,1)), SampleSizes = 4, MCparam = 1,
        TemperedType = "Subordinator", Estimfct = "GMM", saveOutput = FALSE,
        algo = "CueGMM", regularization = "cut-off", WeightingMatrix = "OptAsym",
        t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))

      expect_equal(length(TestObject$outputMat), 10)
      expect_equal(TestObject$outputMat[1],0.5)
      expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                     "lambdaT", "data size",
                                                     "seed", "alphaE", "deltaE",
                                                     "lambdaE", "failure",
                                                     "time"))
})


test_that("TemperedEstim_Simulation_with_Subordinator_Cgmm_IT", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.45,0.55,1)), SampleSizes = 2,
    MCparam = 2, TemperedType = "Subordinator", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "ITCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)

  expect_equal(length(TestObject$outputMat), 20)
  expect_equal(TestObject$outputMat[1],0.45)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                 "lambdaT", "data size",
                                                 "seed", "alphaE", "deltaE",
                                                 "lambdaE", "failure","time"))
})

test_that("TemperedEstim_Simulation_with_Subordinator_Cgmm_Cue", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.45,0.55,1)), SampleSizes = 4,
    MCparam = 1, TemperedType = "Subordinator", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "CueCgmm", alphaReg = 0.01, subdivisions = 20,
    IntegrationMethod = "Uniform", randomIntegrationLaw = "unif", s_min = 0,
    s_max= 1)


  expect_equal(length(TestObject$outputMat), 10)
  expect_equal(TestObject$outputMat[1],0.45)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "deltaT",
                                                 "lambdaT", "data size",
                                                 "seed", "alphaE", "deltaE",
                                                 "lambdaE", "failure","time"))
})


### Normal ###

test_that("TemperedEstim_Simulation_with_Normal_GMM_IT", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.5,1,1,1,1)), SampleSizes = 4, MCparam = 1,
    TemperedType = "Normal", Estimfct = "GMM", saveOutput = FALSE,
    algo = "ITGMM", regularization = "cut-off", WeightingMatrix = "OptAsym",
    t_scheme = "free", alphaReg = 0.005, t_free = seq(0.1,2,length.out=12))

  expect_equal(length(TestObject$outputMat), 14)
  expect_equal(TestObject$outputMat[1],0.5)
  expect_equal(colnames(TestObject$outputMat), c("alphaT", "betaT", "deltaT",
                                                 "lambdaT", "muT",
                                                 "data size", "seed",
                                                 "alphaE", "betaE", "deltaE",
                                                 "lambdaE", "muE", "failure",
                                                 "time"))
})

test_that("TemperedEstim_Simulation_with_Normal_Cgmm_Cue", {

  TestObject <- TemperedEstim_Simulation(
    ParameterMatrix = rbind(c(0.55,0.55,1,1,1)), SampleSizes = 4,
    MCparam = 1, TemperedType = "Normal", Estimfct = "Cgmm",
    saveOutput = FALSE, algo = "CueCgmm", alphaReg = 0.01, subdivisions = 20,
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

## GMC Tests ##

test_that("TemperedEstim_Simulation_with_Subordinator_GMC_IT", {
  TestObject <- TemperedEstim_Simulation(ParameterMatrix = rbind(c(0.5,1,1)),
                                         SampleSizes = c(2), MCparam = 2,
                                         TemperedType = "Subordinator",
                                         Estimfct = "GMC",
                                         saveOutput = FALSE, algo = "ITGMC",
                                         alphaReg = 0.01,
                                         WeightingMatrix = "OptAsym",
                                         regularization = "cut-off",
                                         ncond = 8)
  expect_equal(all(!is.na(TestObject)), TRUE)
})

test_that("TemperedEstim_Simulation_with_Classic_GMC_IT", {
  TestObject <- TemperedEstim_Simulation(ParameterMatrix =
                                           rbind(c(1.45,0.55,1,1,1,0)),
                                         SampleSizes = c(4), MCparam = 1,
                                         TemperedType = "Classic",
                                         Estimfct = "GMC",
                                         saveOutput = FALSE, algo = "ITGMC",
                                         alphaReg = 0.01,
                                         WeightingMatrix = "OptAsym",
                                         regularization = "cut-off",
                                         ncond = 8)
  expect_equal(all(!is.na(TestObject)), TRUE)
})

test_that("TemperedEstim_Simulation_with_Classic_GMC_2S", {
  TestObject <- TemperedEstim_Simulation(ParameterMatrix =
                                           rbind(c(1.45,0.55,1,1,1,0)),
                                         SampleSizes = c(2), MCparam = 1,
                                         TemperedType = "Classic",
                                         Estimfct = "GMC",
                                         saveOutput = FALSE, algo = "2SGMC",
                                         alphaReg = 0.01,
                                         WeightingMatrix = "OptAsym",
                                         regularization = "cut-off",
                                         ncond = 8)
  expect_equal(all(!is.na(TestObject)), TRUE)
})

test_that("TemperedEstim_Simulation_with_Classic_GMC_Cue", {
  TestObject <- TemperedEstim_Simulation(ParameterMatrix =
                                           rbind(c(1.45,0.55,1,1,1,0)),
                                         SampleSizes = c(4), MCparam = 1,
                                         TemperedType = "Classic",
                                         Estimfct = "GMC",
                                         saveOutput = FALSE, algo = "CueGMC",
                                         alphaReg = 0.01,
                                         WeightingMatrix = "OptAsym",
                                         regularization = "cut-off",
                                         ncond = 8)
  expect_equal(all(!is.na(TestObject)), TRUE)
})








