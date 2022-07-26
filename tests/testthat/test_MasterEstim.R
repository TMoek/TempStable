#library(tillsfunctionpackage)

test_that("TemperedEstim_with_Classic_ML_gives_correct_return", {
  testData <- c(1.8873152, -0.4843727,  0.4755897, -0.1257814,  1.3484823,
                -0.3866821, -0.4258380, -0.4658479, -2.9774961,  0.9646364,
                -0.5875601, -2.0316790,  0.3641900,  1.1882307,  1.6635770,
                -0.0554876,  0.4005471,  0.7820444, -0.3786902,  1.5131663)

  suppressWarnings({
    TestObject <- TemperedEstim("Classic","ML",testData)


    expect_equal(TestObject@par[["alpha"]],1e-06)
    expect_equal(TestObject@par[["delta +"]],1e-06)
    expect_equal(TestObject@par[["delta -"]],2.9497572)
    expect_equal(TestObject@par[["lambda +"]],2.3245534)
    expect_equal(TestObject@par[["lambda -"]],1.38320472)
    expect_equal(TestObject@par[["mu"]],0.133412943)

    expect_equal(TestObject@par0,c(1.5,1,1,1,1,0))

    expect_equal(TestObject@others$par, c(0.0000010, 0.0000010, 2.9497572,
                                          2.3245534, 1.3832047, 0.1334129))

    expect_equal(TestObject@others$value, 30.25111)

    expect_equal(TestObject@others$counts[["function"]], 34)
    expect_equal(TestObject@others$counts[["gradient"]], 34)

    expect_equal(TestObject@others$message,
                 "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL")

  })

  #expect_error()

})

test_that("TemperedEstim_with_Subordinator_ML_gives_correct_return", {
  testData <- c(3.0994626, 3.9490938, 1.2214726, 1.1447818, 1.2175446,
                1.0107253, 2.3391506, 3.2869243, 2.0381353, 0.8264428,
                0.9005004, 0.5014983, 1.5031865, 0.6169632, 1.4605593,
                0.9551939, 6.2465224, 2.9974583, 3.3024097, 2.6519853)

  suppressWarnings({
    TestObject <- TemperedEstim("Subordinator","ML",testData)


    expect_equal(TestObject@par[["alpha"]],0.54827989)
    expect_equal(TestObject@par[["delta"]],0.68233765)
    expect_equal(TestObject@par[["lambda"]],0.38312306)

    expect_equal(TestObject@par0,c(0.5,1,1))

    expect_equal(TestObject@others$par, c(0.54827989, 0.68233765, 0.38312306))

    expect_equal(TestObject@others$value, 29.973696)

    expect_equal(TestObject@others$counts[["function"]], 22)
    expect_equal(TestObject@others$counts[["gradient"]], 22)

    expect_equal(TestObject@others$message,
                 "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL")

    expect_equal(TestObject@method,
                 "ML_OptimAlgo=L-BFGS-B")
  })

  expect_error(TemperedEstim("Subordinator","ML"))

})

test_that("TemperedEstim_with_Normal_ML_gives_correct_return", {
  testData <- c(1.20268116, 2.47354907, 1.50248870, 2.53324346, 4.57237410,
                1.31064953, 2.74788769, 3.27900386, 0.15987725, 2.94092597,
                0.06321566, 0.86256026, 0.97488549, 0.20086994, 1.21891900,
                5.60468051, 2.12527944, 0.24420841, 2.08260240, 1.74097067)

  suppressWarnings({
    TestObject <- TemperedEstim("Normal","ML",testData)


    expect_equal(TestObject@par[["alpha"]],1.000000e-06)
    expect_equal(TestObject@par[["beta"]],1006.67285)
    expect_equal(TestObject@par[["delta"]],1.12110825)
    expect_equal(TestObject@par[["lambda"]],632.23181)
    expect_equal(TestObject@par[["mu"]],4.1105327e-02)

    expect_equal(TestObject@par0,c(0.5,0,1,1,0))

    expect_equal(TestObject@others$par, c(1.000000e-06, 1.00667285e+03,
                                          1.12111e+00, 6.3223181e+02,
                                          4.1105e-02))

    expect_equal(TestObject@others$value, 32.302543)

    expect_equal(TestObject@others$counts[["function"]], 67)
    expect_equal(TestObject@others$counts[["gradient"]], 67)

    expect_equal(TestObject@others$message,
                 "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")

    expect_equal(TestObject@method,
                 "ML_OptimAlgo=L-BFGS-B")
  })

  expect_error(TemperedEstim("Normal","ML"))

})

test_that("TemperedEstim_with_CGMY_ML_gives_correct_return", {
  testData <- c(1.20268116, 2.47354907, 1.50248870, 2.53324346, 4.57237410,
                1.31064953, 2.74788769, 3.27900386, 0.15987725, 2.94092597,
                0.06321566, 0.86256026, 0.97488549, 0.20086994, 1.21891900,
                5.60468051, 2.12527944, 0.24420841, 2.08260240, 1.74097067)

  suppressWarnings({
    TestObject <- TemperedEstim("CGMY","ML",testData)


    expect_equal(TestObject@par[["C"]],0.330860591)
    expect_equal(TestObject@par[["G"]],0.140881690)
    expect_equal(TestObject@par[["M"]],0.0000010)
    expect_equal(TestObject@par[["Y"]],1.126007266)

    expect_equal(TestObject@par0,c(1,1,1,1.5))

    expect_equal(TestObject@others$par, c(0.330860591, 0.140881690, 0.0000010,
                                          1.126007266))

    expect_equal(TestObject@others$value, 37.15274)

    expect_equal(TestObject@others$counts[["function"]], 101)
    expect_equal(TestObject@others$counts[["gradient"]], 101)

    expect_equal(TestObject@others$message,
                 "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")

    expect_equal(TestObject@method,
                 "ML_OptimAlgo=L-BFGS-B")
  })

  expect_error(TemperedEstim("CGMY","ML"))

})


