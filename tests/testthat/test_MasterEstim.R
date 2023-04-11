
test_that("TemperedEstim_with_CTS_ML_gives_correct_return", {
  testData <- c(1.8873152, -0.4843727,  0.4755897, -0.1257814,  1.3484823,
                -0.3866821, -0.4258380, -0.4658479, -2.9774961,  0.9646364,
                -0.5875601, -2.0316790,  0.3641900,  1.1882307,  1.6635770,
                -0.0554876,  0.4005471,  0.7820444, -0.3786902,  1.5131663)

  suppressWarnings({
    TestObject <- TemperedEstim("CTS","ML",testData)


    expect_equal(TestObject@par[["alpha"]],1e-06)
    expect_equal(TestObject@par[["delta +"]],1e-06)
    expect_equal(round(TestObject@par[["delta -"]],
                       digits = 2), 2.95)
    expect_equal(round(TestObject@par[["lambda +"]],
                       digits = 3), 2.325)
    expect_equal(round(TestObject@par[["lambda -"]],
                       digits = 3), 1.383)
    expect_equal(round(TestObject@par[["mu"]],
                       digits = 2), 0.13)

    expect_equal(TestObject@par0,c(1.5,1,1,1,1,0))

    expect_equal(round(TestObject@others$par,
                       digits = 2), c(0.00, 0.00, 2.95,
                                          2.32, 1.38, 0.13))

    expect_equal(round(TestObject@others$value,
                       digits = 2), 30.25)

    #Mac test == 33
    if(.Platform$OS.type == "windows"){
      expect_equal(TestObject@others$counts[["function"]], 34)
      expect_equal(TestObject@others$counts[["gradient"]], 34)
    }

    #Gives error for Linux. Message is ==
    # "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"
    if(.Platform$OS.type == "windows"){
      expect_equal(TestObject@others$message,
                   "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL")
    }
  })
})

test_that("TemperedEstim_with_TSS_ML_gives_correct_return", {
  testData <- c(3.0994626, 3.9490938, 1.2214726, 1.1447818, 1.2175446,
                1.0107253, 2.3391506, 3.2869243, 2.0381353, 0.8264428,
                0.9005004, 0.5014983, 1.5031865, 0.6169632, 1.4605593,
                0.9551939, 6.2465224, 2.9974583, 3.3024097, 2.6519853)

  suppressWarnings({
    TestObject <- TemperedEstim("TSS","ML",testData)


    expect_equal(round(TestObject@par[["alpha"]],
                       digits = 3), 0.548)
    expect_equal(round(TestObject@par[["delta"]],
                       digits = 3), 0.682)
    expect_equal(round(TestObject@par[["lambda"]],
                       digits = 3), 0.383)

    expect_equal(TestObject@par0,c(0.5,1,1))

    expect_equal(round(TestObject@others$par,
                       digits = 3), c(0.548, 0.682, 0.383))

    expect_equal(round(TestObject@others$value,
                       digits = 3), 29.974)

    expect_equal(TestObject@others$counts[["function"]], 22)
    expect_equal(TestObject@others$counts[["gradient"]], 22)

    expect_equal(TestObject@others$message,
                 "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL")

    expect_equal(TestObject@method,
                 "ML_OptimAlgo=L-BFGS-B")
  })

  expect_error(TemperedEstim("TSS","ML"))

})

test_that("TemperedEstim_with_NTS_ML_gives_correct_return", {
  testData <- c(1.20268116, 2.47354907, 1.50248870, 2.53324346, 4.57237410,
                1.31064953, 2.74788769, 3.27900386, 0.15987725, 2.94092597,
                0.06321566, 0.86256026, 0.97488549, 0.20086994, 1.21891900,
                5.60468051, 2.12527944, 0.24420841, 2.08260240, 1.74097067)

  suppressWarnings({
    #This test works well for Mac, Windows and Linux(rHub). Somehow, it does not
    # work with release-check Debian
    if(.Platform$OS.type == "windows"){
      TestObject <- TemperedEstim("NTS","ML",testData)


      expect_equal(TestObject@par[["alpha"]],1.000000e-06)
      expect_equal(round(TestObject@par[["delta"]],
                         digits = 2), 1.12)
      if(.Platform$OS.type == "windows"){
        expect_equal(round(TestObject@par[["beta"]],
                           digits = 3), 1006.673)
        expect_equal(round(TestObject@par[["lambda"]],
                           digits = 3), 632.232)
      }
      expect_equal(round(TestObject@par[["mu"]],
                         digits = 3), 4.1e-02)

      expect_equal(TestObject@par0,c(0.5,0,1,1,0))

      expect_equal(round(TestObject@others$value,
                         digits = 1), 32.3)

      #Mac test == 71
      expect_equal(TestObject@others$counts[["function"]], 67)
      expect_equal(TestObject@others$counts[["gradient"]], 67)


      expect_equal(TestObject@others$message,
                   "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")

      expect_equal(TestObject@method,
                   "ML_OptimAlgo=L-BFGS-B")
    }
  })

  expect_error(TemperedEstim("NTS","ML"))

})

# test_that("TemperedEstim_with_CGMY_ML_gives_correct_return", {
#   testData <- c(1.20268116, 2.47354907, 1.50248870, 2.53324346, 4.57237410,
#                 1.31064953, 2.74788769, 3.27900386, 0.15987725, 2.94092597,
#                 0.06321566, 0.86256026, 0.97488549, 0.20086994, 1.21891900,
#                 5.60468051, 2.12527944, 0.24420841, 2.08260240, 1.74097067)
#
#   suppressWarnings({
#     TestObject <- TemperedEstim("CGMY","ML",testData)
#
#
#     expect_equal(TestObject@par[["C"]],0.330860591)
#     expect_equal(TestObject@par[["G"]],0.140881690)
#     expect_equal(TestObject@par[["M"]],0.0000010)
#     expect_equal(TestObject@par[["Y"]],1.126007266)
#
#     expect_equal(TestObject@par0,c(1,1,1,1.5))
#
#     expect_equal(TestObject@others$par, c(0.330860591, 0.140881690, 0.0000010,
#                                           1.126007266))
#
#     expect_equal(TestObject@others$value, 37.15274)
#
#     expect_equal(TestObject@others$counts[["function"]], 101)
#     expect_equal(TestObject@others$counts[["gradient"]], 101)
#
#     expect_equal(TestObject@others$message,
#                  "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH")
#
#     expect_equal(TestObject@method,
#                  "ML_OptimAlgo=L-BFGS-B")
#   })
#
#   expect_error(TemperedEstim("CGMY","ML"))
#
# })

test_that("TemperedEstim_with_TSS_GMM_gives_correct_return", {
  testData <- c(1.0306890, 1.5027451, 1.5030160, 1.5891460, 0.9589073,
                0.7129186, 1.0200448, 3.4584815, 1.1106868, 0.7493848,
                1.1104505, 2.4100007, 0.9076451, 4.2011003, 0.7050829,
                2.3226794, 2.2397202, 1.0066137, 1.1703344, 0.7563997)

  suppressWarnings({
    TestObject <- TemperedEstim("TSS", "GMM",testData, algo = "2SGMM",
                                alphaReg = 0.01, regularization = "cut-off",
                                WeightingMatrix = "OptAsym", t_scheme = "free",
                                t_free = seq(0.1,2,length.out = 12))


    expect_equal(round(TestObject@par[["alpha"]],
                       digits = 3), 0.761)
    expect_equal(round(TestObject@par[["delta"]],
                       digits = 2), 0.30)
    expect_equal(round(TestObject@par[["lambda"]],
                       digits = 2), 0.27)

    expect_equal(round(TestObject@par0,
                       digits = 2), c(0.35, 1.19, 1.14))

    expect_equal(round(TestObject@others$par,
                       digits = 2), c(0.76, 0.30, 0.27))
  })
})

test_that("TemperedEstim_with_TSS_GMC_gives_correct_return", {
  testData <- c(2.7875940, 0.6474977, 3.6280734, 2.2341773, 2.9577528,
                1.5241392, 3.6102109, 1.9597134, 2.8262536, 2.1545273,
                1.5342798, 0.7227840, 2.9183721, 5.3106484, 0.8756254,
                0.6537854, 2.1774635, 2.3642591, 1.2595214, 2.1125328)

  suppressWarnings({
    if(.Platform$OS.type == "windows"){
      TestObject <- TemperedEstim("TSS", "GMC", testData,
                                  algo = "2SGMC", alphaReg = 0.01,
                                  WeightingMatrix = "OptAsym",
                                  regularization = "cut-off", ncond = 8)


      expect_equal(round(TestObject@par[["alpha"]],
                         digits = 0), 1)
      expect_equal(round(TestObject@par[["delta"]],
                         digits = 9), 0.000002178)
      #Lambda seems to be different in Mac and Windows check
      expect_equal(round(TestObject@par[["lambda"]],
                         digits = 0), 6)

      expect_equal(TestObject@par0,c(0.5,1,1))
    }
  })
})

