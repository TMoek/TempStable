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




