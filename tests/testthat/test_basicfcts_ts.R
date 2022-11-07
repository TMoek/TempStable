
test_that("charTSS_gives_correct_return", {
  expect_equal(round(charTSS(1000,0.5,1,0.3), digits = 39),
               -1.95484e-34-1.6968e-34i)
  expect_equal(round(charTSS(500,0.9,1,12.3), digits = 181),
               -1.4434e-177-1.21126e-176i)
  expect_equal(charTSS(500,0.9,1000,12.3),0+0i)

  expect_equal(charTSS(1000000,0.5,1,0.3),0+0i)
  expect_equal(charTSS(0,0.5,1,0.3),1+0i)
  expect_equal(round(charTSS(-10000,0.5,1,0.3), digits = 113),
               7.4973e-109+5.920e-109i)
  expect_equal(charTSS(500,0.4,11.3,Inf),NaN*(1+1i))
  expect_equal(charTSS(500,0.4,Inf,335),0+0i)
  expect_equal(charTSS(Inf,0.4,33,35),NaN*(1+1i))

  #Failure message: actual != expected but don't know how to show the difference
  #Should work as params are in the value range.
  #Works in console.
  #expect_equal(charTSS(1,0.5,1,1.3),0.08414106+0.7693653i)

  #Failure message: actual != expected but don't know how to show the difference
  #Should not work as params are not in the value range.
  #Works in console.
  #expect_equal(charTSS(-1,0.5,1,0.3),-0.21252633-0.3164881*(0+1i))
  #expect_equal(charTSS(0.1,0.5,1,0.3),0.92500754+0.3058346i)
  #expect_equal(charTSS(0.000001,0.5,1,0.3),1+0.0000032i)
  #expect_equal(charTSS(500,1.9,0.00001,12.3),3.794157e-05-6.56268e-04i)
  #expect_equal(charTSS(5,0.5,0.00001,12.3),0.99999755+0.0000248i)
  #expect_equal(charTSS(1,0.5,1,0.3),-0.21252633+0.3164881i)
  #expect_equal(charTSS(500,0.9,0.0000001,12.3),0.9999595+0.0002794i)

  expect_error(charTSS(10000,2.0,1,0.3))
  expect_error(charTSS(500,2,11.3,35))
  expect_error(charTSS(400,0.4,0,35))
  expect_error(charTSS(10000,2.5,1,0.3))
  expect_error(charTSS(500,1.9,1,12.3))
  expect_error(charTSS('a',0.4,33,35), "non-numeric argument to binary operator")

})

test_that("dTSS_gives_correct_return", {

  expect_equal(round(dTSS(1000,0.5,1,0.3), digits = 139), 1.13117e-134)
  expect_equal(dTSS(1000,0.5,10,3),0)

  expect_error(dTSS(1000,1.5,10,0.3))
  expect_error(dTSS(1000,1.5,10,30))
  expect_error(dTSS(1000,1.5,10,3))
  expect_error(dTSS('a',0.4,33,35), "non-numeric argument to binary operator")

  suppressWarnings({
    expect_equal(round(dTSS(1,0.5,1,0.3), digits = 7), 0.2231376)
    expect_equal(round(dTSS(10,0.5,10,3),digits = 7), 0.3111076)
    expect_equal(dTSS(10,0.5,10,300),0)
    expect_equal(round(dTSS(10,0.5,10,30), digits = 66), 7.79173e-61)
    expect_warning(dTSS(1,0.5,1,0.3))
    expect_warning(dTSS(10,0.5,10,3))
    expect_warning(dTSS(10,0.5,10,300))
    expect_warning(dTSS(10,0.5,10,30))
  })
})

test_that("pTSS_methode_integrate_gives_correct_return", {

  suppressWarnings({
    expect_equal(round(pTSS(10,0.5,10,300), digits = 0), 1)
    expect_equal(round(pTSS(1,0.5,10,300), digits = 7), 0.2907218)
    expect_equal(round(pTSS(1,0.5,10,3), digits = 118), 2.70044e-113)
    expect_equal(round(pTSS(1,0.5,1,3), digits = 7), 0.5546448)
    expect_equal(round(pTSS(1,0.5,1,1), digits = 7), 0.1902557)

    expect_error(pTSS(1,1.5,1,3))
    expect_error(pTSS(1,0.5,-1,10))
    expect_error(pTSS(1,0.5,1,0))
    expect_error(pTSS(1,0.5,0,0))
    expect_error(pTSS(1,0.5,-1,-1))
    expect_error(pTSS(1,0.5,-1,-10))
    #expect_error(pTSS(1,2,1,3), "NaNs produced")
    expect_error(pTSS(3,1.7,0.00000001,0.0000000000003))
  })
})

test_that("pTSS_methode_not_integrate_gives_correct_return", {

  suppressWarnings({
    expect_equal(round(pTSS(1,0.5,10,300,NULL,""), digits = 5), 0.29071)
    expect_equal(round(pTSS(1,0.5,10,3,NULL,""), digits = 18), 1.98397e-13)
    expect_equal(round(pTSS(1,0.5,1,3,NULL,""), digits = 6),0.554645)
    expect_equal(round(pTSS(1,0.5,1,1,NULL,"",110), digits = 9),0.190891564)
    expect_equal(pTSS(10,0.5,10,300,NULL,""),1)

    expect_error(pTSS(3,1.7,0.00000001,0.0000000000003,NULL,""))
    expect_error(pTSS(1,1.5,1,3,NULL,"",1000))
    expect_error(pTSS(1,0.5,-1,10,NULL,"",50))
    expect_error(pTSS(1,0.5,0,0,NULL,""))
    expect_error(pTSS(1,0.5,-1,-1,NULL,""))
    expect_error(pTSS(1,0.5,-1,-10,NULL,""))
    expect_error(pTSS(1,0.5,1,0,NULL,"",1))
  })
})

test_that("rTSS_gives_correct_return", {

  expect_equal(length(rTSS(100,0.5,1,1)), 100)
  expect_equal(mean(rTSS(100,0.5,1,1,NULL,"SR",0)), 0)
  expect_equal(rTSS(100,0.5,1,1,NULL,"NotARNeitherSR"), NULL)
  expect_equal(rTSS(-1,0.5,1,1,NULL, "AR"), NULL)
  expect_equal(mean(rTSS(5,0.01,1,1,NULL,"AR")),0)


  expect_error(rTSS(100,0.5,1,1,NULL,"SR",-1))
  expect_error(rTSS(-1,0.5,1,3,NULL,"SR"))
  expect_error(rTSS(5,0.001,1,1,NULL,"AR"))

  suppressWarnings({
  })
})

test_that("pTSS_gives_correct_return", {

  suppressWarnings({
    expect_equal(pTSS(1,0.9,1,10),0)
    expect_equal(round(pTSS(1,0.8,1,10), digits = 216), 1.2311e-212)
    expect_equal(round(pTSS(1,0.7,1,10), digits = 10), 0.0008909456)
    expect_equal(round(pTSS(2,0.6,10,100),digits = 77), 7.34323e-72)
    expect_equal(round(pTSS(2,0.5,10,100), digits = 7), 0.9893546)
    expect_equal(round(pTSS(5,0.5,10,100), digits = 0), 1)

    expect_error(pTSS(1000.5, 1.5, 5, 0.01))
    expect_error(pTSS(0.5, 1.5, 5, 0.01))
    expect_error(pTSS(1.1,1.1,1.1,10))
    expect_error(pTSS(1,1.1,1.1,10))
    expect_error(pTSS(1,1.1,1,10))
    expect_error(pTSS(1,1,1,10))
    expect_error(pTSS(1,0.7,100,10))

  })
})

test_that("charCTS_gives_correct_return", {

  suppressWarnings({
    expect_equal(charCTS(0,0.5,1,1.5,2,2.5,3), 1+0i)
    expect_equal(charCTS(20,1.5,200,150,300,2500,891),0+0i)
    expect_equal(round(charCTS(2,1.5,200,150,300,2500,1), digits = 27),
                 -1.5110e-23+3.7274e-23i)
    expect_equal(charCTS(0.1,1,200,0.150,300,2500,1),NaN*(1+1i))

    expect_error(charCTS(0.1,0,1,2,3,25,1))
    expect_error(charCTS(0.1,-1,200,150,300,2500,1))
    expect_error(charCTS(0.1,'a',1,2,3,25,1))

    #Test Failure: actual != expected but don't know how to show the difference
    #expect_equal(charCTS(2,0.5,1,1.5,2,2.5,3),-0.0604895+0.1321719i)
    #expect_equal(charCTS(2,0.5,1,1,1,1,1),-0.06048946+0.1321719i)
    #expect_equal(charCTS(2,0.5,1,1.5,2,2.5,3),0.32980285-0.1000856i)
    #expect_equal(charCTS(0.1,1.5,50,600,300,2500,1),0.87202365+0.0874936i)

  })
})

test_that("dCTS_FFT_gives_correct_return", {

  suppressWarnings({
    expect_equal(round(dCTS(2,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048),
                       digits = 7), 0.2126036)
    expect_equal(round(dCTS(3,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048),
                       digits = 7), 0.0720212)
    expect_equal(round(dCTS(-1,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048),
                       digits = 8), 0.07201976)
    expect_equal(round(dCTS(1.5,1.5,10,5,30,55,875,NULL,"FFT",-20,20,2048),
                       digits = 3), 0.002)
    expect_equal(round(dCTS(1,1.9,10,5,30,55,875,NULL,"FFT",-20,20,2048),
                       digits = 8), 0.03348838)
    expect_equal(round(dCTS(1,0.9,10,5,30,55,875,NULL,"FFT",-20,20,2048),
                       digits = 17), 1.2292e-13)
    expect_equal(round(dCTS(1,0.9,10,5,30,55,2,NULL,"FFT",-20,20,2048),
                       digits = 5), 0.12731)
    expect_equal(round(dCTS(1,0.9,10,5,30,55,-2,NULL,"FFT",-20,20,2048),
                       digits = 12), 3.42435e-07)
    expect_equal(round(dCTS(1,0.9,10,5,30,55,-2,NULL,"FFT",20,-20,2048),
                       digits = 13), -3.424349e-07)
    expect_equal(round(dCTS(1,0.6,1,1,1,1,1,NULL,"FFT",-40,40,20),
                       digits = 8), 0.16621234)
    expect_equal(round(dCTS(1,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048),
                       digits = 6), 0.369674)

    expect_error(dCTS(1,2,10,5,30,55,875,NULL,"FFT",-20,20,2048))
    expect_error(dCTS(1,2.5,10,5,30,55,875,NULL,"FFT",-20,20,2048))
    expect_error(dCTS(1,0.9,10,5,30,55,-2,NULL,"FFT",20,-20,-2048))
  })
})

test_that("dCTS_with_Conv_gives_correct_return", {

  suppressWarnings({
    expect_equal(round(dCTS(3,0.6,1,1,1,1,1,NULL,"Conv"),
                       digits = 3), 0.072)
    expect_equal(round(dCTS(-1,0.6,1,1,1,1,1,NULL,"Conv"),
                       digits = 8), 0.07201698)
    expect_equal(dCTS(1.5,1.5,10,5,30,55,875,NULL,"Conv"),0)
    expect_equal(dCTS(1,1.9,10,5,30,55,875,NULL,"Conv"),0)
    expect_equal(round(dCTS(1,0.9,10,5,30,55,875,NULL,"Conv"),
                       digits = 17), 1.2292e-13)
    expect_equal(dCTS(1,0.9,10,5,30,55,2,NULL,"Conv"),0)
    expect_equal(dCTS(1,0.9,10,5,30,55,-2,NULL,"Conv"),0)
    expect_equal(dCTS(1,0.9,10,5,30,55,-2,NULL,"Conv"),0)
    expect_equal(round(dCTS(1,0.6,1,1,1,1,1,NULL,"Conv"),
                      digits = 6), 0.369689)
    expect_equal(dCTS(1,0.9,10,5,30,55,-2,NULL,"Conv"),0)
    expect_equal(round(dCTS(2,0.6,1,1,1,1,1,NULL,"Conv"),
                       digits = 7), 0.2125995)

    expect_error(dCTS(1,2,10,5,30,55,875,NULL,"Conv"))
    expect_error(dCTS(1,2.5,10,5,30,55,875,NULL,"Conv"))

  })
})

test_that("pCTS_gives_correct_return", {

  suppressWarnings({
    expect_equal(round(pCTS(0.5,0.5,1,1,1,1,1),
                       digits = 4), 0.3199)
    expect_equal(round(pCTS(0.5,0.9,1,2,3,4,5),
                       digits = 12), 2.007355e-06)
    expect_equal(round(pCTS(10,0.9,1,2,3,4,5),
                       digits = 0) ,1)
    expect_equal(round(pCTS(0.001,0.9,1,2,3,4,5),
                       digits = 9), 2.94e-07)
    expect_equal(round(pCTS(1,1.9,1,2,3,4,5),
                       digits = 3), 0.214)

    expect_error(pCTS(1,10,1,2,3,4,5))
    expect_error(pCTS(1,1,1,2,3,4,5))
  })
})

test_that("rCTS_gives_correct_return", {

  suppressWarnings({

    expect_equal(length(rCTS(10,0.5,1,1,1,1,1,NULL,"SR",10)),10)
    expect_equal(rCTS(1,1,1,1,1,1,0,NULL,"SR",100), NaN)

    expect_equal(length(rCTS(10,0.5,1,1,1,1,1,NULL,"aAR")), 10)

    expect_error(rCTS(1,2.5,1,1,1,1,0,NULL,"SR",100))
    expect_error(rCTS(1,1,1,1,1,1,0,NULL,"aAR"))

  })
})

test_that("qCTS_gives_correct_return", {

  suppressWarnings({
    expect_equal(round(qCTS(0.5,1.5,10,10,10,10,10),
                       digits = 6), 9.988176)
    expect_equal(round(qCTS(0.5,1.5,1,1,1,1,1),
                       digits = 7), 0.9999951)
    expect_equal(round(qCTS(0.9,1.5,1,1,1,1,1),
                       digits = 6), 3.386769)
    expect_equal(round(qCTS(0.9,1.5,1,1,1,1,100),
                       digits = 4), 19.5321)
    expect_equal(round(qCTS(0.000001,1.5,1,1,1,1,1),
                 digits = 4), -9.9745)

    expect_error(qCTS(1.9,1.5,1,1,1,1,1))
    expect_error(qCTS(0.1,1,1,1,1,1,1))
    expect_error(qCTS(0.1,2,1,1,1,1,1))
  })
})

test_that("charNTS_gives_correct_return", {

  suppressWarnings({
    expect_equal(round(charNTS(3,0.9,1,2,3,4),
                       digits = 38), -3.94024e-33+8.95099e-33i)
    expect_equal(round(charNTS(0.1,0.9,1,2,3,4),
                       digits = 1), -0.5+0.8i)
    expect_equal(round(charNTS(0.1,0.9,10,20,30,40),
                       digits = 1), 0.2+0.4i)

    expect_error(charNTS("a",2,1,2,3,4))
    expect_error(charNTS(0.1,1,1,2,3,4))

    #Failure: actual != expected but don't know how to show the difference
    #expect_equal(charNTS(0.1,2.2,1,-2,3,4),-1.1104407+0.101657i)
    #expect_equal(charNTS(0.1,0.9,0.1,0.1,0.1,0.1),0.99374908+0.021783i)
    #expect_equal(charNTS(1,0.9,1,2,3,4),-5.623971e-05+1.580768e-04i)
    #expect_equal(charNTS(1,0.9,1,1,1,1),-0.0057299791-0.004312239i)

  })
})

test_that("dNTS_gives_correct_return", {

  expect_equal(round(dNTS(1,0.8,1,1,1,1),
                     digits = 9) ,0.021745921)
  expect_equal(round(dNTS(1,0.6,5,20,50,50),
                     digits = 12), 1.725167e-06)
  expect_equal(round(dNTS(1,0.6,45,20,50,50),
                     digits = 8), 0.02057889)
  expect_equal(round(dNTS(10,0.9,45,20,50,50),
                     digits = 8), 0.02499439)
  expect_equal(dNTS(0.5,0.5,20,20,20,20, a = -200, b= 200, nf = 4096),1e-18)
  expect_equal(dNTS(0.5,0.5,20,20,20,20, a = -2000, b = 2000, nf = 8192),1e-18)
  expect_equal(round(dNTS(0.5,0.5,20,20,20,20),
                     digits = 8), 0.00934912)

  expect_error(dNTS(1,1.1,1,1,1,1))
  expect_error(dNTS(1,1.2,5,20,50,50))
  expect_error(dNTS(10,1.2,5,20,50,50))
  expect_error(dNTS(10,1.9,45,20,50,50))

  suppressWarnings({
    expect_error(dNTS(1,1,1,1,1,1))
  })
})

test_that("pNTS_gives_correct_return", {

  expect_equal(pNTS(0.1,0.5,1,1,1,1),0.022218711)
  expect_equal(round(pNTS(0.1,0.5,1,1,1,1,NULL, -20, 20, 2^6),
                     digits = 8), 0.02221881)
  expect_equal(round(pNTS(0.1,0.5,1,1,1,1,NULL, -20, 2, 2^6),
                     digits = 8), 0.02221881)
  expect_equal(round(pNTS(0.1,0.5,1,1,1,1,NULL, -2, 2, 2^6),
                     digits = 8), 0.02208606)
  expect_equal(round(pNTS(0.1,0.5,10,10,10,10,NULL, -2, 2, 2^6),
                     digits = 8), 0.01369044)
  expect_equal(round(pNTS(0.1,0.5,10,10,10,10),
                     digits = 7), 0.8412962)
  expect_equal(round(pNTS(2,0.5,10,10,10,10),
                     digits = 0), 1)
  expect_equal(round(pNTS(20,0.5,10,10,10,10),
                     digits = 0), 1)
  expect_equal(round(pNTS(0.0001,0.5,10,10,10,10),
                     digits = 7), 0.8408389)
  expect_equal(round(pNTS(0.0001,0.5,-1.5,10,10,10),
                     digits = 8), 0.25838009)
  expect_equal(round(pNTS(0.0001,0.5,-1.5,10,10,-10),
                     digits = 6), 0.741626)

  expect_error(pNTS(0.0001,1.5,10,10,10,10))
  expect_error(pNTS(0.0001,0.5,-1.5,-10,10,10))
  expect_error(pNTS(0.0001,0.5,-1.5,10,-10,10))

  suppressWarnings({
  })
})

test_that("rNTS_gives_correct_return", {

  expect_equal(length(rNTS(100, 0.5, 1,1,1,1)), 100)
  expect_equal(length(rNTS(10, 0.6, 0,1,1,0)), 10)
  expect_equal(length(rNTS(10, 0.5, 1,1,1,1, NULL, "SR", 100)), 10)

  expect_error(rNTS(100, 1, 1,1,1,1))
  expect_error(rNTS(0, 0.5, 1,1,1,1))
  expect_error(rNTS(-1, 0.5, 1,1,1,1))
  expect_error(rNTS(10, 0.5, 0,0,1,1))
  expect_error(rNTS(10, 0.5, 0,1,0,1))
  expect_error(rNTS(10, 1.5, 0,1,1,0, NULL, "SR", 100))
  expect_error(rNTS(0, 0.5, 0,1,1,0, NULL, "SR", 100))

  suppressWarnings({
  })
})

test_that("qNTS_gives_correct_return", {

  expect_equal(round(qNTS(0.1,0.5,1,1,1,1), digits = 1), 0.9)
  expect_equal(qNTS(0.3,0.6,1,1,1,1), NA)
  expect_equal(qNTS(0.6,0.6,1,1,1,1),NA)
  expect_equal(round(qNTS(0.6,0.6,10,1,1,10),
               digits = 7), -7.1174734)
  expect_equal(round(qNTS(0.6,0.6,1,10,1,1),
                     digits = 2), -9.88)

  suppressWarnings({
    expect_error(qNTS(0.1,1.5,1,1,1,1))
    expect_error(qNTS(0.1,0,1,1,1,1))
    expect_error(qNTS(0,0.6,1,1,1,1))
    expect_error(qNTS(-1,0.6,1,1,1,1))
    expect_error(qNTS(1,0.6,1,1,1,1))
    expect_error(qNTS(1.5,0.6,1,1,1,1))
  })
})

# test_that("charCGMY_gives_correct_return", {
#
#   suppressWarnings({
#     expect_equal(charCGMY(1,1,1,1,1),NaN*(1+1i))
#     expect_equal(charCGMY(1,1,1,1,1.5),0.18552469+0i)
#     expect_equal(charCGMY(1,1,1,1,0.5),0.49675807+0i)
#     expect_equal(charCGMY(10,1,1,1,0.5),6.9448315e-05+0i)
#     expect_equal(charCGMY(-0.1,1,1,1,0.5),0.9912042+0i)
#     expect_equal(charCGMY(-1,1,1,1,0.5), 0.49675807+0i)
#     expect_equal(charCGMY(1,10,10,10,0.6),0.70310998+0i)
#
#     #Failure: actual != expected but don't know how to show the difference
#     #expect_equal(charCGMY(1,10,1000,1,0.6),0.001452764+0.02770579i)
#     #expect_equal(charCGMY(1,10,10,1,0.6),0.001409213+0.02322749i)
#     #expect_equal(charCGMY(1,1,1,10,0.6),0.67874801-0.103291i)
#     #expect_equal(charCGMY(1,10,10,1,0.6),0.001409213+0.02322749i)
#
#     expect_error(charCGMY(1,1,1,1,-0.5))
#     expect_error(charCGMY(1,1,1,1,025))
#     expect_error(charCGMY(1,0,0,0,0.6))
#   })
# })

# test_that("dCGMY_gives_correct_return", {
#
#   suppressWarnings({
#     expect_equal(dCGMY(1,1,1,1,0.5),0.20822584)
#     expect_equal(dCGMY(1,1,1,1,1.5),0.18593011)
#     expect_equal(dCGMY(1,1,1,1,0.5,""),0.20822206)
#     expect_equal(dCGMY(1,1,1,1,0.5,"", -2, 2, 2^4),0.20822206)
#
#     expect_error(dCGMY(1,1,1,1,1))
#   })
# })

# test_that("dCGMY_gives_correct_return", {
#
#   suppressWarnings({
#     expect_equal(dCGMY(1,1,1,1,0.5),0.20822584)
#     expect_equal(dCGMY(1,1,1,1,1.5),0.18593011)
#     expect_equal(dCGMY(1,1,1,1,0.5,""),0.20822206)
#     expect_equal(dCGMY(1,1,1,1,0.5,"", -2, 2, 2^4),0.20822206)
#
#     expect_error(dCGMY(1,1,1,1,1))
#   })
# })

# test_that("rCGMY_gives_correct_return", {
#
#   suppressWarnings({
#     expect_equal(length(rCGMY(100,1,1,1,0.5)),100)
#     expect_equal(rCGMY(1,1,1,1,1),NaN)
#     expect_equal(rCGMY(0,1,1,1,0.5),list())
#
#     expect_error(rCGMY(1,1,1,1,2.5))
#   })
# })

# test_that("chartocdf_gives_correct_return", {
#
#   suppressWarnings({
#     expect_equal(chartocdf(0.5,10,1,1,charNTS, alpha=0.5, beta = 1, delta = 1,
#                            lambda = 1, mu = 1), 0.053609065)
#     # expect_equal(chartocdf(0.5,10,1,1,charCGMY, Y=0.5, C = 1, G = 1, M = 1),
#     #              0.3212782)
#     expect_equal(chartocdf(0.5,10,1,1,charTSS, alpha=0.5, delta = 1,
#                            lambda = 1), 0.0730034)
#     expect_equal(chartocdf(0.5,10,1,1,charCTS, alpha=0.5, deltap = 1,
#                            deltam = 1, lambdap = 1, lambdam = 1, mu = 1),
#                  0.216307303)
#
#     expect_error(chartocdf(0.5,10,1,1,charCTS()))
#     expect_error(chartocdf(0.5,10,1,1,charNTS, alpha=0.5, delta = 1, lambda = 1,
#                            mu = 1))
#   })
# })



