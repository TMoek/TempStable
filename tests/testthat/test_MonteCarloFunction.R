#library(tillsfunctionpackage)

test_that(
  "ComputeMCSimForTempered_parallel_with_Subordinator_gives_correct_return",
  {
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm",
                                                  IntegrationMethod = "Simpson",
                                                  randomIntegrationLaw = "unif")
                 [1]
                 ,1)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm",
                                                  IntegrationMethod = "Simpson",
                                                  randomIntegrationLaw = "unif")
                 [2]
                 ,0.5)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm",
                                                  IntegrationMethod = "Simpson",
                                                  randomIntegrationLaw = "unif")
                 [3]
                 ,1)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm",
                                                  IntegrationMethod = "Simpson",
                                                  randomIntegrationLaw = "unif")
                 [4]
                 ,1)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm",
                                                  IntegrationMethod = "Simpson",
                                                  randomIntegrationLaw = "unif")
                 [5]
                 ,10)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm")[6],
                 NaN)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm")[7],
                 NaN)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm")[8],
                 NaN)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm")[9],
                 1)
    expect_equal(ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,
                                                  "Subordinator","Cgmm")[10],
                 0)


    expect_warning(
      ComputeMCSimForTempered_parallel(20,c(0.5, 1, 1,1),10,
                                       "Subordinator","Cgmm"),
      "number of items to replace is not a multiple of replacement length")

  }
)




