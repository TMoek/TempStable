
#' Monte Carlo Simulation
#'
#' @description
#' Runs Monte Carlo simulation for a selected estimation method. The function
#' can save results in a file.
#'
#' @details
#' \strong{Error Handling} It is advisable to set it to TRUE when user is
#' planning to launch long simulations as it will prevent the procedure to stop
#' if an error occurs for one sample data. The estimation function will produce
#' a vector of NA as estimated parameters related to this (error generating)
#' sample data and move on to the next Monte Carlo step.
#'
#' \strong{Output file} Setting \code{saveOutput} to \code{TRUE} will have the
#' side effect of saving a csv file in the working directory. This file will
#' have \code{MCparam*length(SampleSizes)} lines and its columns will be:
#' \describe{
#'   \item{alphaT, ...:}{the true value of the parameters.}
#'   \item{data size:}{the sample size used to generate the simulated data.}
#'   \item{seed:}{the seed value used to generate the simulated data.}
#'   \item{alphaE, ...:}{the estimate of the parameters.}
#'   \item{failure:}{binary: 0 for success, 1 for failure.}
#'   \item{time:}{estimation running time in seconds.}
#' }
#' The file name is informative to let the user identify the value of the true
#' parameters, the MC parameters as well as the options selected for the
#' estimation method. The csv file is updated after each MC estimation which is
#' useful when the simulation stops before it finishes.
#'
#' \strong{SeedOptions} If user does not want to control the seed generation,
#' he could ignore this argument (default value NULL). This argument can be
#' more useful when one wants to cut the simulation (even for one parameter
#' value) into pieces. In that case, he can control which part of the seed
#' vector he wants to use.
#' \describe{
#'   \item{MCtot:}{total values of MC simulations in the entire process.}
#'   \item{seedStart:}{starting index in the seed vector. The vector extracted
#'   will be of size MCparam.}
#' }
#'
#' \strong{Estimfct} Detailed documentation of the individual parameters
#' of the various estimate functions is available in the package
#' \code{StableEstim}.
#' \describe{
#'   \item{For ML:}{use \code{?StableEstim::MLParametersEstim}. See usage of
#'   Maximum likelihood estimation in Kim et al. (2008)}
#'   \item{For GMM:}{use \code{?StableEstim::GMMParametersEstim}. Generalized
#'   Method of Moments by Hansen (1982)}
#'   \item{For Cgmm:}{use \code{?StableEstim::CgmmParametersEstim}. Continuum
#'   Generalized Methods of Moments by Carrasco & Kotchoni (2017)}
#'   \item{For GMC:}{We also use a method of moment approach which follows
#'    Kuechler & Tappe (2013). They match empirical cumulants with their
#'    theoretical counterparts. We extend this by using Hansen's (1982) GMM
#'    framework. We call the approach generalized method of cumulants (GMC) to
#'    distinguish it from the GMM method using characteristic function moment
#'    conditions. However, it fits well into Hansen's (1982) framework allowing
#'    for standard asymptotic theory.}
#' }
#'
#' @seealso
#' \url{https://github.com/GeoBosh/StableEstim/blob/master/R/Simulation.R}
#'
#' @references
#' Massing, T. (2022), 'Parametric Estimation of Tempered Stable Laws';
#'
#' Kim, Y. s., Rachev, S. T., Bianchi, M. L. & Fabozzi, F. J. (2008), 'Financial
#' market models with lévy processes and time-varying volatility'
#' \url{https://doi.org/10.1016/j.jbankfin.2007.11.004};
#'
#' Hansen, L. P. (1982), 'Large sample properties of generalized method of
#' moments estimators' \url{https://doi.org/10.2307/1912775};
#'
#' Carrasco, M. & Kotchoni, R. (2017), 'Efficient estimation using the
#' characteristic function' \url{https://doi.org/10.1017/S0266466616000025};
#'
#' Kuechler, U. & Tappe, S. (2013), 'Tempered stable distribution and processes'
#' \url{https://doi.org/10.1016/j.spa.2013.06.012};
#'
#' @param ParameterMatrix The matrix is to be composed of vectors, row by row.
#' Each vector must fit the pattern of theta of the \code{TemperedType}.
#' @param SampleSizes Sample sizes to be used to simulate the data. By default,
#'  we use 200 (small sample size) and 1600 (large sample size);
#'  vector of integer.
#' @param MCparam Number of Monte Carlo simulation for each couple of parameter,
#'  default=100; integer
#' @param TemperedType A String. Either "Classic", "Subordinator", or "Normal".
#' @param Estimfct The estimation function to be used. A String.
#'  Either "ML", "GMM", "Cgmm", or "GMC".
#' @param HandleError Logical flag: if set to TRUE, the simulation doesn't stop
#'  when an error in the estimation function is encountered. A vector of
#'  (size 4) NA is saved and the the simulation carries on. See details.
#' @param saveOutput Logical flag: if set to TRUE, a csv file (for each couple
#'  of parameter) with the the estimation
#'  information is saved in the current directory. See details.
#' @param SeedOptions List to control the seed generation. See details.
#' @param eps A gap holder.
#' @param ... Other arguments to be passed to the estimation function.
#'
#' @return If \code{saveOutput == FALSE}, the return object is a list of 2.
#' Results of the simulation are listed in \code{$outputMat}. If \code{
#' saveOutput == TRUE}, only a csv file is saved and nothing is returned.
#'
#' @examples
#' \donttest{
#' TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.5,1,1,1,1,0),
#'                                                  c(0.5,1,1,1,1,0)),
#'                          SampleSizes = 10, MCparam = 10,
#'                          TemperedType = "Classic", Estimfct = "ML",
#'                          saveOutput = FALSE)
#'
#' TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
#'                          SampleSizes = 40, MCparam = 40,
#'                          TemperedType = "Classic", Estimfct = "GMM",
#'                          saveOutput = FALSE, algo = "2SGMM",
#'                          regularization = "cut-off",
#'                          WeightingMatrix = "OptAsym", t_scheme = "free",
#'                          alphaReg = 0.005,
#'                          t_free = seq(0.1,2,length.out=12))
#'
#' TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)),
#'                          SampleSizes = 4, MCparam = 4,
#'                          TemperedType = "Classic", Estimfct = "Cgmm",
#'                          saveOutput = FALSE, algo = "2SCgmm",
#'                          alphaReg = 0.01, subdivisions = 20,
#'                          IntegrationMethod = "Uniform",
#'                          randomIntegrationLaw = "unif",
#'                          s_min = 0, s_max= 1)
#'
#' TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)),
#'                          SampleSizes = 4, MCparam = 4,
#'                          TemperedType = "Classic", Estimfct = "GMC",
#'                          saveOutput = FALSE, algo = "2SGMC",
#'                          alphaReg = 0.01, WeightingMatrix = "OptAsym",
#'                          regularization = "cut-off", ncond = 8)
#' }
#'
#' @export
TemperedEstim_Simulation <- function(ParameterMatrix,
                                     SampleSizes = c(200, 1600),
                                     MCparam = 100,
                                     TemperedType = c("Classic", "Subordinator",
                                                      "Normal"),
                                     Estimfct = c("ML", "GMM", "Cgmm", "GMC"),
                                     HandleError = TRUE, saveOutput = TRUE,
                                     SeedOptions = NULL, eps = 1e-06,
                                     parallelization = FALSE, ...) {
    #seeAlso: https://github.com/GeoBosh/StableEstim/blob/master/R/Simulation.R
    SeedVector <- getSeedVector(MCparam, SeedOptions)
    Estimfct <- match.arg(Estimfct)
    TemperedType <- match.arg(TemperedType)
    nab <- nrow(ParameterMatrix)
    npar <- ncol(ParameterMatrix)
    lS <- length(SampleSizes)
    nRowOutput <- nab * lS
    OutputCollection <- empty_list <- vector(mode = "list", length = nab)
    returnList <- empty_list <- vector(mode = "list")
    indexStatOutput <- 1

    CheckPointValues <- readCheckPoint(ParameterMatrix, TemperedType, Estimfct,
                                       nab, npar, lS, MCparam, ...)
    updatedCheckPointValues <- updateCheckPointValues(CheckPointValues, MCparam,
                                                      lS, nab)

    # Ist im Zielordner bereits eine csv-Datei mit dem gleichen Namen (thetaT
    # und MCparam sind gleich), wird die Datei aktuell um die weiteren
    # Ergebnisse aktualisiert. Mention in Details
    # Wird dieser Code aktiviert, würde das Überschreiben einer Datei verboten
    # werden.
    #
    #if (updatedCheckPointValues$mc_start != 1){
    #  print("'Can't Compute Stat summary when the process doesn't
    #        start from the beginning!!")
    #}

    for (ab in updatedCheckPointValues$ab_start:nab) {
        thetaT <- ParameterMatrix[ab, ]

        outputString <- switch(TemperedType,
                               Classic = paste("Alpha=", thetaT[1] ,
                                               " *** DeltaP=", thetaT[2],
                                               " *** DeltaM=", thetaT[3],
                                               " *** LambdaP=", thetaT[4],
                                               " *** LambdaM=", thetaT[5],
                                               " *** mu=", thetaT[6], sep = ""),
                               Subordinator = paste("Alpha=", thetaT[1] ,
                                                    " *** Delta=", thetaT[2],
                                                    " *** Lambda=", thetaT[3],
                                                    sep = ""),
                               Normal = paste("Alpha=", thetaT[1] ,
                                              " *** Beta=", thetaT[2],
                                              " *** Delta=", thetaT[3],
                                              " *** Lambda=", thetaT[4],
                                              " *** mu=", thetaT[5], sep = "")
                               # ,CGMY = paste("C=", thetaT[1] ,
                               #              " *** G=", thetaT[2],
                               #              " *** M=", thetaT[3],
                               #              " *** Y=", thetaT[4], sep = "")
                               )

        cat("---------------- ", outputString, " --------------- \n", sep = "")

        if (saveOutput) initOutputFile(thetaT, MCparam, TemperedType,
                                      Estimfct, ...)

        EstimOutput <- ComputeMCSimForTempered(thetaT = thetaT,
                                               MCparam = MCparam,
                                               SampleSizes =
                                                 as.vector(SampleSizes),
                                               SeedVector = SeedVector,
                                               TemperedType = TemperedType,
                                               Estimfct = Estimfct,
                                               HandleError = HandleError,
                                               ab_current = ab,
                                               nab = nab,
                                               npar = npar, ParameterMatrix,
                                               CheckPointValues =
                                                 updatedCheckPointValues,
                                               saveOutput = saveOutput, eps,
                                               parallelization,
                                               ...)

        OutputCollection <- EstimOutput

        if (saveOutput == FALSE){
          if (length(returnList) == 0) returnList <- EstimOutput
          else returnList <- Map(rbind,returnList, EstimOutput)
        }
    }

    deleteCheckPoint(ParameterMatrix, TemperedType, Estimfct, nab, npar, lS,
                     MCparam, ...)

    if (saveOutput == FALSE){
      return(returnList)
    }

}

# No export.
getSeedVector <- function(Outputsize, SeedOptions = NULL) {
    set.seed(345)
    if (is.null(SeedOptions))
        vec <- as.vector(sample.int(n = 3 * Outputsize, size = Outputsize))
 else {
        MCtot <- SeedOptions$MCtot
        seedStart <- SeedOptions$seedStart
        seedEnd <- seedStart + Outputsize
        vec <-
          as.vector(sample.int(n = 3 * MCtot, size = MCtot))[seedStart:seedEnd]
    }
    vec
}

# No export.
#' @importFrom foreach %dopar%
ComputeMCSimForTempered <- function(thetaT, MCparam, SampleSizes, SeedVector,
                                    TemperedType, Estimfct, HandleError,
                                    ab_current, nab, npar, ParameterMatrix,
                                    CheckPointValues = NULL, saveOutput, eps,
                                    parallelization, ...) {

    if (TemperedType == "Classic") {
        Ncol <- 16
    } else if (TemperedType == "Subordinator") {
        Ncol <- 10
    } else if (TemperedType == "Normal") {
        Ncol <- 14
    }
  # else {
  #       Ncol <- 12
  #   }
    nSS <- length(SampleSizes)
    Nrow <- nSS * MCparam
    Output <- matrix(data = NA, ncol = Ncol, nrow = Nrow)

    if (TemperedType == "Classic") {
        colnames(Output) <- c("alphaT", "delta+T", "delta-T", "lambda+T",
                              "lambda-T", "muT", "data size", "seed", "alphaE",
                              "delta+E", "delta-E", "lambda+E", "lambda-E",
                              "muE", "failure", "time")
    } else if (TemperedType == "Subordinator") {
        colnames(Output) <- c("alphaT", "deltaT", "lambdaT", "data size",
                              "seed", "alphaE", "deltaE", "lambdaE",
                              "failure", "time")
    } else if (TemperedType == "Normal") {
        colnames(Output) <- c("alphaT", "betaT", "deltaT", "lambdaT", "muT",
                              "data size", "seed", "alphaE", "betaE", "deltaE",
                              "lambdaE", "muE", "failure", "time")
    }
    # else {
    #     colnames(Output) <- c("C.T", "G.T", "M.T", "Y.T", "data size", "seed",
    #                           "C.E", "G.E", "M.E", "Y.E", "failure", "time")
    # }
    if (ab_current == CheckPointValues$ab_start) {
        sample_start = CheckPointValues$sample_start
        mc_start = CheckPointValues$mc_start
    } else {
        sample_start = 1
        mc_start = 1
    }


    #Test 17.10.22
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(5)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl,list('setClassesForeach', 'rCTS',
                                    'rCTS_aAR', 'rCTS_aARp',
                                    'getTempEstimation',
                                    'getTempEstimFcts',
                                    'TemperedEstim', 'NameParamsObjects',
                                    'MLParametersEstim_CTS',
                                    '.asymptoticVarianceEstimML_CTS',
                                    '.methodDesML_CTS',
                                    '.initResTemp',
                                    'NameParamsObjectsTemp',
                                    'writeCheckPoint',
                                    'Estim_Des_Temp',
                                    'get_filename_checkPoint_Temp',
                                    'setClassesForeach',
                                    'updateOutputFile',
                                    'get_filename',

                                    'Output', 'SeedVector', 'TemperedType',
                                    'saveOutput',
                                    'SeedVector','thetaT'),
                            envir = environment())
    parallel::clusterExport(cl,
                            varlist = ls(),
                            envir = environment())
    doRNG::registerDoRNG(1234)

    Output <- foreach::foreach(sample = sample_start:nSS,
                     .combine = "rbind",
                     .packages = c("StableEstim")
    ) %dopar%{

    #for (sample in sample_start:nSS) {

        #Predefine unknows values in this environment
        setClassesForeach()

        size <- SampleSizes[sample]
        if (sample != sample_start) mc_start = 1

        for (mc in mc_start:MCparam) {
          tIter <- getTime_()
          iter <- mc + (sample - 1) * MCparam
          set.seed(seed <- SeedVector[mc])
          if (TemperedType == "Classic") {
            x <- rCTS(n = size, alpha = thetaT[1], deltap = thetaT[2],
                      deltam = thetaT[3], lambdap = thetaT[4],
                      lambdam = thetaT[5], mu = thetaT[6])
          } else if (TemperedType == "Subordinator") {
            x <- rSTS(n = size, alpha = thetaT[1], delta = thetaT[2],
                      lambda = thetaT[3])
          } else if (TemperedType == "Normal") {
            x <- rNTS(n = size, alpha = thetaT[1], beta = thetaT[2],
                      delta = thetaT[3], lambda = thetaT[4], mu = thetaT[5])
          }

          # else {
          #     x <- rCGMY(n = size, C = thetaT[1], M = thetaT[2],
          #                G = thetaT[3], Y = thetaT[4])
          # }
          Estim <- getTempEstimation(thetaT = thetaT, x = x, seed = seed,
                                     size = size, Ncol = Ncol,
                                     TemperedType = TemperedType,
                                     Estimfct = Estimfct,
                                     HandleError = HandleError, eps, ...)

          Output[iter, ] <- Estim$outputMat
          file <- Estim$file

          if (!is.null(CheckPointValues)) {
            writeCheckPoint(ParameterMatrix, TemperedType, Estimfct, ab_current,
                            nab, npar, sample, nSS, mc,
                            MCparam, ...)
          }

          if (saveOutput) updateOutputFile(thetaT, MCparam, TemperedType,
                                           Estim)

          StableEstim::PrintEstimatedRemainingTime(iter, tIter, Nrow)
        }

        Output

        #parallelization == FALSE
        # if (isFALSE(parallelization)){
        #   #Here Code from above
        # }

        # #parallelization == TRUE
        # if (isTRUE(parallelization)){
        #   #predefine rfunction
        #   if (TemperedType == "Classic") {
        #     rfunction <- function(size, thetaT){
        #       rCTS(n = size, alpha = thetaT[1], deltap = thetaT[2],
        #            deltam = thetaT[3], lambdap = thetaT[4],
        #            lambdam = thetaT[5], mu = thetaT[6])
        #     }
        #   } else if (TemperedType == "Subordinator") {
        #     rfunction <- function(size, thetaT){
        #       rSTS(n = size, alpha = thetaT[1], delta = thetaT[2],
        #            lambda = thetaT[3])
        #     }
        #   } else { #"Normal"
        #     rfunction <- function(size, thetaT){
        #       rNTS(n = size, alpha = thetaT[1], beta = thetaT[2],
        #            delta = thetaT[3], lambda = thetaT[4], mu = thetaT[5])
        #     }
        #   }
        #
        #   cores <- parallel::detectCores()
        #   cl <- parallel::makeCluster(5)
        #   doParallel::registerDoParallel(cl)
        #   parallel::clusterExport(cl,list('setClassesForeach', 'rCTS',
        #                                   'rCTS_aAR', 'rCTS_aARp',
        #                                   'getTempEstimation',
        #                                   'getTempEstimFcts',
        #                                   'TemperedEstim', 'NameParamsObjects',
        #                                   'MLParametersEstim_CTS',
        #                                   '.asymptoticVarianceEstimML_CTS',
        #                                   '.methodDesML_CTS',
        #                                   '.initResTemp',
        #                                   'NameParamsObjectsTemp',
        #
        #                                   'sample'
        #   ))
        #   doRNG::registerDoRNG(1234)
        #
        #   #Start Loop
        #   OutputForeach <- foreach::foreach(mc = mc_start:MCparam,
        #                                     .combine = "rbind",
        #                                     .packages = c("StableEstim")
        #   ) %dopar%{
        #     #Define classes and objects again to make them reachable
        #     size <- SampleSizes[sample]
        #     setClassesForeach()
        #
        #
        #     set.seed(seed <- SeedVector[mc])
        #     x <- rfunction(size = size, thetaT = thetaT)
        #
        #     Estim <- getTempEstimation(thetaT = thetaT, x = x, seed = seed,
        #                                size = size, Ncol = Ncol,
        #                                TemperedType = TemperedType,
        #                                Estimfct = Estimfct,
        #                                HandleError = HandleError, eps, ...)
        #
        #     file <- Estim$file
        #
        #     if (saveOutput) updateOutputFile(thetaT, MCparam, TemperedType,
        #                                      Estim)
        #     Estim$outputMat
        #   }
        #   #End Loop
        #   parallel::stopCluster(cl)
        #
        #   #Translate OutputForeach in Output
        #   for(i in 0:(length(attributes(OutputForeach)$rng)-1)){
        #     mc <- i + 1
        #     iter <- mc + (sample - 1) * MCparam
        #     Output[iter, 1:(length(thetaT) + 2)] <- c(
        #       thetaT, size, SeedVector[mc])
        #     for(a in 1:length(attributes(OutputForeach)$rng[[i+1]])){
        #       Output[iter,(a+length(thetaT)+2)] <- attributes(
        #         OutputForeach)$rng[[i+1]][a]
        #     }
        #   }
        #
        # }
        # #End Parallelization
    }
    #End Sample

    #Test Today
    parallel::stopCluster(cl)

    #return(OutputForeach)

    return(list(outputMat = Output, file = file))
}

# No export.
getTempEstimation <- function(thetaT, x, seed, size, Ncol, TemperedType,
                              Estimfct, HandleError, eps, ...) {
    output <- vector(length = Ncol)
    if (TemperedType == "Classic") {
        output[1:8] <- c(thetaT, size, seed)
    } else if (TemperedType == "Subordinator") {
        output[1:5] <- c(thetaT, size, seed)
    } else if (TemperedType == "Normal") {
        output[1:7] <- c(thetaT, size, seed)
    }
    # else {
    #     output[1:6] <- c(thetaT, size, seed)
    # }
    theta0 <- thetaT - 0.1  #noise
    EstimRes <- TemperedEstim(TemperedType = TemperedType,
                              EstimMethod = Estimfct,
                              data = x, theta0 = theta0, ComputeCov = FALSE,
                              HandleError = HandleError, eps = eps, ...)
    if (TemperedType == "Classic") {
        output[9:14] <- EstimRes@par
    } else if (TemperedType == "Subordinator") {
        output[6:8] <- EstimRes@par
    } else if (TemperedType == "Normal") {
        output[8:12] <- EstimRes@par
    }
    # else {
    #     output[7:10] <- EstimRes@par
    # }
    if (TemperedType == "Classic") {
        output[15:16] <- c(EstimRes@failure, EstimRes@duration)
    } else if (TemperedType == "Subordinator") {
        output[9:10] <- c(EstimRes@failure, EstimRes@duration)
    } else if (TemperedType == "Normal") {
        output[13:14] <- c(EstimRes@failure, EstimRes@duration)
    }
    # else {
    #     output[11:12] <- c(EstimRes@failure, EstimRes@duration)
    # }
    list(outputMat = output, file = EstimRes@method)
}

# Function title
# Merge with EstimSimulation
#
# @examples
# ComputeMCSimForTempered_parallel(1,c(1.5, 1, 1, 1, 1, 0),10,"Classic","ML")
# ComputeMCSimForTempered_parallel(1,c(0.5, 1, 1),10,"Subordinator","Cgmm",
#                                  IntegrationMethod = "Simpson",
#                                  randomIntegrationLaw = "unif")
#
# No export.
# ComputeMCSimForTempered_parallel <- function(MCparam, thetaT, size,
#                                              TemperedType = c("Classic",
#                                                               "Subordinator",
#                                                               "Normal","CGMY"),
#                                              Estimfct = c("ML", "GMM", "Cgmm",
#                                                           "GMC"),
#                                              HandleError = TRUE, eps, ...) {
#     Estimfct <- match.arg(Estimfct)
#     TemperedType <- match.arg(TemperedType)
#     Ncol <- ifelse(TemperedType == "Classic", 15, 9)
#     if (TemperedType == "Classic") {
#         Ncol <- 15
#     } else if (TemperedType == "Subordinator") {
#         Ncol <- 9
#     } else if (TemperedType == "Normal") {
#         Ncol <- 13
#     } else {
#         Ncol <- 11
#     }
#     Output <- numeric(Ncol)
#     if (TemperedType == "Classic") {
#         names(Output) <- c("alphaT", "delta+T", "delta-T", "lambda+T",
#                            "lambda-T", "muT", "data size", "alphaE", "delta+E",
#                            "delta-E", "lambda+E", "lambda-E", "muE", "failure",
#                            "time")
#     } else if (TemperedType == "Subordinator") {
#         names(Output) <- c("alphaT", "deltaT", "lambdaT", "data size", "alphaE",
#                            "deltaE", "lambdaE", "failure", "time")
#     } else if (TemperedType == "Normal") {
#         names(Output) <- c("alphaT", "betaT", "deltaT", "lambdaT", "muT",
#                            "data size", "alphaE", "betaE", "deltaE", "lambdaE",
#                            "muE", "failure", "time")
#     } else {
#         names(Output) <- c("C.T", "G.T", "M.T", "Y.T", "data size", "C.E",
#                            "G.E", "M.E", "Y.E", "failure", "time")
#     }
#
#     if (TemperedType == "Classic") {
#         x <- rCTS(n = size, alpha = thetaT[1], deltap = thetaT[2],
#                   deltam = thetaT[3], lambdap = thetaT[4], lambdam = thetaT[5],
#                   mu = thetaT[6])
#     } else if (TemperedType == "Subordinator") {
#         x <- rSTS(n = size, alpha = thetaT[1], delta = thetaT[2],
#                   lambda = thetaT[3])
#     } else if (TemperedType == "Normal") {
#         x <- rNTS(n = size, alpha = thetaT[1], beta = thetaT[2],
#                   delta = thetaT[3], lambda = thetaT[4], mu = thetaT[5])
#     } else {
#         x <- rCGMY(n = size, C = theta[1], G = theta[2], G = theta[3],
#                    Y = theta[4])
#     }
#     Estim <- getTempEstimation_parallel(thetaT = thetaT, x = x, size = size,
#                                         Ncol = Ncol,
#                                         TemperedType = TemperedType,
#                                         Estimfct = Estimfct,
#                                         HandleError = HandleError,
#                                         eps = eps, ...)
#     Output <- c(MCparam, Estim)
#
#     return(Output)
# }


# No export.
# getTempEstimation_parallel <- function(thetaT, x, size, Ncol, TemperedType,
#                                        Estimfct, HandleError, eps, ...) {
#     output <- vector(length = Ncol)
#     if (TemperedType == "Classic") {
#         output[1:7] <- c(thetaT, size)
#     } else if (TemperedType == "Subordinator") {
#         output[1:4] <- c(thetaT, size)
#     } else if (TemperedType == "Normal") {
#         output[1:6] <- c(thetaT, size)
#     } else {
#         output[1:5] <- c(thetaT, size)
#     }
#     theta0 <- thetaT - 0.1  #noise
#     EstimRes <- TemperedEstim_v2(TemperedType = TemperedType,
#                                  EstimMethod = Estimfct, data = x,
#                                  theta0 = theta0, ComputeCov = FALSE,
#                                  HandleError = HandleError, eps = eps, ...)
#     if (TemperedType == "Classic") {
#         output[8:13] <- EstimRes$par
#     } else if (TemperedType == "Subordinator") {
#         output[5:7] <- EstimRes$par
#     } else if (TemperedType == "Normal") {
#         output[7:11] <- EstimRes$par
#     } else {
#         output[6:9] <- EstimRes$par
#     }
#     if (TemperedType == "Classic") {
#         output[14:15] <- c(EstimRes$failure, EstimRes$duration)
#     } else if (TemperedType == "Subordinator") {
#         output[8:9] <- c(EstimRes$failure, EstimRes$duration)
#     } else if (TemperedType == "Normal") {
#         output[12:13] <- c(EstimRes$failure, EstimRes$duration)
#     } else {
#         output[10:11] <- c(EstimRes$failure, EstimRes$duration)
#     }
#     return(output)
# }


##### for statistical summary#####

# Function title
#
# Gap holder for description.
#
# Gap holder for details.
#
# @param EstimOutput A gap holder.
# @param FctsToApply A gap holder.
# @param SampleSizes A gap holder.
# @param CheckMat A gap holder.
# @param tolFailCheck A gap holder.
# @param MCparam A gap holder.
#
# @return Gap holder for return.
#
# @export
# ComputeStatOutput <- function(EstimOutput, FctsToApply, SampleSizes, CheckMat,
#                               tolFailCheck, MCparam, ...) {
#     list(alpha = ComputeStatOutputPar(EstimOutput = EstimOutput,
#                                       FctsToApply = FctsToApply, par = "alpha",
#                                       SampleSizes = SampleSizes,
#                                       CheckMat = CheckMat,
#                                       tolFailCheck = tolFailCheck,
#                                       MCparam = MCparam, ...),
#          beta = ComputeStatOutputPar(EstimOutput = EstimOutput,
#                                      FctsToApply = FctsToApply, par = "beta",
#                                      SampleSizes = SampleSizes,
#                                      CheckMat = CheckMat,
#                                      tolFailCheck = tolFailCheck,
#                                      MCparam = MCparam, ...),
#          gamma = ComputeStatOutputPar(EstimOutput = EstimOutput,
#                                       FctsToApply = FctsToApply, par = "gamma",
#                                       SampleSizes = SampleSizes,
#                                       CheckMat = CheckMat,
#                                       tolFailCheck = tolFailCheck,
#                                       MCparam = MCparam, ...),
#          delta = ComputeStatOutputPar(EstimOutput = EstimOutput,
#                                       FctsToApply = FctsToApply, par = "delta",
#                                       SampleSizes = SampleSizes,
#                                       CheckMat = CheckMat,
#                                       tolFailCheck = tolFailCheck,
#                                       MCparam = MCparam, ...))
# }


##### for Output File#####


#  No export.
initOutputFile <- function(thetaT, MCparam, TemperedType, Estimfct, ...) {

  method <- Estim_Des_Temp(TemperedType, Estimfct, ...)
  fileName <- get_filename(thetaT, MCparam, TemperedType, method)
  if (!file.exists(fileName)) {
    x <- switch(TemperedType,
                Classic = paste("alphaT", "delta+T", "delta-T", "lambda+T",
                                "lambda-T", "muT", "data size", "seed",
                                "alphaE", "delta+E", "delta-E", "lambda+E",
                                "lambda-E", "muE", "failure", "time",
                                sep = ","),
                Subordinator =  paste("alphaT", "deltaT", "lambdaT",
                                      "data size", "seed", "alphaE", "deltaE",
                                      "lambdaE", "failure", "time",
                                      sep = ","),
                Normal = paste("alphaT", "betaT", "deltaT", "lambdaT", "muT",
                               "data size", "seed", "alphaE", "betaE", "deltaE",
                               "lambdaE", "muE", "failure", "time",
                               sep = ",")
                # ,CGMY = paste("CT", "GT", "MT", "YT", "data size", "seed", "CE",
                #              "GE", "ME", "YE","failure", "time",
                #              sep = ",")
                )

    write(x, file = fileName, sep = "\n")
  }
}


# No export.
Estim_Des_Temp <- function(TemperedType = c("Classic", "Subordinator",
                                            "Normal"),
                           EstimMethod = c("ML", "GMM", "Cgmm", "GMC"), ...) {
    TemperedType <- match.arg(TemperedType)
    EstimMethod <- match.arg(EstimMethod)
    EstimFcts <- getTempEstimFcts(TemperedType, EstimMethod)
    EstimFcts$methodDes(...)
}

##### for Checkpoints#####


# No export.
readCheckPoint <- function(ParameterMatrix, TemperedType, Estimfct, nab, npar,
                           nSS, MCparam, ...) {
    method <- Estim_Des_Temp(TemperedType, Estimfct, ...)
    fileName <- get_filename_checkPoint_Temp(ParameterMatrix, nab, npar,
                                             MCparam, method)
    if (!file.exists(fileName)) {
        write(x = "## ab;nab;npar;sample;nSS;mc;MCparam", file = fileName,
              sep = "\n")
        ab <- 1
        sample <- 1
        mc <- 0
        write(x = paste("--", ab, nab, npar, sample, nSS, mc, MCparam,
                        sep = ";"), file = fileName, sep = "\n", append = TRUE)
    } else {
        tab <- as.numeric(utils::read.table(file = fileName,
                                            header = F, sep = ";"))
        ab <- tab[2]
        sample <- tab[5]
        mc <- tab[7]
        n_ab <- tab[3]
        n_par <- tab[4]
        n_SS <- tab[6]
        mc_Param <- tab[8]
        stopifnot(nab == n_ab, npar == n_par, nSS == n_SS, mc_Param == MCparam)
    }
    list(ab = ab, nab = nab, npar = npar, sample = sample, nSS = nSS, mc = mc,
         MCparam = MCparam)
}


# No export.
# This function writes every value of the current checkpoint in a string as
# return.
get_filename_checkPoint_Temp <- function(ParameterMatrix, nab, npar, MCparam,
                                         method) {

    # This should be adapted for every Tempered Type.
    # BUT: The filenames are getting too long and this results in errors.
    # Version 0.1.0 will not feature checkpoints during calculation.
    #
    # case: StableEstim
    # if (npar == 3) {
    #     MC <- paste(paste("alpha0=", ParameterMatrix[1, 1], sep = ""),
    #                 paste("delta0=", ParameterMatrix[1, 2], sep = ""),
    #                 paste("lambda0=", ParameterMatrix[1, 3], sep = ""),
    #                 paste("alphan=", ParameterMatrix[nab, 1], sep = ""),
    #                 paste("deltan=", ParameterMatrix[nab, 2], sep = ""),
    #                 paste("lambdan=", ParameterMatrix[nab, 3], sep = ""),
    #                 paste("MCparam", MCparam, sep = ""), sep = "_")
    # } else {
    # case: "Classic"
    #     MC <- paste(paste("alpha0=", ParameterMatrix[1, 1], sep = ""),
    #                 paste("delta+0=", ParameterMatrix[1, 2], sep = ""),
    #                 paste("delta-0=", ParameterMatrix[1, 3], sep = ""),
    #                 paste("lambda+0=", ParameterMatrix[1, 4], sep = ""),
    #                 paste("lambda-0=",ParameterMatrix[1, 5], sep = ""),
    #                 paste("mu0=", ParameterMatrix[1, 6], sep = ""),
    #                 paste("alphan=", ParameterMatrix[nab,1], sep = ""),
    #                 paste("delta+n=", ParameterMatrix[nab, 2], sep = ""),
    #                 paste("delta-n=", ParameterMatrix[nab,3], sep = ""),
    #                 paste("lambda+n=", ParameterMatrix[nab, 4], sep = ""),
    #                 paste("lambda-n=", ParameterMatrix[nab,5], sep = ""),
    #                 paste("mun=", ParameterMatrix[nab, 5], sep = ""),
    #                 paste("MCparam", MCparam, sep = ""), sep = "_")
    # }

    MC <- paste("Test", sep = "")

    fileName <- paste(MC, method, "_CHECKPOINT.txt", sep = "")
    fileName
}


# No export.
updateCheckPointValues <- function(CheckPointValues, MCparam, lS, nab) {
    ab_start <- CheckPointValues$ab
    sample_start <- CheckPointValues$sample
    mc_start <- CheckPointValues$mc
    if (CheckPointValues$mc == MCparam) {
        mc_start = 1
        if (CheckPointValues$sample == lS) {
            sample_start = 1
            if (CheckPointValues$ab == nab)
                stop("Simulation finished already! check your output file")
            else ab_start = CheckPointValues$ab + 1
        } else sample_start = CheckPointValues$sample + 1
    } else mc_start = mc_start + 1
    list(ab_start = ab_start, sample_start = sample_start, mc_start = mc_start)
}

# No export.
deleteCheckPoint <- function(ParameterMatrix, TemperedType, Estimfct, nab, npar,
                             nSS, MCparam, ...) {
    method <- Estim_Des_Temp(TemperedType, Estimfct, ...)
    fileName <- get_filename_checkPoint_Temp(ParameterMatrix, nab, npar,
                                             MCparam, method)
    unlink(x = fileName, force = TRUE)
}

# No export.
writeCheckPoint <- function(ParameterMatrix, TemperedType, Estimfct, ab, nab,
                            npar, sample, nSS, mc, MCparam, ...) {
    method <- Estim_Des_Temp(TemperedType, Estimfct, ...)
    fileName <- get_filename_checkPoint_Temp(ParameterMatrix, nab, npar,
                                             MCparam, method)
    line = readLines(fileName, -1)
    line[2] = paste("--", ab, nab, npar, sample, nSS, mc, MCparam, sep = ";")
    writeLines(line, fileName)
}

#Added by Cedric 20220726
# Currently not necessary. Not even adatpted
#NameStatOutput <- function(FctsToApply, StatOutput) {
#  Names <- c("alpha", "beta", "n", names(FctsToApply), "failure", "time")
#  lapply(StatOutput, function(x) {
#    colnames(x) <- Names
#    return(x)
#  })
#}

#Added by Cedric 20220729
# No export.
get_filename <- function(thetaT, MCparam, TemperedType, method,
                         extension = ".csv") {
  MC <- switch(TemperedType,
               Classic = paste(paste("Alpha=", thetaT[1]),
                               paste("DeltaP=", thetaT[2]),
                               paste("DeltaM=", thetaT[3]),
                               paste("LambdaP=", thetaT[4]),
                               paste("LambdaM=", thetaT[5]),
                               paste("mu=", thetaT[6]),
                               "MCparam", MCparam, sep = "_"),
               Subordinator = paste(paste("Alpha=", thetaT[1]),
                                    paste("Delta=", thetaT[2]),
                                    paste("Lambda=", thetaT[3]),
                                    "MCparam", MCparam,  sep = "_"),
               Normal = paste(paste("Alpha=", thetaT[1]),
                              paste("Beta=", thetaT[2]),
                              paste("Delta=", thetaT[3]),
                              paste("Lambda=", thetaT[4]),
                              paste("mu=", thetaT[5]),
                              "MCparam", MCparam, sep = "_")
               # ,CGMY = paste(paste("C=", thetaT[1]),
               #              paste("G=", thetaT[2]),
               #              paste("M=", thetaT[3]),
               #              paste("Y=", thetaT[4]),
               #              "MCparam", MCparam, sep = "_")
               )

  fileName <- paste(MC, method, extension, sep = "")
  fileName
}


#Added by Cedric 20220805
# No export.
updateOutputFile <- function(thetaT, MCparam, TemperedType, Output){
  method <- Output$file
  fileName <- get_filename(thetaT, MCparam, TemperedType, method)

  if (!file.exists(fileName)) {
    x <- switch(TemperedType,
                Classic = paste("alphaT", "delta+T", "delta-T", "lambda+T",
                                "lambda-T", "muT", "data size", "seed",
                                "alphaE", "delta+E", "delta-E", "lambda+E",
                                "lambda-E", "muE", "failure", "time",
                                sep = ","),
                Subordinator =  paste("alphaT", "deltaT", "lambdaT",
                                      "data size", "seed", "alphaE", "deltaE",
                                      "lambdaE", "failure", "time",
                                      sep = ","),
                Normal = paste("alphaT", "betaT", "deltaT", "lambdaT", "muT",
                               "data size", "seed", "alphaE", "betaE", "deltaE",
                               "lambdaE", "muE", "failure", "time",
                               sep = ",")
                # ,CGMY = paste("CT", "GT", "MT", "YT", "data size", "seed", "CE",
                #              "GE", "ME", "YE","failure", "time",
                #              sep = ",")
                )

    write(x, file = fileName, sep = "\n")
  }

  write(x = paste(as.character(Output$outputMat), collapse=","),
        file = fileName, sep="\n", append=TRUE)
}


#Added by Cedric 20221011
# No export.
getTime_ <- function() proc.time()[3]


# Added by Cedric 20221012
# No export.
setClassesForeach <- function(){
  EstimClassicClass <- setClass(
    "EstimClassicClass",
    slots = list(
      par = "numeric",
      par0 = "numeric",
      vcov = "matrix",
      confint = "matrix",
      data = "numeric",
      sampleSize = "numeric",
      others = "list",
      duration = "numeric",
      failure = "numeric",
      method = "character"
    ),
    contains = list(),
    validity = function(object) {
      par <- object@par
      if (length(par) == 6)
        ansp <- TRUE
      else
        ansp <- "Parameter of length different of 6"
      par0 <- object@par0
      if (length(par0) == 6)
        ansp0 <- TRUE
      else
        ansp0 <- "Initial Parameter of length different of 6"
      vcov <- object@vcov
      if (ncol(vcov) == 6 &&
          nrow(vcov) == 6)
        anscov <- TRUE
      else
        anscov <- "covariance matrix of length different of 6x6"
      confint <- object@confint
      if (ncol(confint) == 2 &&
          nrow(confint) == 6)
        ansconfint <- TRUE
      else
        ansconfint <-
        "confidance intervall matrix of length different of 6x2"
      if (ansp == TRUE &&
          ansp0 == TRUE && anscov == TRUE &&
          ansconfint == TRUE)
        res <- TRUE
      else if (is.character(ansp))
        res <- ansp
      else if (is.character(ansp0))
        res <- ansp0
      else if (is.character(anscov))
        res <- anscov
      else if (is.character(ansconfint))
        res <- ansconfint
      res
    }
  )

  ## Init method

  setMethod("initialize", "EstimClassicClass",
            function(.Object,
                     par,
                     par0,
                     vcov,
                     confint,
                     method,
                     level,
                     others,
                     data,
                     duration,
                     failure,
                     ...) {
              ## handle missing
              if (missing(par))
                par        <- numeric(6)
              if (missing(par0))
                par0       <- numeric(6)
              if (missing(vcov))
                vcov       <- matrix(0, nrow = 6, ncol = 6)
              if (missing(confint))
                confint    <- matrix(0, nrow = 6, ncol = 2)
              if (missing(data))
                data       <- numeric(100)
              sampleSize <- length(data)
              if (missing(method))
                method     <- "Default"
              if (missing(others))
                others     <- list()
              if (missing(level))
                level      <- 0
              if (missing(duration))
                duration   <- 0
              if (missing(failure))
                failure    <- 0

              ## set up names
              NameParamsObjects(par, "Classic")
              NameParamsObjects(par0, "Classic")
              NameParamsObjects(vcov, "Classic")
              NameParamsObjects(confint, "Classic")
              attr(confint, "level") <- level

              methods::callNextMethod(
                .Object,
                par = par,
                par0 = par0,
                vcov = vcov,
                confint = confint,
                data = data,
                sampleSize = sampleSize,
                method = method,
                others = others,
                duration = duration,
                failure = failure,
                ...
              )
            })

  setMethod("show", "EstimClassicClass",
            function(object) {
              cat("*** Tempered Estim Classic, method Show *** \n")
              cat("** Method ** \n")
              print(object@method)
              cat("** Parameters Estimation ** \n")
              print(object@par)
              cat("** Covariance Matrix Estimation ** \n")
              print(object@vcov)
              cat("** Confidence interval Estimation ** \n")
              print(paste("Confidence level=", attributes(object@confint)$level))
              print(paste("data length=", object@sampleSize))
              print(object@confint)
              cat("** Estimation time ** \n")
              PrintDuration(object@duration)
              cat("** Estimation status ** \n")
              if (object@failure == 0)
                cat("success")
              else
                cat("failure")
              cat("\n ******* End Show (Tempered Estim Classic) ******* \n")
            })

}

