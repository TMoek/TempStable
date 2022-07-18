
#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param ParameterMatrix A gap holder.
#' @param SampleSizes A gap holder.
#' @param MCparam A gap holder.
#' @param TemperedType A String. Either "Classic", "Subordinator", "Normal", or
#' "CGMY".
#' @param Estimfct A String. Either "ML", "GMM", "Cgmm", or "GMC".
#' @param HandleError A Boolean. \code{TRUE} by default.
#' @param FctsToApply A gap holder.
#' @param saveOutput A gap holder.
#' @param StatSummary A gap holder.
#' @param CheckMat A gap holder.
#' @param tolFailCheck A gap holder.
#' @param SeedOptions A gap holder.
#' @param eps A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
TemperedEstim_Simulation <- function(ParameterMatrix,
                                     SampleSizes = c(200, 1600),
                                     MCparam = 100,
                                     TemperedType = c("Classic", "Subordinator",
                                                      "Normal", "CGMY"),
                                     Estimfct = c("ML", "GMM", "Cgmm", "GMC"),
                                     HandleError = TRUE, FctsToApply = StatFcts,
                                     saveOutput = TRUE, StatSummary = FALSE,
                                     CheckMat = TRUE, tolFailCheck = tolFailure,
                                     SeedOptions = NULL, eps, ...) {
    SeedVector <- getSeedVector(MCparam, SeedOptions)
    Estimfct <- match.arg(Estimfct)
    TemperedType <- match.arg(TemperedType)
    nab <- nrow(ParameterMatrix)
    npar <- ncol(ParameterMatrix)
    lS <- length(SampleSizes)
    nRowOutput <- nab * lS
    OutputCollection <- empty_list <- vector(mode = "list", length = nab)
    indexStatOutput <- 1
    CheckPointValues <- readCheckPoint(ParameterMatrix, TemperedType, Estimfct,
                                       nab, npar, lS, MCparam, ...)
    updatedCheckPointValues <- updateCheckPointValues(CheckPointValues, MCparam,
                                                      lS, nab)
    # if (updatedCheckPointValues$mc_start != 1 && StatSummary) { print('Can't Compute Stat summary when the process
    # doesn't start from the beginning!!') StatSummary = FALSE } if (StatSummary) { if(npar == 6){ StatOutputLength
    # <- length(FctsToApply) + 5 } else { } StatOutputLength <- length(FctsToApply) + 5 StatOutput <- list(alpha =
    # matrix(data = NA, ncol = StatOutputLength, nrow = nRowOutput), beta = matrix(data = NA, ncol =
    # StatOutputLength, nrow = nRowOutput), gamma = matrix(data = NA, ncol = StatOutputLength, nrow = nRowOutput),
    # delta = matrix(data = NA, ncol = StatOutputLength, nrow = nRowOutput)) }
    for (ab in updatedCheckPointValues$ab_start:nab) {
        thetaT <- ParameterMatrix[ab, ]
        cat("---------------- theta=", thetaT, sep = "")
        # if (saveOutput) initOutputFile(thetaT, MCparam, TemperedType, Estimfct, ...)
        EstimOutput <- ComputeMCSimForTempered(thetaT = thetaT,
                                               MCparam = MCparam,
                                               SampleSizes =
                                                 as.vector(SampleSizes),
                                               SeedVector = SeedVector,
                                               TemperedType = TemperedType,
                                               Estimfct = Estimfct,
                                               HandleError = HandleError,
                                               ab_current = ab, nab = nab,
                                               npar = npar, ParameterMatrix,
                                               CheckPointValues =
                                                 updatedCheckPointValues,
                                               saveOutput = saveOutput,eps, ...)
        # if (StatSummary) { res <- ComputeStatOutput(EstimOutput = EstimOutput$outputMat, FctsToApply = FctsToApply,
        # SampleSizes = SampleSizes, CheckMat = CheckMat, tolFailCheck = tolFailCheck, MCparam = MCparam, ...)
        # IndexSec <- seq(indexStatOutput, indexStatOutput + (lS - 1), 1) StatOutput$alpha[IndexSec, ] <- res$alpha
        # StatOutput$beta[IndexSec, ] <- res$beta StatOutput$gamma[IndexSec, ] <- res$gamma
        # StatOutput$delta[IndexSec, ] <- res$delta indexStatOutput <- indexStatOutput + lS }
        OutputCollection <- EstimOutput
    }
    deleteCheckPoint(ParameterMatrix, TemperedType, Estimfct, nab, npar, lS,
                     MCparam, ...)
    if (StatSummary)
        return(NameStatOutput(FctsToApply, StatOutput))
    #TODO: Woher kommt die Funktion NameStatOutput(...)? Es gibt die Funktion weder im Package, noch in den anderen

}

#' No export.
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

#' No export.
ComputeMCSimForTempered <- function(thetaT, MCparam, SampleSizes, SeedVector,
                                    TemperedType, Estimfct, HandleError,
                                    ab_current,nab, npar, ParameterMatrix,
                                    CheckPointValues = NULL, SaveOutput = TRUE,
                                    eps, ...) {
    if (TemperedType == "Classic") {
        Ncol <- 16
    } else if (TemperedType == "Subordinator") {
        Ncol <- 10
    } else if (TemperedType == "Normal") {
        Ncol <- 14
    } else {
        Ncol <- 12
    }
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
    } else {
        colnames(Output) <- c("C.T", "G.T", "M.T", "Y.T", "data size", "seed",
                              "C.E", "G.E", "M.E", "Y.E", "failure", "time")
    }
    if (ab_current == CheckPointValues$ab_start) {
        sample_start = CheckPointValues$sample_start
        mc_start = CheckPointValues$mc_start
    } else {
        sample_start = 1
        mc_start = 1
    }
    for (sample in sample_start:nSS) {
        size <- SampleSizes[sample]
        if (sample != sample_start)
            mc_start = 1
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
            } else {
                x <- rCGMY(n = size, C = theta[1], G = theta[2], G = theta[3],
                           Y = theta[4])
            }
            Estim <- getTempEstimation(thetaT = thetaT, x = x, seed = seed,
                                       size = size, Ncol = Ncol,
                                       TemperedType = TemperedType,
                                       Estimfct = Estimfct,
                                       HandleError = HandleError,
                                       eps = eps, ...)
            Output[iter, ] <- Estim$outputMat
            file <- Estim$file
            if (!is.null(CheckPointValues)) {
                writeCheckPoint(ParameterMatrix, Estimfct, ab_current, nab,
                                npar, sample, nSS, mc, MCparam, ...)
            }

            # if (SaveOutput) updateOutputFile(alphaT, betaT, MCparam, Estim)
            StableEstim::PrintEstimatedRemainingTime(iter, tIter, Nrow)
        }
    }
    list(outputMat = Output, file = file)
}

#' No export.
getTempEstimation <- function(thetaT, x, seed, size, Ncol, TemperedType,
                              Estimfct, HandleError, eps, ...) {
    output <- vector(length = Ncol)
    if (TemperedType == "Classic") {
        output[1:8] <- c(thetaT, size, seed)
    } else if (TemperedType == "Subordinator") {
        output[1:5] <- c(thetaT, size, seed)
    } else if (TemperedType == "Normal") {
        output[1:7] <- c(thetaT, size, seed)
    } else {
        output[1:6] <- c(thetaT, size, seed)
    }
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
    } else {
        output[7:10] <- EstimRes@par
    }
    if (TemperedType == "Classic") {
        output[15:16] <- c(EstimRes@failure, EstimRes@duration)
    } else if (TemperedType == "Subordinator") {
        output[9:10] <- c(EstimRes@failure, EstimRes@duration)
    } else if (TemperedType == "Normal") {
        output[13:14] <- c(EstimRes@failure, EstimRes@duration)
    } else {
        output[11:12] <- c(EstimRes@failure, EstimRes@duration)
    }
    list(outputMat = output, file = EstimRes@method)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param MCparam A gap holder.
#' @param thetaT A gap holder.
#' @param size A gap holder.
#' @param TemperedType A String. Either "Classic", "Subordinator", "Normal", or
#' "CGMY".
#' @param Estimfct A String. Either "ML", "GMM", "Cgmm", or "GMC".
#' @param HandleError A Boolean. \code{TRUE} by default.
#' @param eps A gap holder.
#'
#' @return Gap holder for return.
#'
#' @example
#' Look in Archive for examples.
#'
#' @export
ComputeMCSimForTempered_parallel <- function(MCparam, thetaT, size,
                                             TemperedType = c("Classic",
                                                              "Subordinator",
                                                              "Normal","CGMY"),
                                             Estimfct = c("ML", "GMM", "Cgmm",
                                                          "GMC"),
                                             HandleError = TRUE, eps, ...) {
    Estimfct <- match.arg(Estimfct)
    TemperedType <- match.arg(TemperedType)
    Ncol <- ifelse(TemperedType == "Classic", 15, 9)
    if (TemperedType == "Classic") {
        Ncol <- 15
    } else if (TemperedType == "Subordinator") {
        Ncol <- 9
    } else if (TemperedType == "Normal") {
        Ncol <- 13
    } else {
        Ncol <- 11
    }
    Output <- numeric(Ncol)
    if (TemperedType == "Classic") {
        names(Output) <- c("alphaT", "delta+T", "delta-T", "lambda+T",
                           "lambda-T", "muT", "data size", "alphaE", "delta+E",
                           "delta-E", "lambda+E", "lambda-E", "muE", "failure",
                           "time")
    } else if (TemperedType == "Subordinator") {
        names(Output) <- c("alphaT", "deltaT", "lambdaT", "data size", "alphaE",
                           "deltaE", "lambdaE", "failure", "time")
    } else if (TemperedType == "Normal") {
        names(Output) <- c("alphaT", "betaT", "deltaT", "lambdaT", "muT",
                           "data size", "alphaE", "betaE", "deltaE", "lambdaE",
                           "muE", "failure", "time")
    } else {
        names(Output) <- c("C.T", "G.T", "M.T", "Y.T", "data size", "C.E",
                           "G.E", "M.E", "Y.E", "failure", "time")
    }

    if (TemperedType == "Classic") {
        x <- rCTS(n = size, alpha = thetaT[1], deltap = thetaT[2],
                  deltam = thetaT[3], lambdap = thetaT[4], lambdam = thetaT[5],
                  mu = thetaT[6])
    } else if (TemperedType == "Subordinator") {
        x <- rSTS(n = size, alpha = thetaT[1], delta = thetaT[2],
                  lambda = thetaT[3])
    } else if (TemperedType == "Normal") {
        x <- rNTS(n = size, alpha = thetaT[1], beta = thetaT[2],
                  delta = thetaT[3], lambda = thetaT[4], mu = thetaT[5])
    } else {
        x <- rCGMY(n = size, C = theta[1], G = theta[2], G = theta[3],
                   Y = theta[4])
    }
    Estim <- getTempEstimation_parallel(thetaT = thetaT, x = x, size = size,
                                        Ncol = Ncol,
                                        TemperedType = TemperedType,
                                        Estimfct = Estimfct,
                                        HandleError = HandleError,
                                        eps = eps, ...)
    Output <- c(MCparam, Estim)

    return(Output)
}

#' No export.
getTempEstimation_parallel <- function(thetaT, x, size, Ncol, TemperedType,
                                       Estimfct, HandleError, eps, ...) {
    output <- vector(length = Ncol)
    if (TemperedType == "Classic") {
        output[1:7] <- c(thetaT, size)
    } else if (TemperedType == "Subordinator") {
        output[1:4] <- c(thetaT, size)
    } else if (TemperedType == "Normal") {
        output[1:6] <- c(thetaT, size)
    } else {
        output[1:5] <- c(thetaT, size)
    }
    theta0 <- thetaT - 0.1  #noise
    EstimRes <- TemperedEstim_v2(TemperedType = TemperedType,
                                 EstimMethod = Estimfct, data = x,
                                 theta0 = theta0, ComputeCov = FALSE,
                                 HandleError = HandleError, eps = eps, ...)
    if (TemperedType == "Classic") {
        output[8:13] <- EstimRes$par
    } else if (TemperedType == "Subordinator") {
        output[5:7] <- EstimRes$par
    } else if (TemperedType == "Normal") {
        output[7:11] <- EstimRes$par
    } else {
        output[6:9] <- EstimRes$par
    }
    if (TemperedType == "Classic") {
        output[14:15] <- c(EstimRes$failure, EstimRes$duration)
    } else if (TemperedType == "Subordinator") {
        output[8:9] <- c(EstimRes$failure, EstimRes$duration)
    } else if (TemperedType == "Normal") {
        output[12:13] <- c(EstimRes$failure, EstimRes$duration)
    } else {
        output[10:11] <- c(EstimRes$failure, EstimRes$duration)
    }
    return(output)
}


##### for statistical summary#####

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param EstimOutput A gap holder.
#' @param FctsToApply A gap holder.
#' @param SampleSizes A gap holder.
#' @param CheckMat A gap holder.
#' @param tolFailCheck A gap holder.
#' @param MCparam A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
ComputeStatOutput <- function(EstimOutput, FctsToApply, SampleSizes, CheckMat,
                              tolFailCheck, MCparam, ...) {
    list(alpha = ComputeStatOutputPar(EstimOutput = EstimOutput,
                                      FctsToApply = FctsToApply, par = "alpha",
                                      SampleSizes = SampleSizes,
                                      CheckMat = CheckMat,
                                      tolFailCheck = tolFailCheck,
                                      MCparam = MCparam, ...),
         beta = ComputeStatOutputPar(EstimOutput = EstimOutput,
                                     FctsToApply = FctsToApply, par = "beta",
                                     SampleSizes = SampleSizes,
                                     CheckMat = CheckMat,
                                     tolFailCheck = tolFailCheck,
                                     MCparam = MCparam, ...),
         gamma = ComputeStatOutputPar(EstimOutput = EstimOutput,
                                      FctsToApply = FctsToApply, par = "gamma",
                                      SampleSizes = SampleSizes,
                                      CheckMat = CheckMat,
                                      tolFailCheck = tolFailCheck,
                                      MCparam = MCparam, ...),
         delta = ComputeStatOutputPar(EstimOutput = EstimOutput,
                                      FctsToApply = FctsToApply, par = "delta",
                                      SampleSizes = SampleSizes,
                                      CheckMat = CheckMat,
                                      tolFailCheck = tolFailCheck,
                                      MCparam = MCparam, ...))
}


##### for Output File#####

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param thetaT A gap holder.
#' @param MCparam A gap holder.
#' @param TemperedType A gap holder.
#' @param Estimfct A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
initOutputFile <- function(thetaT, MCparam, TemperedType, Estimfct, ...) {
    method <- Estim_Des_Temp(TemperedType, Estimfct, ...)
    #TODO: get_filename(...) is not defined
    fileName <- get_filename(thetaT, MCparam, method)
    if (!file.exists(fileName)) {
        write(x = paste("alphaT", "delta+T", "delta-T", "lambda+T", "lambda-T",
                        "muT", "data size", "seed", "alphaE", "delta+E",
                        "delta-E", "lambda+E", "lambda-E", "muE", "failure",
                        "time", sep = ","),
              file = fileName, sep = "\n")
    }
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param TemperedType A String. Either "Classic", "Subordinator", or "Normal".
#' @param EstimMethod A String. Either "ML", "GMM", "Cgmm", or "Kout".
#'
#' @return Gap holder for return.
#'
#' @export
Estim_Des_Temp <- function(TemperedType = c("Classic", "Subordinator"),
                           EstimMethod = c("ML", "GMM", "Cgmm", "Kout"), ...) {
    type <- match.arg(TemperedType)
    method <- match.arg(EstimMethod)
    EstimFcts <- getTempEstimFcts(type, method)
    EstimFcts$methodDes(...)
}

##### for Checkpoints#####

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param ParameterMatrix A gap holder.
#' @param TemperedType A gap holder.
#' @param Estimfct A gap holder.
#' @param nab A gap holder.
#' @param npar A gap holder.
#' @param nSS A gap holder.
#' @param MCparam A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
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
        tab <- as.numeric(read.table(file = fileName, header = F, sep = ";"))
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

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param ParameterMatrix A gap holder.
#' @param nab A gap holder.
#' @param npar A gap holder.
#' @param MCparam A gap holder.
#' @param method A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
get_filename_checkPoint_Temp <- function(ParameterMatrix, nab, npar, MCparam,
                                         method) {
    if (npar == 3) {
        MC <- paste(paste("alpha0=", ParameterMatrix[1, 1], sep = ""),
                    paste("delta0=", ParameterMatrix[1, 2], sep = ""),
                    paste("lambda0=", ParameterMatrix[1, 3], sep = ""),
                    paste("alphan=", ParameterMatrix[nab, 1], sep = ""),
                    paste("deltan=",ParameterMatrix[nab, 2], sep = ""),
                    paste("lambdan=", ParameterMatrix[nab, 3], sep = ""),
                    paste("MCparam",MCparam, sep = ""), sep = "_")
    } else {
        MC <- paste(paste("alpha0=", ParameterMatrix[1, 1], sep = ""),
                    paste("delta+0=", ParameterMatrix[1, 2], sep = ""),
                    paste("delta-0=", ParameterMatrix[1, 3], sep = ""),
                    paste("lambda+0=", ParameterMatrix[1, 4], sep = ""),
                    paste("lambda-0=",ParameterMatrix[1, 5], sep = ""),
                    paste("mu0=", ParameterMatrix[1, 6], sep = ""),
                    paste("alphan=", ParameterMatrix[nab,1], sep = ""),
                    paste("delta+n=", ParameterMatrix[nab, 2], sep = ""),
                    paste("delta-n=", ParameterMatrix[nab,3], sep = ""),
                    paste("lambda+n=", ParameterMatrix[nab, 4], sep = ""),
                    paste("lambda-n=", ParameterMatrix[nab,5], sep = ""),
                    paste("mun=", ParameterMatrix[nab, 5], sep = ""),
                    paste("MCparam", MCparam, sep = ""), sep = "_")
    }
    fileName <- paste(MC, method, "_CHECKPOINT.txt", sep = "")
    fileName
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param CheckPointValues A gap holder.
#' @param MCparam A gap holder.
#' @param lS A gap holder.
#' @param nab A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
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

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param ParameterMatrix A gap holder.
#' @param TemperedType A gap holder.
#' @param Estimfct A gap holder.
#' @param nab A gap holder.
#' @param npar A gap holder.
#' @param nSS A gap holder.
#' @param MCparam A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
deleteCheckPoint <- function(ParameterMatrix, TemperedType, Estimfct, nab, npar,
                             nSS, MCparam, ...) {
    method <- Estim_Des_Temp(TemperedType, Estimfct, ...)
    fileName <- get_filename_checkPoint_Temp(ParameterMatrix, nab, npar,
                                             MCparam, method)
    unlink(x = fileName, force = TRUE)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param ParameterMatrix A gap holder.
#' @param Estimfct A gap holder.
#' @param ab A gap holder.
#' @param nab A gap holder.
#' @param npar A gap holder.
#' @param sample A gap holder.
#' @param nSS A gap holder.
#' @param mc A gap holder.
#' @param MCparam A gap holder.
#'
#' @return Gap holder for return.
#'
#' @export
writeCheckPoint <- function(ParameterMatrix, Estimfct, ab, nab, npar, sample,
                            nSS, mc, MCparam, ...) {
    method <- Estim_Des_Temp(TemperedType, Estimfct, ...)
    fileName <- get_filename_checkPoint_Temp(ParameterMatrix, nab, npar,
                                             MCparam, method)
    line = readLines(fileName, -1)
    line[2] = paste("--", ab, nab, npar, sample, nSS, mc, MCparam, sep = ";")
    writeLines(line, fileName)
}
