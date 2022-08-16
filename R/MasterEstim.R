##### Master function#####

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param TemperedType A String. Either "Classic", "Subordinator", "Normal", or
#' "CGMY".
#' @param EstimMethod A String. Either "ML", "GMM", "Cgmm", or "GMC".
#' @param data A gap holder.
#' @param theta0 A gap holder. \code{NULL} by default.
#' @param ComputeCov A Boolean. \code{FALSE} by default.
#' @param HandleError A Boolean. \code{TRUE} by default.
#' @param eps A gap holder. \code{1e-06} by default.
#' @param ... A gap holder.
#'
#'
#' @return Gap holder for return.
#'
#' @examples
#' testData <- c(1.8873152, -0.4843727,  0.4755897, -0.1257814,  1.3484823,
#'              -0.3866821, -0.4258380, -0.4658479, -2.9774961,  0.9646364,
#'              -0.5875601, -2.0316790,  0.3641900,  1.1882307,  1.6635770,
#'              -0.0554876,  0.4005471,  0.7820444, -0.3786902,  1.5131663)
#' TemperedEstim("Classic","ML",testData)
#'
#' @export
TemperedEstim <- function(TemperedType = c("Classic", "Subordinator", "Normal",
                                           "CGMY"),
                          EstimMethod = c("ML", "GMM", "Cgmm", "GMC"), data,
                          theta0 = NULL, ComputeCov = FALSE, HandleError = TRUE,
                          eps = 1e-06, ...) {
    if (missing(data))
        stop("data not provided !")
    if (is.null(theta0)) {
        if (TemperedType == "Classic") {
            theta0 <- MoC_CTS(x <- data, c(1.5, 1, 1, 1, 1, 0), eps = eps)
        } else if (TemperedType == "Subordinator") {
            theta0 <- MoC_STS(x <- data, c(0.5, 1, 1), eps = eps)
        } else if (TemperedType == "Normal") {
            theta0 <- MoC_NTS(x <- data, c(0.5, 0, 1, 1, 0), eps = eps)
        } else {
            theta0 <- MoC_CGMY(x <- data, c(1, 1, 1, 1.5), eps = eps)
        }
    }
    if (TemperedType == "Classic") {
        OutputObj <- methods::new(Class="EstimClassicClass",par = numeric(6),
                                  par0 = theta0, vcov = matrix(0, 6, 6),
                                  confint = matrix(0, 6,2), data = data,
                                  failure = 1)
    } else if (TemperedType == "Subordinator") {
        OutputObj <- methods::new(Class = "EstimSubClass", par = numeric(3),
                                  par0 = theta0, vcov = matrix(0, 3, 3),
                                  confint = matrix(0, 3, 2), data = data,
                                  failure = 1)
    } else if (TemperedType == "Normal") {
        OutputObj <- methods::new(Class = "EstimNormalClass", par = numeric(5),
                                  par0 = theta0, vcov = matrix(0, 5, 5),
                                  confint = matrix(0, 5, 2), data = data,
                                  failure = 1)
    } else {
        OutputObj <- methods::new(Class = "EstimCGMYClass", par = numeric(4),
                                  par0 = theta0, vcov = matrix(0, 4, 4),
                                  confint = matrix(0, 4, 2), data = data,
                                  failure = 1)
    }
    type <- match.arg(TemperedType)
    method <- match.arg(EstimMethod)
    EstimFcts <- getTempEstimFcts(type, method)
    res <- .initResTemp(type, method)
    if (HandleError) {
        tr <- tryCatch(EstimFcts$Params(x = data, theta0 = theta0, ...),
                       error = function(e) e)
        err <- inherits(tr, "error")
        if (!err) {
            res <- tr
            OutputObj@failure <- 0
        }
    } else {
        res <- EstimFcts$Params(x = data, theta0 = theta0, ...)
        OutputObj@failure <- 0
    }
    OutputObj@par <- NameParamsObjectsTemp(res$Estim$par, type)
    OutputObj@others <- res$Estim
    OutputObj@duration <- as.numeric(res$duration)
    OutputObj@method <- res$method
    if (ComputeCov) {
        OutputObj@vcov <- EstimFcts$CovarianceMat(data = OutputObj@data,
                                                  EstimObj = res, ...)
        OutputObj@confint <- AsymptoticConfidenceInterval(
          thetaEst = OutputObj@par, n_sample = OutputObj@sampleSize,
          Cov = OutputObj@vcov, qLaw = stats::qnorm, type = type, ...)
    }
    OutputObj
}

# Function title
#
# Gap holder for description.
#
# Gap holder for details.
#
# @param TemperedType A String. Either "Classic", "Subordinator", "Normal", or
# "CGMY".
# @param EstimMethod A String. Either "ML", "GMM", "Cgmm", or "GMC".
# @param data A gap holder.
# @param theta0 A gap holder. \code{NULL} by default.
# @param ComputeCov A Boolean. \code{FALSE} by default.
# @param HandleError A Boolean. \code{TRUE} by default.
# @param eps A gap holder. \code{1e-06} by default.
#
# @return Gap holder for return.
#
# @export
# TemperedEstim_v2 <- function(TemperedType = c("Classic", "Subordinator",
#                                               "Normal", "CGMY"),
#                              EstimMethod = c("ML", "GMM", "Cgmm", "GMC"), data,
#                              theta0 = NULL, ComputeCov = FALSE,
#                              HandleError = TRUE, eps = 1e-06, ...) {
#     if (missing(data))
#         stop("data not provided !")
#     if (is.null(theta0)) {
#         if (TemperedType == "Classic") {
#             theta0 <- MoC_CTS(x = data, c(1.5, 1, 1, 1, 1, 0), eps = eps)
#         } else if (TemperedType == "Subordinator") {
#             theta0 <- MoC_STS(x = data, c(0.5, 1, 1), eps = eps)
#         } else if (TemperedType == "Normal") {
#             theta0 <- MoC_NTS(x = data, c(0.5, 0, 1, 1, 0), eps = eps)
#         } else {
#             theta0 <- MoC_CGMY(x = data, c(1, 1, 1, 1.5), eps = eps)
#         }
#     }
#     if (TemperedType == "Classic") {
#         OutputObj <- list(par = numeric(6), par0 = theta0,
#                           vcov = matrix(0, 6, 6), confint = matrix(0, 6, 2),
#                           data = data, failure = 1)
#     } else if (TemperedType == "Subordinator") {
#         OutputObj <- list(par = numeric(3), par0 = theta0,
#                           vcov = matrix(0, 3, 3), confint = matrix(0, 3, 2),
#                           data = data, failure = 1)
#     } else if (TemperedType == "Normal") {
#         OutputObj <- list(par = numeric(5), par0 = theta0,
#                           vcov = matrix(0, 5, 5), confint = matrix(0, 5, 2),
#                           data = data, failure = 1)
#     } else {
#         OutputObj <- list(par = numeric(4), par0 = theta0,
#                           vcov = matrix(0, 4, 4), confint = matrix(0, 4, 2),
#                           data = data, failure = 1)
#     }
#     type <- match.arg(TemperedType)
#     method <- match.arg(EstimMethod)
#     EstimFcts <- getTempEstimFcts(type, method)
#     res <- .initResTemp(type, method)
#     if (HandleError) {
#         tr <- tryCatch(EstimFcts$Params(x = data, theta0 = theta0, ...),
#                        error = function(e) e)
#         err <- inherits(tr, "error")
#         if (!err) {
#             res <- tr
#             OutputObj$failure <- 0
#         }
#     } else {
#         res <- EstimFcts$Params(x = data, theta0 = theta0, ...)
#         OutputObj$failure <- 0
#     }
#     OutputObj$par <- NameParamsObjectsTemp(res$Estim$par, type)
#     OutputObj$others <- res$Estim
#     OutputObj$duration <- as.numeric(res$duration)
#     OutputObj$method <- res$method
#     if (ComputeCov) {
#         OutputObj$vcov <- EstimFcts$CovarianceMat(data = OutputObj$data,
#                                                   EstimObj = res, ...)
#         OutputObj$confint <- AsymptoticConfidenceInterval(
#           thetaEst = OutputObj$par, n_sample = OutputObj$sampleSize,
#           Cov = OutputObj$vcov, qLaw = qnorm, type = type, ...)
#     }
#     OutputObj
# }

##### auxiliaries#####

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param type A String. Either "Classic", "Subordinator", "Normal", or
#' "CGMY".
#' @param method A String. Either "ML", "GMM", "Cgmm", or "GMC".
#'
#' @return Gap holder for return.
#'
#' @export
getTempEstimFcts <- function(
    type = c("Classic", "Subordinator", "Normal", "CGMY"),
    method = c("ML", "GMM", "Cgmm", "GMC")) {
    if (type == "Classic") {
        Output <- switch(method, ML = {
            list(Params = MLParametersEstim_CTS,
                 CovarianceMat = .asymptoticVarianceEstimML_CTS,
                 methodDes = .methodDesML_CTS)
        }, GMM = {
            list(Params = GMMParametersEstim_CTS,
                 CovarianceMat = .asymptoticVarianceEstimGMM_CTS,
                 methodDes = getGMMmethodName_CTS)
        }, Cgmm = {
            list(Params = CgmmParametersEstim_CTS,
                 CovarianceMat = .asymptoticVarianceEstimCgmm_CTS,
                 methodDes = getCgmmMethodName_CTS)
        }, GMC = {
            list(Params = GMCParametersEstim_CTS
                 , CovarianceMat = .asymptoticVarianceEstimGMC_CTS,
                 methodDes = getGMCmethodName_CTS)
        }, stop(paste(method, " not taken into account !")))
        Output
    } else if (type == "Subordinator") {
        Output <- switch(method, ML = {
            list(Params = MLParametersEstim_STS,
                 CovarianceMat = .asymptoticVarianceEstimML_STS,
                 methodDes = .methodDesML_STS)
        }, GMM = {
            list(Params = GMMParametersEstim_STS,
                 CovarianceMat = .asymptoticVarianceEstimGMM_STS,
                 methodDes = getGMMmethodName_STS)
        }, Cgmm = {
            list(Params = CgmmParametersEstim_STS,
                 CovarianceMat = .asymptoticVarianceEstimCgmm_STS,
                 methodDes = getCgmmMethodName_STS)
        }, GMC = {
            list(Params = GMCParametersEstim_STS,
                 CovarianceMat = .asymptoticVarianceEstimGMC_STS,
                 methodDes = getGMCmethodName_STS)
        }, stop(paste(method, " not taken into account !")))
        Output
    } else if (type == "Normal") {
        Output <- switch(method, ML = {
            list(Params = MLParametersEstim_NTS,
                 CovarianceMat = .asymptoticVarianceEstimML_NTS,
                 methodDes = .methodDesML_NTS)
        }, GMM = {
            list(Params = GMMParametersEstim_NTS,
                 CovarianceMat = .asymptoticVarianceEstimGMM_NTS,
                 methodDes = getGMMmethodName_NTS)
        }, Cgmm = {
            list(Params = CgmmParametersEstim_NTS,
                 CovarianceMat = .asymptoticVarianceEstimCgmm_NTS,
                 methodDes = getCgmmMethodName_NTS)
        }, stop(paste(method, " not taken into account !")))
        Output
    } else {
        Output <- switch(method, ML = {
            list(Params = MLParametersEstim_CGMY,
                 CovarianceMat = .asymptoticVarianceEstimML_CGMY,
                 methodDes = .methodDesML_CGMY)
        }, GMM = {
            list(Params = GMMParametersEstim_CGMY,
                 CovarianceMat = .asymptoticVarianceEstimGMM_CGMY,
                 methodDes = getGMMmethodName_CGMY)
        }, Cgmm = {
            list(Params = CgmmParametersEstim_CGMY,
                 CovarianceMat = .asymptoticVarianceEstimCgmm_CGMY,
                 methodDes = getCgmmMethodName_CGMY)
        }, GMC = {
            list(Params = GMCParametersEstim_CGMY,
                 CovarianceMat = .asymptoticVarianceEstimGMC_CGMY,
                 methodDes = getGMCmethodName_CGMY)
        }, stop(paste(method, " not taken into account !")))
        Output
    }
}

# No export.
.initResTemp <- function(type, method) {
    if (type == "Classic") {
        npar <- 6
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    } else if (type == "Subordinator") {
        npar <- 3
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    } else if (type == "Normal") {
        npar <- 5
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    } else {
        npar <- 4
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    }
}

# No export.
NameParamsObjectsTemp <- function(mat, type = c("Classic", "Subordinator",
                                                "Normal")) {
    if (type == "Classic") {
        parNames <- c("alpha", "delta +", "delta -", "lambda +", "lambda -",
                      "mu")
        minMaxCol <- c("min", "max")
        if (length(mat) == 6) {
            names(mat) <- parNames
        } else if (is.matrix(mat) && nrow(mat) == 6) {
            rownames(mat) <- parNames
            if (ncol(mat) == 2)
                colnames(mat) <- minMaxCol else if (ncol(mat) == 6)
                colnames(mat) <- parNames
        }
    } else if (type == "Subordinator") {
        parNames <- c("alpha", "delta", "lambda")
        minMaxCol <- c("min", "max")
        if (length(mat) == 3) {
            names(mat) <- parNames
        } else if (is.matrix(mat) && nrow(mat) == 3) {
            rownames(mat) <- parNames
            if (ncol(mat) == 2)
                colnames(mat) <- minMaxCol else if (ncol(mat) == 3)
                colnames(mat) <- parNames
        }
    } else if (type == "Normal") {
        parNames <- c("alpha", "beta", "delta", "lambda", "mu")
        minMaxCol <- c("min", "max")
        if (length(mat) == 5) {
            names(mat) <- parNames
        } else if (is.matrix(mat) && nrow(mat) == 5) {
            rownames(mat) <- parNames
            if (ncol(mat) == 2)
                colnames(mat) <- minMaxCol else if (ncol(mat) == 5)
                colnames(mat) <- parNames
        }
    } else {
        parNames <- c("C", "G", "M", "Y")
        minMaxCol <- c("min", "max")
        if (length(mat) == 4) {
            names(mat) <- parNames
        } else if (is.matrix(mat) && nrow(mat) == 4) {
            rownames(mat) <- parNames
            if (ncol(mat) == 2)
                colnames(mat) <- minMaxCol else if (ncol(mat) == 4)
                colnames(mat) <- parNames
        }

    }
    mat
}

# No Export.
CheckParametersRange_STS <- function(theta) {
    alpha <- theta[1]
    delta <- theta[2]
    lambda <- theta[3]
    checkParams <- list(alpha = checkRange(alpha, 0, 1, "alpha"),
                        delta = checkRange(delta, 0, Inf, "delta"),
                        lambda = checkRange(lambda, 0, Inf,
                                            ParamName = "lambda"))
    .printErr <- function(errList) if (!errList$bool)
        stop(errList$msg)
    lapply(checkParams, .printErr)
}

# No Export.
CheckParametersRange_CTS <- function(theta) {
    alpha <- theta[1]
    deltap <- theta[2]
    deltam <- theta[3]
    lambdap <- theta[4]
    lambdam <- theta[5]
    mu <- theta[6]
    checkParams <- list(alpha = checkRange(alpha, 0, 2, "alpha"),
                        deltap = checkRange(deltap, 0, Inf, "delta+"),
                        deltam = checkRange(deltam, 0, Inf, "delta-"),
                        lambdap = checkRange(lambdap, 0, Inf, "lambda+"),
                        lambdam = checkRange(lambdam, 0, Inf, "lambda-"),
                        mu = checkRange(mu, -Inf, Inf, "mu"))
    .printErr <- function(errList) if (!errList$bool)
        stop(errList$msg)
    lapply(checkParams, .printErr)
}

# No Export.
CheckParametersRange_NTS <- function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    delta <- theta[3]
    lambda <- theta[4]
    mu <- theta[5]
    checkParams <- list(alpha = checkRange(alpha, 0, 2, "alpha"),
                        beta = checkRange(beta, -Inf, Inf, "beta"),
                        delta = checkRange(delta, 0, Inf, "delta"),
                        lambda = checkRange(lambda, 0, Inf, "lambda"),
                        mu = checkRange(mu, -Inf, Inf, "mu"))
    .printErr <- function(errList) if (!errList$bool)
        stop(errList$msg)
    lapply(checkParams, .printErr)
}

# No Export.
CheckParametersRange_CGMY <- function(theta) {
    C <- theta[1]
    G <- theta[2]
    M <- theta[3]
    Y <- theta[4]
    checkParams <- list(C = checkRange(C, 0, Inf, "C"),
                        G = checkRange(G, 0, Inf, "G"),
                        M = checkRange(M, 0, Inf, "M"),
                        Y = checkRange(Y, 0, 2, "Y"))
    .printErr <- function(errList) if (!errList$bool)
        stop(errList$msg)
    lapply(checkParams, .printErr)
}

# No export.
checkRange <- function(Parameter, min = -Inf, max = Inf, ParamName) {
  if ((Parameter >= min) && (Parameter <= max))
    return(list(bool = TRUE, msg = "valid"))
  else return(list(bool = FALSE,
                   msg = paste(ParamName, " = ", Parameter,
                               " sould be in the interval [", min, max, "]")))
}

##### Asymptotic Confidence Interval#####

# No export.
AsymptoticConfidenceInterval <- function(thetaEst, n_sample, Cov,
                                         qLaw = stats::qnorm, level = 0.95,
                                         type, ...) {
    if (type == "Classic") {
        nr <- 6
    } else if (type == "Subordinator") {
        nr <- 3
    } else if (type == "Normal") {
        nr <- 5
    } else {
        nr <- 4
    }
    mat <- matrix(NaN, ncol = 2, nrow = nr)
    attr(mat, "level") <- level
    z <- qLaw(level)
    V <- sqrt(diag(Cov)) * z
    mat[, 1] <- thetaEst - V
    mat[, 2] <- thetaEst + V
    NameParamsObjectsTemp(mat, type = type)
}


# No export.
# Added by Cedric 20220811
NameParamsObjects <- function(mat, type = NULL) {

  parNames <- c("alpha", "beta", "gamma", "delta")

  if (!is.null(type)){
    outputString <- switch(type,
                           Classic = c("Alpha", "DeltaP", "DeltaM", "LambdaP",
                                       "LambdaM", "mu"),
                           Subordinator = c("Alpha=", "Delta", "Lambda"),
                           Normal = c("Alpha", "Beta", "Delta", "Lambda",
                                      "mu"),
                           CGMY = c("C", "G", "M", "Y"))
  }
  else {
    if(length(mat) == 6) parNames = c("Alpha", "DeltaP", "DeltaM", "LambdaP",
                                      "LambdaM", "mu")
    else if(length(mat) == 3) parNames = c("Alpha=", "Delta", "Lambda")
    else if(length(mat) == 5) parNames  = c("Alpha", "Beta", "Delta", "Lambda",
                                            "mu")
    else if(length(mat) == 4) parNames = c("C", "G", "M", "Y")
  }

  minMaxCol <- c("min", "max")

  if (length(mat) > 2 && length(mat) < 7) {
    names(mat) <- parNames
  }
  else if (is.matrix(mat) && nrow(mat) > 2 && is.matrix(mat) && nrow(mat) < 7) {
    rownames(mat) <- parNames
    if (ncol(mat) == 2)
      colnames(mat) <- minMaxCol
    else if (ncol(mat) > 2 && ncol(mat) < 7)
      colnames(mat) <- parNames
  }
  mat
}


##### Classes#####

# No export.
EstimSubClass <- setClass("EstimSubClass",
                          slots = list(par = "numeric", par0 = "numeric",
                                       vcov = "matrix", confint = "matrix",
                                       data = "numeric", sampleSize = "numeric",
                                       others = "list", duration = "numeric",
                                       failure = "numeric",
                                       method = "character"),
                          contains = list(), validity = function(object) {
        par <- object@par
        if (length(par) == 3)
            ansp <- TRUE
        else ansp <- "Parameter of length different of 3"
        par0 <- object@par0
        if (length(par0) == 3)
            ansp0 <- TRUE
        else ansp0 <- "Initial Parameter of length different of 3"
        vcov <- object@vcov
        if (ncol(vcov) == 3 && nrow(vcov) == 3)
            anscov <- TRUE
        else anscov <- "covariance matrix of length different of 3x3"
        confint <- object@confint
        if (ncol(confint) == 2 && nrow(confint) == 3)
            ansconfint <- TRUE
        else ansconfint <-
          "confidance intervall matrix of length different of 3x2"
        if (ansp == TRUE && ansp0 == TRUE && anscov == TRUE &&
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
    })

# No export.
EstimClassicClass <- setClass("EstimClassicClass",
                              slots = list(par = "numeric", par0 = "numeric",
                                           vcov = "matrix", confint = "matrix",
                                           data = "numeric",
                                           sampleSize = "numeric",
                                           others = "list",
                                           duration = "numeric",
                                           failure = "numeric",
                                           method = "character"),
                              contains = list(), validity = function(object) {
        par <- object@par
        if (length(par) == 6)
            ansp <- TRUE
        else ansp <- "Parameter of length different of 6"
        par0 <- object@par0
        if (length(par0) == 6)
            ansp0 <- TRUE
        else ansp0 <- "Initial Parameter of length different of 6"
        vcov <- object@vcov
        if (ncol(vcov) == 6 && nrow(vcov) == 6)
            anscov <- TRUE
        else anscov <- "covariance matrix of length different of 6x6"
        confint <- object@confint
        if (ncol(confint) == 2 && nrow(confint) == 6)
            ansconfint <- TRUE
        else ansconfint <-
          "confidance intervall matrix of length different of 6x2"
        if (ansp == TRUE && ansp0 == TRUE && anscov == TRUE &&
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
    })

# No export.
EstimNormalClass <- setClass("EstimNormalClass",
                             slots = list(par = "numeric", par0 = "numeric",
                                          vcov = "matrix", confint = "matrix",
                                          data = "numeric",
                                          sampleSize = "numeric",
                                          others = "list", duration = "numeric",
                                          failure = "numeric",
                                          method = "character"),
                             contains = list(), validity = function(object) {
        par <- object@par
        if (length(par) == 5)
            ansp <- TRUE
        else ansp <- "Parameter of length different of 5"
        par0 <- object@par0
        if (length(par0) == 5)
            ansp0 <- TRUE
        else ansp0 <- "Initial Parameter of length different of 5"
        vcov <- object@vcov
        if (ncol(vcov) == 5 && nrow(vcov) == 5)
            anscov <- TRUE
        else anscov <- "covariance matrix of length different of 5x5"
        confint <- object@confint
        if (ncol(confint) == 2 && nrow(confint) == 5)
            ansconfint <- TRUE
        else ansconfint <-
          "confidance intervall matrix of length different of 5x2"
        if (ansp == TRUE && ansp0 == TRUE && anscov == TRUE &&
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
    })

# No export.
EstimCGMYClass <- setClass("EstimCGMYClass",
                           slots = list(par = "numeric", par0 = "numeric",
                                        vcov = "matrix", confint = "matrix",
                                        data = "numeric",
                                        sampleSize = "numeric", others = "list",
                                        duration = "numeric",
                                        failure = "numeric",
                                        method = "character"),
                           contains = list(), validity = function(object) {
        par <- object@par
        if (length(par) == 4)
            ansp <- TRUE
        else ansp <- "Parameter of length different of 4"
        par0 <- object@par0
        if (length(par0) == 4)
            ansp0 <- TRUE
        else ansp0 <- "Initial Parameter of length different of 4"
        vcov <- object@vcov
        if (ncol(vcov) == 4 && nrow(vcov) == 4)
            anscov <- TRUE
        else anscov <- "covariance matrix of length different of 4x4"
        confint <- object@confint
        if (ncol(confint) == 2 && nrow(confint) == 4)
            ansconfint <- TRUE
        else ansconfint <-
          "confidance intervall matrix of length different of 4x2"
        if (ansp == TRUE && ansp0 == TRUE && anscov == TRUE &&
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
    })
