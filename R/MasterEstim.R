##### Master function#####

#' Estimation function
#'
#' Main estimation function for the tempered stable subordinator distribution
#' (TSS), the classical tempered stable distribution (CTS), and the normal
#' tempered stable distribution (NTS). It allows the user to select the
#' preferred estimation method and several related options.
#'
#' \strong{TemperedType} Detailed documentation of the individual tempered
#' stable distributions can be viewed in the respective characteristic function.
#' Use [charTSS()], [charCTS()], or [charNTS()].
#'
#' \strong{Estimfct} Additional parameters are needed for different estimation
#' functions. These are listed below for each function. The list of additional
#' parameters starts after the parameter \code{eps} in the parameter list.
#' \describe{
#'   \item{For ML:}{ See usage of Maximum likelihood estimation in Kim et al.
#'   (2008). No additional parameters are needed.}
#'   \item{For GMM:}{Generalized Method of Moments by Feuerverger (1981).
#'   The parameters \code{algo, alphaReg, regularization, WeightingMatrix, and
#'   t_scheme} must be specified.
#'
#'   Parameter \code{t_scheme}: One of the most important features of this
#'   method is that it allows the user to choose how to place the points where
#'   the moment conditions are evaluated. One can choose among 6 different
#'   options. Depending on the option, further parameters have to be passed.
#'   \describe{
#'     \item{"equally":}{equally placed points in [min_t,max_t]. When provided,
#'      user's \code{min_t} and \code{max_t} will be used (when
#'      \code{Coinstrained == FALSE}).
#'     }
#'     \item{"NonOptAr":}{non optimal arithmetic placement.
#'     }
#'     \item{"uniformOpt":}{uniform optimal placement.
#'     }
#'     \item{"ArithOpt":}{arithmetic optimal placement.
#'     }
#'     \item{"Var Opt":}{optimal variance placement as explained above.
#'     }
#'     \item{"free":}{user needs to pass his own set of points in \code{t_free}.
#'     }
#'   }
#'
#'   Parameter \code{WeightingMatrix}: One can choose among 3 different options:
#'   \describe{
#'     \item{"OptAsym":}{the optimal asymptotic choice.
#'     }
#'     \item{"DataVar":}{the covariance matrix of the data provided.
#'     }
#'     \item{"Id":}{the identity matrix.
#'     }
#'   }
#'   }
#'   \item{For Cgmm:}{Continuum Generalized Methods of Moments by Carrasco &
#'   Kotchoni (2017). The parameters \code{algo, alphaReg, subdivisions,
#'   IntegrationMethod, randomIntegrationLaw, s_min, and s_max} must be
#'   specified.
#'   }
#'   \item{For GMC:}{Generalized Method of Cumulants (GMC) by Massing, T.
#'    (2022). The parameters \code{algo, alphaReg, regularization,
#'    WeightingMatrix, and ncond} must be specified.
#'    }
#' }
#'
#' \strong{Estim-Class} Class storing all the information about the estimation
#' method; output of this function.
#' \strong{Slots of the return class}
#' \describe{
#'   \item{par:}{Object of class "\code{numeric}"; Value of the estimated
#'   parameters.}
#'   \item{par0:}{Object of class "\code{numeric}"; Initial guess for the
#'   parameters.}
#'   \item{vcov:}{Object of class "\code{matrix}" representing the covariance
#'   matrix.}
#'   \item{confint:}{Object of class "\code{matrix}" representing the confidence
#'   interval computed at a specific level (attribute of the object).}
#'   \item{data:}{Object of class "\code{numeric}" used to compute the
#'   estimation.}
#'   \item{sampleSize:}{Object of class "\code{numeric}" ; length of the data.}
#'   \item{others:}{Object of class "\code{list}" ; more information about the
#'   estimation method.}
#'   \item{duration:}{Object of class "\code{numeric}" ; duration in seconds.}
#'   \item{failure:}{Object of class "\code{numeric}" representing the status of
#'    the procedure: 0 failure or 1 success.}
#'   \item{method:}{Object of class "\code{character}" description of the
#'   parameter used in the estimation.}
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
#' \doi{10.1016/j.jbankfin.2007.11.004};
#'
#' Hansen, L. P. (1982), 'Large sample properties of generalized method of
#' moments estimators' \doi{10.2307/1912775};
#'
#' Hansen, L. P.; Heaton, J. & Yaron, A. (1996), 'Finite-Sample Properties of
#' Some Alternative GMM Estimators' \doi{10.1080/07350015.1996.10524656}
#'
#' Carrasco, M. & Kotchoni, R. (2017), 'Efficient estimation using the
#' characteristic function' \doi{10.1017/S0266466616000025};
#'
#' Kuechler, U. & Tappe, S. (2013), 'Tempered stable distribution and processes'
#' \doi{10.1016/j.spa.2013.06.012};
#'
#' Feuerverger, A. & McDunnough, P. (1981), 'On the efficiency of empirical
#' characteristic function procedures'
#' \doi{10.1111/j.2517-6161.1981.tb01143.x};
#'
#' @param TemperedType A String. Either "Classic", "Subordinator", or "Normal"
#' @param EstimMethod A String. Either "ML", "GMM", "Cgmm", or "GMC".
#' @param data Data used to perform the estimation: numeric vector of length n.
#' @param theta0 A vector of numeric values corresponding to the pattern of the
#' \code{TemperedType}.
#' @param ComputeCov 	Logical flag: If set to TRUE, the asymptotic covariance
#' matrix is computed. \code{FALSE} by default.
#' @param HandleError Logical flag: If set to \code{TRUE} and if an error occurs
#' during the estimation procedure, the computation will carry on and NA will be
#' returned. Useful for Monte Carlo simulations.\code{TRUE} by default.
#' @param eps Error tolerance. \code{1e-06} by default.
#' @param algo algorithm: For GMM: \code{"2SGMM"} is the two step GMM proposed
#' by Hansen (1982). \code{"CueGMM"} and \code{"ITGMM"} are respectively the
#' continuous updated and the iterative GMM proposed by Hansen, Eaton et Yaron
#' (1996) and adapted to the continuum case. For Cgmm: \code{"2SCgmm",
#' "CueCgmm", ...}. Same for GMC.
#' @param regularization regularization scheme to be used, one of
#' \code{"Tikhonov"} (Tikhonov), \code{"LF"} (Landweber-Fridmann) and
#' \code{"cut-off"} (spectral cut-off).
#' @param WeightingMatrix type of weighting matrix used to compute the
#' objective function, one of \code{"OptAsym"} (the optimal asymptotic),
#' \code{"DataVar"} (the data driven) and \code{"Id"} (the identity matrix).
#' @param t_scheme scheme used to select the points where the moment conditions
#' are evaluated, one of \code{"equally"} (equally placed), \code{"NonOptAr"}
#' (non optimal arithmetic placement), \code{"uniformOpt"}
#' (uniform optimal placement), \code{"ArithOpt"} (arithmetic optimal
#' placement), \code{"Var Opt"} (optimal variance placement) and \code{"free"}
#' (users need to pass their own set of points in ...).
#' @param alphaReg value of the regularisation parameter; numeric.
#' @param t_free sequence, if \code{t_scheme=="free"}.
#' @param subdivisions 	Number of subdivisions used to compute the different
#' integrals involved in the computation of the objective function (to
#' minimise); numeric.
#' @param IntegrationMethod Numerical integration method to be used to
#' approximate the (vectorial) integrals. Users can choose between "Uniform"
#' discretization or the "Simpson"'s rule (the 3-point Newton-Cotes quadrature
#' rule).
#' @param randomIntegrationLaw Probability measure associated to the Hilbert
#' space spanned by the moment conditions.
#' @param s_min,s_max Lower and Upper bounds of the interval where the moment
#' conditions are considered; numeric.
#' @param ncond TODO
#' @param ... Other arguments to be passed to the estimation function or the
#' asymptotic confidence level.
#'
#'
#' @return Object of a estim-class. See details for more information.
#'
#' @examples
#' \donttest{
#' TemperedEstim(TemperedType = "Classic", EstimMethod = "ML",
#'                data = rCTS(2,1.5,1,1,1,1,0),
#'                theta0 = c(1.5,1,1,1,1,0) - 0.1);
#' TemperedEstim("Subordinator", "GMM", rTSS(20,0.5,1,1), algo = "2SGMM",
#'               alphaReg = 0.01, regularization = "cut-off",
#'               WeightingMatrix = "OptAsym", t_scheme = "free",
#'               t_free = seq(0.1,2,length.out = 12));
#' TemperedEstim("Normal", "Cgmm", rNTS(20,0.5,1,1,1,0), algo = "2SCgmm",
#'               alphaReg = 0.01, subdivisions = 20,
#'               IntegrationMethod = "Uniform", randomIntegrationLaw = "unif",
#'               s_min = 0, s_max= 1);
#' TemperedEstim("Subordinator", "GMC", rTSS(20, 0.5, 1, 1), algo = "2SGMC",
#'               alphaReg = 0.01, WeightingMatrix = "OptAsym",
#'               regularization = "cut-off", ncond = 8);
#' }
#' @export
TemperedEstim <- function(TemperedType = c("Classic", "Subordinator", "Normal"),
                          EstimMethod = c("ML", "GMM", "Cgmm", "GMC"), data,
                          theta0 = NULL, ComputeCov = FALSE, HandleError = TRUE,
                          eps = 1e-06, algo = NULL, regularization = NULL,
                          WeightingMatrix = NULL, t_scheme = NULL,
                          alphaReg = NULL, t_free = NULL, subdivisions = NULL,
                          IntegrationMethod = NULL, randomIntegrationLaw = NULL,
                          s_min = NULL, s_max = NULL, ncond = NULL, ...) {
    if (missing(data))
        stop("data not provided !")
    if (is.null(theta0)) {
        if (TemperedType == "Classic") {
            theta0 <- MoC_CTS(x <- data, c(1.5, 1, 1, 1, 1, 0))
        } else if (TemperedType == "Subordinator") {
            theta0 <- MoC_TSS(x <- data, c(0.5, 1, 1))
        } else if (TemperedType == "Normal") {
            theta0 <- MoC_NTS(x <- data, c(0.5, 0, 1, 1, 0))
        }
      # else {
      #       theta0 <- MoC_CGMY(x <- data, c(1, 1, 1, 1.5))
      #   }
    }
    if (TemperedType == "Classic") {
        OutputObj <- methods::new(Class="EstimClassicClass",par = numeric(6),
                                  par0 = theta0, vcov = matrix(0, 6, 6),
                                  confint = matrix(0, 6, 2), data = data,
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
    }
    # else {
    #     OutputObj <- methods::new(Class = "EstimCGMYClass", par = numeric(4),
    #                               par0 = theta0, vcov = matrix(0, 4, 4),
    #                               confint = matrix(0, 4, 2), data = data,
    #                               failure = 1)
    # }
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
#             theta0 <- MoC_TSS(x = data, c(0.5, 1, 1), eps = eps)
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

# No export.
getTempEstimFcts <- function(
    type = c("Classic", "Subordinator", "Normal"),
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
            list(Params = MLParametersEstim_TSS,
                 CovarianceMat = .asymptoticVarianceEstimML_TSS,
                 methodDes = .methodDesML_TSS)
        }, GMM = {
            list(Params = GMMParametersEstim_TSS,
                 CovarianceMat = .asymptoticVarianceEstimGMM_TSS,
                 methodDes = getGMMmethodName_TSS)
        }, Cgmm = {
            list(Params = CgmmParametersEstim_TSS,
                 CovarianceMat = .asymptoticVarianceEstimCgmm_TSS,
                 methodDes = getCgmmMethodName_TSS)
        }, GMC = {
            list(Params = GMCParametersEstim_TSS,
                 CovarianceMat = .asymptoticVarianceEstimGMC_TSS,
                 methodDes = getGMCmethodName_TSS)
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
    }
  # else {
  #       Output <- switch(method, ML = {
  #           list(Params = MLParametersEstim_CGMY,
  #                CovarianceMat = .asymptoticVarianceEstimML_CGMY,
  #                methodDes = .methodDesML_CGMY)
  #       }, GMM = {
  #           list(Params = GMMParametersEstim_CGMY,
  #                CovarianceMat = .asymptoticVarianceEstimGMM_CGMY,
  #                methodDes = getGMMmethodName_CGMY)
  #       }, Cgmm = {
  #           list(Params = CgmmParametersEstim_CGMY,
  #                CovarianceMat = .asymptoticVarianceEstimCgmm_CGMY,
  #                methodDes = getCgmmMethodName_CGMY)
  #       }, GMC = {
  #           list(Params = GMCParametersEstim_CGMY,
  #                CovarianceMat = .asymptoticVarianceEstimGMC_CGMY,
  #                methodDes = getGMCmethodName_CGMY)
  #       }, stop(paste(method, " not taken into account !")))
  #       Output
  #   }
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
    }
    # else {
    #     npar <- 4
    #     list(Estim = list(par = rep(NaN, npar)), duration = 0,
    #          method = paste(type, method, "failed", sep = "_"))
    # }
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
    }
    # else {
    #     parNames <- c("C", "G", "M", "Y")
    #     minMaxCol <- c("min", "max")
    #     if (length(mat) == 4) {
    #         names(mat) <- parNames
    #     } else if (is.matrix(mat) && nrow(mat) == 4) {
    #         rownames(mat) <- parNames
    #         if (ncol(mat) == 2)
    #             colnames(mat) <- minMaxCol else if (ncol(mat) == 4)
    #             colnames(mat) <- parNames
    #     }
    #
    # }
    mat
}

# No Export.
CheckParametersRange_TSS <- function(theta) {
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
    }
    # else {
    #     nr <- 4
    # }
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
    parNames <- switch(type,
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
#' @importFrom methods new
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

## Init method

setMethod("initialize", "EstimSubClass",
          function(.Object, par, par0, vcov, confint, method, level, others,
                   data, duration, failure, ...) {
            ## handle missing
            if (missing(par))
              par        <- numeric(3)
            if (missing(par0))
              par0       <- numeric(3)
            if (missing(vcov))
              vcov       <- matrix(nrow = 3, ncol = 3)
            if (missing(confint))
              confint    <- matrix(nrow = 3, ncol = 2)
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
            NameParamsObjects(par, "Sub")
            NameParamsObjects(par0, "Sub")
            NameParamsObjects(vcov, "Sub")
            NameParamsObjects(confint, "Sub")
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

setMethod("show", "EstimSubClass",
          function(object) {
            cat("*** Tempered Estim Sub, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim Sub) ******* \n")
          })

# No export.
#' @importFrom methods new
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

## Init method

setMethod("initialize", "EstimClassicClass",
          function(.Object, par, par0, vcov, confint, method, level, others,
                   data, duration, failure, ...) {
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


# No export.
#' @importFrom methods new
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

## Init method

setMethod("initialize", "EstimNormalClass",
          function(.Object, par, par0, vcov, confint, method, level, others,
                   data, duration, failure, ...) {
            ## handle missing
            if (missing(par))
              par        <- numeric(5)
            if (missing(par0))
              par0       <- numeric(5)
            if (missing(vcov))
              vcov       <- matrix(nrow = 5, ncol = 5)
            if (missing(confint))
              confint    <- matrix(nrow = 5, ncol = 2)
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
            NameParamsObjects(par, "Normal")
            NameParamsObjects(par0, "Normal")
            NameParamsObjects(vcov, "Normal")
            NameParamsObjects(confint, "Normal")
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

setMethod("show", "EstimNormalClass",
          function(object) {
            cat("*** Tempered Estim Normal, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim Normal) ******* \n")
          })


# No export.
#' @importFrom methods new
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

## Init method

setMethod("initialize", "EstimCGMYClass",
          function(.Object, par, par0, vcov, confint, method, level, others,
                   data, duration, failure, ...) {
            ## handle missing
            if (missing(par))
              par        <- numeric(4)
            if (missing(par0))
              par0       <- numeric(4)
            if (missing(vcov))
              vcov       <- matrix(nrow = 4, ncol = 4)
            if (missing(confint))
              confint    <- matrix(nrow = 4, ncol = 2)
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
            NameParamsObjects(par, "CGMY")
            NameParamsObjects(par0, "CGMY")
            NameParamsObjects(vcov, "CGMY")
            NameParamsObjects(confint, "CGMY")
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

setMethod("show", "EstimCGMYClass",
          function(object) {
            cat("*** Tempered Estim CGMY, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim CGMY) ******* \n")
          })
