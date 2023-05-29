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
#' Use [charTSS()], [charCTS()], or [charNTS()]. TODO: For all other
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
#'     \item{"equally":}{equally placed points in \code{min_t,max_t}. When
#'     provided, user's \code{min_t} and \code{max_t} will be used (when
#'     \code{Coinstrained == FALSE}).
#'     }
#'     \item{"NonOptAr":}{non optimal arithmetic placement.
#'     }
#'     \item{"uniformOpt":}{uniform optimal placement.
#'     }
#'     \item{"ArithOpt":}{arithmetic optimal placement.
#'     }
#'     \item{"Var Opt":}{optimal variance placement as explained above.
#'     }
#'     \item{"free":}{user needs to pass own set of points in \code{t_free}.
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
#'
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
#' \strong{IterationControl} If \code{algo = "IT..."} or \code{algo =
#' "Cue..."} the user can control each iteration by setting up the list
#' IterationControl which contains the following elements:
#' \describe{
#'   \item{NbIter}{maximum number of iteration. The loop stops when NBIter is
#'   reached; default = 10.}
#'   \item{PrintIterlogical}{if set to TRUE, the value of the current parameter
#'   estimation is printed to the screen at each iteration; default = TRUE.}
#'   \item{RelativeErrMax}{the loop stops if the relative error between two
#'   consecutive estimation steps is smaller than RelativeErrMax;
#'   default = 1e-3.}
#' }
#'
#' Since this package is structurally based on the \strong{"StableEstim"
#' package by Tarak Kharrat and Georgi N. Boshnakov}, more detailed
#' documentation can be found in their documentation.
#'
#' @seealso
#' \url{https://github.com/GeoBosh/StableEstim/blob/master/R/Simulation.R}
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Kim, Y. s., Rachev, S. T., Bianchi, M. L. & Fabozzi, F. J. (2008), 'Financial
#' market models with lévy processes and time-varying volatility'
#' \doi{10.1016/j.jbankfin.2007.11.004}
#'
#' Hansen, L. P. (1982), 'Large sample properties of generalized method of
#' moments estimators' \doi{10.2307/1912775}
#'
#' Hansen, L. P.; Heaton, J. & Yaron, A. (1996), 'Finite-Sample Properties of
#' Some Alternative GMM Estimators' \doi{10.1080/07350015.1996.10524656}
#'
#' Carrasco, M. & Kotchoni, R. (2017), 'Efficient estimation using the
#' characteristic function' \doi{10.1017/S0266466616000025}
#'
#' Kuechler, U. & Tappe, S. (2013), 'Tempered stable distribution and processes'
#' \doi{10.1016/j.spa.2013.06.012}
#'
#' Feuerverger, A. & McDunnough, P. (1981), 'On the efficiency of empirical
#' characteristic function procedures'
#' \doi{10.1111/j.2517-6161.1981.tb01143.x}
#'
#' @param TemperedType A String. Either "CTS", "TSS", "NTS", "MTS", "GTS",
#' "KRTS", "RDTS".
#' @param EstimMethod A String. Either "ML", "GMM", "Cgmm", or "GMC".
#' @param data Data used to perform the estimation: numeric vector of length n.
#' @param theta0 A vector of numeric values corresponding to the pattern of the
#' \code{TemperedType}.
#' @param ComputeCov 	Logical flag: If set to TRUE, the asymptotic covariance
#' matrix is computed. \code{FALSE} by default.
#' @param HandleError Logical flag: If set to \code{TRUE} and if an error occurs
#' during the estimation procedure, the computation will carry on and NA will be
#' returned. Useful for Monte Carlo simulations.\code{TRUE} by default.
#' @param eps Numerical error tolerance. \code{1e-06} by default.
#' @param algo algorithm: For GMM: \code{"2SGMM"} is the two step GMM proposed
#' by Hansen (1982). \code{"CueGMM"} and \code{"ITGMM"} are respectively the
#' continuous updated and the iterative GMM proposed by Hansen, Eaton et Yaron
#' (1996) and adapted to the continuum case. For GMC: \code{"2SGMC", "CueGMC"}.
#' For Cgmm: \code{"2SCgmm", "CueCgmm", ...}.
#' @param regularization regularization scheme to be used for moment methods,
#' one of \code{"Tikhonov"} (Tikhonov), \code{"LF"} (Landweber-Fridmann) and
#' \code{"cut-off"} (spectral cut-off).
#' @param WeightingMatrix type of weighting matrix used to compute the
#' objective function for the GMM and GMC methods, one of \code{"OptAsym"} (the
#' optimal asymptotic), \code{"DataVar"} (the data driven, only for GMM) and
#' \code{"Id"} (the identity matrix).
#' @param t_scheme scheme used to select the points for the GMM method where the
#' moment conditions are evaluated, one of \code{"equally"} (equally placed),
#' \code{"NonOptAr"} (non optimal arithmetic placement), \code{"uniformOpt"}
#' (uniform optimal placement), \code{"ArithOpt"} (arithmetic optimal
#' placement), \code{"Var Opt"} (optimal variance placement) and \code{"free"}
#' (users need to pass their own set of points in ...).
#' @param alphaReg value of the regularisation parameter; numeric. Example Value
#' could be ==0.01.
#' @param t_free sequence, if \code{t_scheme=="free"}.
#' @param subdivisions 	Number of subdivisions used to compute the different
#' integrals involved in the computation of the objective function for the Cgmm
#' method (to minimise); numeric.
#' @param IntegrationMethod Numerical integration method to be used to
#' approximate the (vectorial) integrals for the Cgmm method. Users can choose
#' between "Uniform" discretization or the "Simpson"'s rule (the 3-point
#' Newton-Cotes quadrature rule).
#' @param randomIntegrationLaw Probability measure associated to the Hilbert
#' space spanned by the moment conditions for the Cgmm method.
#' @param s_min,s_max Lower and Upper bounds of the interval where the moment
#' conditions are considered for the Cgmm method; numeric.
#' @param ncond Integer. Number of moment conditions (until order \code{ncond})
#' for the GMC method. Must not be less than 3 for TSS, 6 for CTS, 5 for NTS.
#' @param IterationControl only used if algo = "IT..." or algo = "Cue..."
#' to control the iterations. See Details.
#' @param ... Other arguments to be passed to the estimation function or the
#' asymptotic confidence level.
#'
#'
#' @return Object of a estim-class. See details for more information.
#'
#' @examples
#' \donttest{
#' TemperedEstim(TemperedType = "CTS", EstimMethod = "ML",
#'                data = rCTS(2,1.5,1,1,1,1,0),
#'                theta0 = c(1.5,1,1,1,1,0) - 0.1);
#' TemperedEstim("TSS", "GMM", rTSS(20,0.5,1,1), algo = "2SGMM",
#'               alphaReg = 0.01, regularization = "cut-off",
#'               WeightingMatrix = "OptAsym", t_scheme = "free",
#'               t_free = seq(0.1,2,length.out = 12));
#' TemperedEstim("NTS", "Cgmm", rNTS(20,0.5,1,1,1,0), algo = "2SCgmm",
#'               alphaReg = 0.01, subdivisions = 50,
#'               IntegrationMethod = "Uniform", randomIntegrationLaw = "unif",
#'               s_min = 0, s_max= 1);
#' TemperedEstim("TSS", "GMC", rTSS(20, 0.5, 1, 1), algo = "2SGMC",
#'               alphaReg = 0.01, WeightingMatrix = "OptAsym",
#'               regularization = "cut-off", ncond = 8, theta0 = c(0.5,1,1));
#' }
#' @export
TemperedEstim <- function(TemperedType = c("CTS", "TSS", "NTS", "MTS", "GTS",
                                           "KRTS", "RDTS"),
                          EstimMethod = c("ML", "GMM", "Cgmm", "GMC"), data,
                          theta0 = NULL, ComputeCov = FALSE, HandleError = TRUE,
                          eps = 1e-06, algo = NULL, regularization = NULL,
                          WeightingMatrix = NULL, t_scheme = NULL,
                          alphaReg = NULL, t_free = NULL, nb_t = NULL,
                          subdivisions = NULL, IntegrationMethod = NULL,
                          randomIntegrationLaw = NULL, s_min = NULL,
                          s_max = NULL, ncond = NULL, IterationControl = NULL,
                          ...) {
    if (missing(data))
        stop("data not provided !")
    if (is.null(theta0)) {
        if (TemperedType == "CTS") {
            theta0 <- MoC_CTS(x <- data, c(1.5, 1, 1, 1, 1, 0))
        } else if (TemperedType == "TSS") {
            theta0 <- MoC_TSS(x <- data, c(0.5, 1, 1))
        } else if (TemperedType == "NTS") {
            theta0 <- MoC_NTS(x <- data, c(0.5, 0, 1, 1, 0))
        } else if (TemperedType == "MTS") {
            theta0 <- MoC_MTS(x <- data, c(0.6, 1, 1, 1, 0))
        } else if (TemperedType == "GTS") {
            theta0 <- MoC_GTS(x <- data, c(1.5, 1.5, 1, 1, 1, 1, 0))
        } else if (TemperedType == "KRTS") {
            theta0 <- MoC_KRTS(x <- data, c(1.5, 1, 1, 1, 1, 1, 1, 0))
        } else if (TemperedType == "RDTS") {
          theta0 <- MoC_RDTS(x <- data, c(0.5, 1, 1, 1, 0))
        }

      # else {
      #       theta0 <- MoC_CGMY(x <- data, c(1, 1, 1, 1.5))
      #   }
    }
    if (TemperedType == "CTS") {
        OutputObj <- methods::new(Class="EstimCTSClass",par = numeric(6),
                                  par0 = theta0, vcov = matrix(0, 6, 6),
                                  confint = matrix(0, 6, 2), data = data,
                                  failure = 1)
    } else if (TemperedType == "TSS") {
        OutputObj <- methods::new(Class = "EstimSubClass", par = numeric(3),
                                  par0 = theta0, vcov = matrix(0, 3, 3),
                                  confint = matrix(0, 3, 2), data = data,
                                  failure = 1)
    } else if (TemperedType == "NTS") {
        OutputObj <- methods::new(Class = "EstimNTSClass", par = numeric(5),
                                  par0 = theta0, vcov = matrix(0, 5, 5),
                                  confint = matrix(0, 5, 2), data = data,
                                  failure = 1)
    } else if (TemperedType == "MTS") {
      OutputObj <- methods::new(Class = "EstimMTSClass", par = numeric(5),
                                par0 = theta0, vcov = matrix(0, 5, 5),
                                confint = matrix(0, 5, 2), data = data,
                                failure = 1)
    } else if (TemperedType == "GTS") {
      OutputObj <- methods::new(Class = "EstimGTSClass", par = numeric(7),
                                par0 = theta0, vcov = matrix(0, 7, 7),
                                confint = matrix(0, 7, 2), data = data,
                                failure = 1)
    } else if (TemperedType == "KRTS") {
      OutputObj <- methods::new(Class = "EstimKRTSClass", par = numeric(8),
                                par0 = theta0, vcov = matrix(0, 8, 8),
                                confint = matrix(0, 8, 2), data = data,
                                failure = 1)
    } else if (TemperedType == "RDTS") {
      OutputObj <- methods::new(Class = "EstimRDTSClass", par = numeric(5),
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
    EstimFcts <- getTempEstimFcts(type, method,
                                  eps = eps,
                                  algo = algo,
                                  regularization = regularization,
                                  WeightingMatrix =
                                    WeightingMatrix,
                                  t_scheme = t_scheme,
                                  alphaReg = alphaReg,
                                  t_free = t_free,
                                  nb_t = nb_t,
                                  subdivisions = subdivisions,
                                  IntegrationMethod =
                                    IntegrationMethod,
                                  randomIntegrationLaw =
                                    randomIntegrationLaw,
                                  s_min = s_min,
                                  s_max = s_max,
                                  ncond = ncond,
                                  IterationControl = IterationControl,
                                  x = data,
                                  ...)
    res <- .initResTemp(type, method)
    if (HandleError) {
        tr <- tryCatch(EstimFcts$Params(x = data, theta0 = theta0,
                                        eps = eps,
                                        algo = algo,
                                        regularization = regularization,
                                        WeightingMatrix =
                                          WeightingMatrix,
                                        t_scheme = t_scheme,
                                        alphaReg = alphaReg,
                                        t_free = t_free,
                                        nb_t = nb_t,
                                        subdivisions = subdivisions,
                                        IntegrationMethod =
                                          IntegrationMethod,
                                        randomIntegrationLaw =
                                          randomIntegrationLaw,
                                        s_min = s_min,
                                        s_max = s_max,
                                        ncond = ncond,
                                        IterationControl = IterationControl,
                                        ...),
                       error = function(e) e)
        err <- inherits(tr, "error")
        if (!err) {
            res <- tr
            OutputObj@failure <- 0
        }
    } else {
        res <- EstimFcts$Params(x = data, theta0 = theta0,
                                eps = eps,
                                algo = algo,
                                regularization = regularization,
                                WeightingMatrix =
                                  WeightingMatrix,
                                t_scheme = t_scheme,
                                alphaReg = alphaReg,
                                t_free = t_free,
                                nb_t = nb_t,
                                subdivisions = subdivisions,
                                IntegrationMethod =
                                  IntegrationMethod,
                                randomIntegrationLaw =
                                  randomIntegrationLaw,
                                s_min = s_min,
                                s_max = s_max,
                                ncond = ncond,
                                IterationControl = IterationControl,
                                ...)
        OutputObj@failure <- 0
    }
    OutputObj@par <- NameParamsObjectsTemp(res$Estim$par, type)
    OutputObj@others <- res$Estim
    OutputObj@duration <- as.numeric(res$duration)
    OutputObj@method <- res$method
    if (ComputeCov) {
        OutputObj@vcov <- EstimFcts$CovarianceMat(data = OutputObj@data,
                                                  EstimObj = res,
                                                  eps = eps,
                                                  algo = algo,
                                                  regularization =
                                                    regularization,
                                                  WeightingMatrix =
                                                    WeightingMatrix,
                                                  t_scheme = t_scheme,
                                                  alphaReg = alphaReg,
                                                  t_free = t_free,
                                                  nb_t = nb_t,
                                                  subdivisions = subdivisions,
                                                  IntegrationMethod =
                                                    IntegrationMethod,
                                                  randomIntegrationLaw =
                                                    randomIntegrationLaw,
                                                  s_min = s_min,
                                                  s_max = s_max,
                                                  ncond = ncond,
                                                  IterationControl =
                                                    IterationControl,
                                                  ...)
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
# @param TemperedType A String. Either "CTS", "TSS", "NTS", or
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
# TemperedEstim_v2 <- function(TemperedType = c("CTS", "TSS",
#                                               "NTS", "CGMY"),
#                              EstimMethod = c("ML", "GMM", "Cgmm", "GMC"), data,
#                              theta0 = NULL, ComputeCov = FALSE,
#                              HandleError = TRUE, eps = 1e-06, ...) {
#     if (missing(data))
#         stop("data not provided !")
#     if (is.null(theta0)) {
#         if (TemperedType == "CTS") {
#             theta0 <- MoC_CTS(x = data, c(1.5, 1, 1, 1, 1, 0), eps = eps)
#         } else if (TemperedType == "TSS") {
#             theta0 <- MoC_TSS(x = data, c(0.5, 1, 1), eps = eps)
#         } else if (TemperedType == "NTS") {
#             theta0 <- MoC_NTS(x = data, c(0.5, 0, 1, 1, 0), eps = eps)
#         } else {
#             theta0 <- MoC_CGMY(x = data, c(1, 1, 1, 1.5), eps = eps)
#         }
#     }
#     if (TemperedType == "CTS") {
#         OutputObj <- list(par = numeric(6), par0 = theta0,
#                           vcov = matrix(0, 6, 6), confint = matrix(0, 6, 2),
#                           data = data, failure = 1)
#     } else if (TemperedType == "TSS") {
#         OutputObj <- list(par = numeric(3), par0 = theta0,
#                           vcov = matrix(0, 3, 3), confint = matrix(0, 3, 2),
#                           data = data, failure = 1)
#     } else if (TemperedType == "NTS") {
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
    type = c("CTS", "TSS", "NTS", "MTS", "GTS", "KRTS", "RDTS"),
    method = c("ML", "GMM", "Cgmm", "GMC"),
    eps,
    algo,
    regularization,
    WeightingMatrix,
    t_scheme,
    alphaReg,
    t_free,
    nb_t,
    subdivisions,
    IntegrationMethod,
    randomIntegrationLaw,
    s_min,
    s_max,
    ncond,
    IterationControl,
    x,
    ...){

    if (type == "CTS") {
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
    } else if (type == "TSS") {
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
    } else if (type == "NTS") {
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
    } else if (type == "MTS") {
      Output <- switch(method, Cgmm = {
        list(Params = CgmmParametersEstim_MTS,
             CovarianceMat = .asymptoticVarianceEstimCgmm_MTS,
             methodDes = getCgmmMethodName_MTS)
      }, stop(paste(method, " not taken into account ! For now, only Cgmm works
                    with this TS.")))
      Output
    } else if (type == "GTS") {
      Output <- switch(method, ML = {
        list(Params = MLParametersEstim_GTS,
             CovarianceMat = .asymptoticVarianceEstimML_GTS,
             methodDes = .methodDesML_GTS)
      }, GMM = {
        list(Params = GMMParametersEstim_GTS,
             CovarianceMat = .asymptoticVarianceEstimGMM_GTS,
             methodDes = getGMMmethodName_GTS)
      }, Cgmm = {
        list(Params = CgmmParametersEstim_GTS,
             CovarianceMat = .asymptoticVarianceEstimCgmm_GTS,
             methodDes = getCgmmMethodName_GTS)
      }, GMC = {
        list(Params = GMCParametersEstim_GTS
             , CovarianceMat = .asymptoticVarianceEstimGMC_GTS,
             methodDes = getGMCmethodName_GTS)
      }, stop(paste(method, " not taken into account !")))
      Output
    } else if (type == "KRTS") {
      Output <- switch(method, ML = {
        list(Params = MLParametersEstim_KRTS,
             CovarianceMat = .asymptoticVarianceEstimML_KRTS,
             methodDes = .methodDesML_KRTS)
      }, GMM = {
        list(Params = GMMParametersEstim_KRTS,
             CovarianceMat = .asymptoticVarianceEstimGMM_KRTS,
             methodDes = getGMMmethodName_KRTS)
      }, Cgmm = {
        list(Params = CgmmParametersEstim_KRTS,
             CovarianceMat = .asymptoticVarianceEstimCgmm_KRTS,
             methodDes = getCgmmMethodName_KRTS)
      }, stop(paste(method, " not taken into account !")))
      Output
    } else if (type == "RDTS") {
      Output <- switch(method, ML = {
        list(Params = MLParametersEstim_RDTS,
             CovarianceMat = .asymptoticVarianceEstimML_RDTS,
             methodDes = .methodDesML_RDTS)
      }, GMM = {
        list(Params = GMMParametersEstim_RDTS,
             CovarianceMat = .asymptoticVarianceEstimGMM_RDTS,
             methodDes = getGMMmethodName_RDTS)
      }, Cgmm = {
        list(Params = CgmmParametersEstim_RDTS,
             CovarianceMat = .asymptoticVarianceEstimCgmm_RDTS,
             methodDes = getCgmmMethodName_RDTS)
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
    if (type == "CTS") {
        npar <- 6
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    } else if (type == "TSS") {
        npar <- 3
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    } else if (type == "NTS") {
        npar <- 5
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    } else if (type == "MTS") {
        npar <- 5
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
              method = paste(type, method, "failed", sep = "_"))
    } else if (type == "GTS") {
        npar <- 7
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
              method = paste(type, method, "failed", sep = "_"))
    } else if (type == "´KRTS") {
        npar <- 8
        list(Estim = list(par = rep(NaN, npar)), duration = 0,
             method = paste(type, method, "failed", sep = "_"))
    } else if (type == "RDTS") {
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
NameParamsObjectsTemp <- function(mat, type = c("CTS", "TSS",
                                                "NTS", "MTS", "GTS", "KRTS",
                                                "RDTS")) {
    if (type == "CTS") {
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
    } else if (type == "TSS") {
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
    } else if (type == "NTS") {
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
    } else if (type == "MTS") {
      parNames <- c("alpha", "delta", "lambda +", "lambda -", "mu")
      minMaxCol <- c("min", "max")
      if (length(mat) == 5) {
        names(mat) <- parNames
      } else if (is.matrix(mat) && nrow(mat) == 5) {
        rownames(mat) <- parNames
        if (ncol(mat) == 2)
          colnames(mat) <- minMaxCol else if (ncol(mat) == 5)
            colnames(mat) <- parNames
      }
    } else if (type == "GTS") {
      parNames <- c("alpha +", "alpha -", "delta +", "delta -", "lambda +",
                    "lambda -", "mu")
      minMaxCol <- c("min", "max")
      if (length(mat) == 7) {
        names(mat) <- parNames
      } else if (is.matrix(mat) && nrow(mat) == 7) {
        rownames(mat) <- parNames
        if (ncol(mat) == 2)
          colnames(mat) <- minMaxCol else if (ncol(mat) == 7)
            colnames(mat) <- parNames
      }
    } else if (type == "KRTS") {
      parNames <- c("alpha", "k +", "k -", "r +", "r -", "p +", "p -", "mu")
      minMaxCol <- c("min", "max")
      if (length(mat) == 8) {
        names(mat) <- parNames
      } else if (is.matrix(mat) && nrow(mat) == 8) {
        rownames(mat) <- parNames
        if (ncol(mat) == 2)
          colnames(mat) <- minMaxCol else if (ncol(mat) == 8)
            colnames(mat) <- parNames
      }
    } else if (type == "RDTS") {
      parNames <- c("alpha", "delta", "lambda +", "lambda -", "mu")
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
CheckParametersRange_MTS <- function(theta) {
  alpha <- theta[1]
  delta <- theta[2]
  lambdap <- theta[3]
  lambdam <- theta[4]
  mu <- theta[5]
  checkParams <- list(alpha = checkRange(alpha, -Inf, 1, "alpha"),
                      delta = checkRange(delta, 0, Inf, "delta"),
                      lambdap = checkRange(lambdap, 0, Inf, "lambda+"),
                      lambdam = checkRange(lambdam, 0, Inf, "lambda-"),
                      mu = checkRange(mu, -Inf, Inf, "mu"))
  .printErr <- function(errList) if (!errList$bool)
    stop(errList$msg)
  lapply(checkParams, .printErr)
}

# No Export.
CheckParametersRange_GTS <- function(theta) {
  alphap <- theta[1]
  alpham <- theta[2]
  deltap <- theta[3]
  deltam <- theta[4]
  lambdap <- theta[5]
  lambdam <- theta[6]
  mu <- theta[7]
  checkParams <- list(alphap = checkRange(alpha, 0, 2, "alpha+"),
                      alpham = checkRange(alpha, 0, 2, "alpha-"),
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
CheckParametersRange_KRTS <- function(theta, alpha0, ...) {
  alpha <- theta[1]
  kp <- theta[2]
  km <- theta[3]
  rp <- theta[4]
  rm <- theta[5]
  pp <- theta[6]
  pm <- theta[7]
  mu <- theta[8]

  #TODO: eigentlich sollte der alpha0 kennen. Es klappt aber einfach nicht...
  # Händisch 0.5 eingesetzt
  checkParams <- list(alpha = checkRange(alpha, 0, 2, "alpha"),
                      kp = checkRange(kp, 0, Inf, "k+"),
                      km = checkRange(km, 0, Inf, "k-"),
                      rp = checkRange(rp, 0, Inf, "r+"),
                      rm = checkRange(rm, 0, Inf, "r-"),
                      pp = checkRange(pp, -0.5, Inf, "p+"),
                      pm = checkRange(pm, -0.5, Inf, "p-"),
                      mu = checkRange(mu, -Inf, Inf, "mu"))
  .printErr <- function(errList) if (!errList$bool)
    stop(errList$msg)
  lapply(checkParams, .printErr)
}

# No Export.
CheckParametersRange_RDTS <- function(theta) {
  alpha <- theta[1]
  delta <- theta[2]
  lambdap <- theta[3]
  lambdam <- theta[4]
  mu <- theta[5]
  checkParams <- list(alpha = checkRange(alpha, 0, 2, "alpha"),
                      delta = checkRange(delta, 0, Inf, "delta"),
                      lambdap = checkRange(lambdap, 0, Inf, "lambda+"),
                      lambdam = checkRange(lambdam, 0, Inf, "lambda-"),
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
                               " should be in the interval [", min, max, "]")))
}

##### Asymptotic Confidence Interval#####

# No export.
AsymptoticConfidenceInterval <- function(thetaEst, n_sample, Cov,
                                         qLaw = stats::qnorm, level = 0.95,
                                         type, ...) {
    if (type == "CTS") {
        nr <- 6
    } else if (type == "TSS") {
        nr <- 3
    } else if (type == "NTS") {
        nr <- 5
    } else if (type == "MTS") {
      nr <- 5
    } else if (type == "GTS") {
      nr <- 7
    } else if (type == "KRTS") {
      nr <- 8
    } else if (type == "RDTS") {
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

  parNames <- switch(type,
                     CTS = c("Alpha", "DeltaP", "DeltaM", "LambdaP",
                                 "LambdaM", "mu"),
                     TSS = c("Alpha", "Delta", "Lambda"),
                     NTS = c("Alpha", "Beta", "Delta", "Lambda",
                                "mu"),
                     MTS = c("Alpha", "Delta", "LambdaP", "LambdaM",
                             "mu"),
                     GTS = c("AlphaP", "AlphaM", "DeltaP", "DeltaM",
                             "LambdaP", "LambdaM", "mu"),
                     KRTS = c("Alpha", "kP", "kM", "rP", "rM",
                             "pP", "pM", "mu"),
                     CGMY = c("C", "G", "M", "Y"))


  minMaxCol <- c("min", "max")

  if (length(mat) > 2 && length(mat) < 9) {
    names(mat) <- parNames
  }
  else if (is.matrix(mat) && nrow(mat) > 2 && is.matrix(mat) && nrow(mat) < 9) {
    rownames(mat) <- parNames
    if (ncol(mat) == 2)
      colnames(mat) <- minMaxCol
    else if (ncol(mat) > 2 && ncol(mat) < 9)
      colnames(mat) <- parNames
  }
  mat
}


##### Classes#####

#### Sub Class ####

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


#### CTS Class ####

# No export.
#' @importFrom methods new
EstimCTSClass <- setClass("EstimCTSClass",
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

setMethod("initialize", "EstimCTSClass",
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
            NameParamsObjects(par, "CTS")
            NameParamsObjects(par0, "CTS")
            NameParamsObjects(vcov, "CTS")
            NameParamsObjects(confint, "CTS")
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

setMethod("show", "EstimCTSClass",
          function(object) {
            cat("*** Tempered Estim CTS, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim CTS) ******* \n")
          })


#### NTS Class ####

# No export.
#' @importFrom methods new
EstimNTSClass <- setClass("EstimNTSClass",
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

setMethod("initialize", "EstimNTSClass",
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
            NameParamsObjects(par, "NTS")
            NameParamsObjects(par0, "NTS")
            NameParamsObjects(vcov, "NTS")
            NameParamsObjects(confint, "NTS")
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

setMethod("show", "EstimNTSClass",
          function(object) {
            cat("*** Tempered Estim NTS, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim NTS) ******* \n")
          })


#### MTS Class ####

# No export.
#' @importFrom methods new
EstimMTSClass <- setClass("EstimMTSClass",
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

setMethod("initialize", "EstimMTSClass",
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
            NameParamsObjects(par, "MTS")
            NameParamsObjects(par0, "MTS")
            NameParamsObjects(vcov, "MTS")
            NameParamsObjects(confint, "MTS")
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

setMethod("show", "EstimMTSClass",
          function(object) {
            cat("*** Tempered Estim MTS, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim MTS) ******* \n")
          })


#### GTS Class ####

# No export.
#' @importFrom methods new
EstimGTSClass <- setClass("EstimGTSClass",
                          slots = list(par = "numeric", par0 = "numeric",
                                       vcov = "matrix", confint = "matrix",
                                       data = "numeric",
                                       sampleSize = "numeric",
                                       others = "list", duration = "numeric",
                                       failure = "numeric",
                                       method = "character"),
                          contains = list(), validity = function(object) {
                            par <- object@par
                            if (length(par) == 7)
                              ansp <- TRUE
                            else ansp <- "Parameter of length different of 7"
                            par0 <- object@par0
                            if (length(par0) == 7)
                              ansp0 <- TRUE
                            else ansp0 <- "Initial Parameter of length different of 7"
                            vcov <- object@vcov
                            if (ncol(vcov) == 7 && nrow(vcov) == 7)
                              anscov <- TRUE
                            else anscov <- "covariance matrix of length different of 7x7"
                            confint <- object@confint
                            if (ncol(confint) == 2 && nrow(confint) == 7)
                              ansconfint <- TRUE
                            else ansconfint <-
                              "confidance intervall matrix of length different of 7x2"
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

setMethod("initialize", "EstimGTSClass",
          function(.Object, par, par0, vcov, confint, method, level, others,
                   data, duration, failure, ...) {
            ## handle missing
            if (missing(par))
              par        <- numeric(7)
            if (missing(par0))
              par0       <- numeric(7)
            if (missing(vcov))
              vcov       <- matrix(nrow = 7, ncol = 7)
            if (missing(confint))
              confint    <- matrix(nrow = 7, ncol = 2)
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
            NameParamsObjects(par, "GTS")
            NameParamsObjects(par0, "GTS")
            NameParamsObjects(vcov, "GTS")
            NameParamsObjects(confint, "GTS")
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

setMethod("show", "EstimGTSClass",
          function(object) {
            cat("*** Tempered Estim GTS, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim GTS) ******* \n")
          })


#### KRTS ####

# No export.
#' @importFrom methods new
EstimKRTSClass <- setClass("EstimKRTSClass",
                          slots = list(par = "numeric", par0 = "numeric",
                                       vcov = "matrix", confint = "matrix",
                                       data = "numeric",
                                       sampleSize = "numeric",
                                       others = "list", duration = "numeric",
                                       failure = "numeric",
                                       method = "character"),
                          contains = list(), validity = function(object) {
                            par <- object@par
                            if (length(par) == 8)
                              ansp <- TRUE
                            else ansp <- "Parameter of length different of 8"
                            par0 <- object@par0
                            if (length(par0) == 8)
                              ansp0 <- TRUE
                            else ansp0 <- "Initial Parameter of length different of 8"
                            vcov <- object@vcov
                            if (ncol(vcov) == 8 && nrow(vcov) == 8)
                              anscov <- TRUE
                            else anscov <- "covariance matrix of length different of 8x8"
                            confint <- object@confint
                            if (ncol(confint) == 2 && nrow(confint) == 8)
                              ansconfint <- TRUE
                            else ansconfint <-
                              "confidance intervall matrix of length different of 8x2"
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

setMethod("initialize", "EstimKRTSClass",
          function(.Object, par, par0, vcov, confint, method, level, others,
                   data, duration, failure, ...) {
            ## handle missing
            if (missing(par))
              par        <- numeric(8)
            if (missing(par0))
              par0       <- numeric(8)
            if (missing(vcov))
              vcov       <- matrix(nrow = 8, ncol = 8)
            if (missing(confint))
              confint    <- matrix(nrow = 8, ncol = 2)
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
            NameParamsObjects(par, "KRTS")
            NameParamsObjects(par0, "KRTS")
            NameParamsObjects(vcov, "KRTS")
            NameParamsObjects(confint, "KRTS")
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

setMethod("show", "EstimKRTSClass",
          function(object) {
            cat("*** Tempered Estim KRTS, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim KRTS) ******* \n")
          })



#### RDTS Class ####

# No export.
#' @importFrom methods new
EstimRDTSClass <- setClass("EstimRDTSClass",
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

setMethod("initialize", "EstimRDTSClass",
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
            NameParamsObjects(par, "RDTS")
            NameParamsObjects(par0, "RDTS")
            NameParamsObjects(vcov, "RDTS")
            NameParamsObjects(confint, "RDTS")
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

setMethod("show", "EstimRDTSClass",
          function(object) {
            cat("*** Tempered Estim RDTS, method Show *** \n")
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
            cat("\n ******* End Show (Tempered Estim RDTS) ******* \n")
          })


#### CGMY Class ####

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
