##### ML#####
.asymptoticVarianceEstimML_NTS <- function(data, EstimObj,
                                           type = "NTS",
                                           eps,
                                           algo,
                                           regularization,
                                           WeightingMatrix,
                                           t_scheme,
                                           alphaReg,
                                           t_free,
                                           subdivisions,
                                           IntegrationMethod,
                                           randomIntegrationLaw,
                                           s_min,
                                           s_max,
                                           ncond,
                                           IterationControl,
                                           ...) {
    asymptoticVarianceEstimML_NTS(thetaEst = EstimObj$Estim$par,
                                  n_sample = length(data), type = type,
                                  eps = eps,
                                  algo = algo,
                                  regularization = regularization,
                                  WeightingMatrix =
                                    WeightingMatrix,
                                  t_scheme = t_scheme,
                                  alphaReg = alphaReg,
                                  t_free = t_free,
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
}

asymptoticVarianceEstimML_NTS <- function(thetaEst, n_sample,
                                          type = "NTS",
                                          subdivisions = 100,
                                          eps,
                                          algo,
                                          regularization,
                                          WeightingMatrix,
                                          t_scheme,
                                          alphaReg,
                                          t_free,
                                          IntegrationMethod,
                                          randomIntegrationLaw,
                                          s_min,
                                          s_max,
                                          ncond,
                                          IterationControl,
                                          ...) {
    NameParamsObjectsTemp(invFisherMatrix_NTS(as.numeric(thetaEst),
                                              subdivisions)/n_sample,
                          type = type)
}


invFisherMatrix_NTS <- function(theta, subdivisions = 100) {
    mat <- matrix(NA, 5, 5)
    integrand <- function(x, i, j) {
        invf <- 1/VectorialDensity_NTS(theta, x)
        df <- jacVectorialDensity_NTS(theta, x)
        y <- invf * df[, i] * df[, j]
    }
    for (i in 1:5) {
        for (j in 1:i) {
            mat[i, j] <- stats::integrate(f = integrand, lower = -Inf, upper = Inf,
                                   i = i, j = j,
                                   subdivisions = subdivisions)$value
            mat[j, i] <- mat[i, j]
        }
    }
    solve(mat)
}

VectorialDensity_NTS <- function(theta, xi) {
    dNTS(xi, theta[1], theta[2], theta[3], theta[4], theta[5])
}

jacVectorialDensity_NTS <- function(theta, xi) {
    NumDeriv_jacobian_NTS(fctToDeriv = VectorialDensity_NTS,
                          WhereFctIsEvaluated = theta, xi = xi)
}

NumDeriv_jacobian_NTS <- function(fctToDeriv, WhereFctIsEvaluated, ...) {
    numDeriv::jacobian(fctToDeriv, WhereFctIsEvaluated, method = "Richardson",
             method.args = list(), ...)
}

##### GMM#####
.asymptoticVarianceEstimGMM_NTS <- function(data, EstimObj,
                                            type = "NTS", eps,
                                            algo,
                                            regularization,
                                            WeightingMatrix,
                                            t_scheme,
                                            alphaReg,
                                            t_free,
                                            subdivisions,
                                            IntegrationMethod,
                                            randomIntegrationLaw,
                                            s_min,
                                            s_max,
                                            ncond,
                                            IterationControl,
                                            ...) {
    V <- solve(GMMasymptoticVarianceEstim_NTS(theta = EstimObj$Estim$par,
                                              t = EstimObj$tEstim, x = data,
                                              eps = eps,
                                              algo = algo,
                                              regularization = regularization,
                                              WeightingMatrix =
                                                WeightingMatrix,
                                              t_scheme = t_scheme,
                                              alphaReg = alphaReg,
                                              t_free = t_free,
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
                                              ...))/length(data)
    NameParamsObjects(V, type = type)
}

##### CGMM#####
.asymptoticVarianceEstimCgmm_NTS <- function(data, EstimObj,
                                             type = "NTS",
                                             eps,
                                             algo,
                                             regularization,
                                             WeightingMatrix,
                                             t_scheme,
                                             alphaReg,
                                             t_free,
                                             subdivisions,
                                             IntegrationMethod,
                                             randomIntegrationLaw,
                                             s_min,
                                             s_max,
                                             ncond,
                                             IterationControl,
                                             ...) {
    V <- ComputeCovarianceCgmm_NTS(theta = EstimObj$Estim$par,
                                   thetaHat = EstimObj$Estim$par, x = data,
                                   eps = eps,
                                   algo = algo,
                                   regularization = regularization,
                                   WeightingMatrix =
                                     WeightingMatrix,
                                   t_scheme = t_scheme,
                                   alphaReg = alphaReg,
                                   t_free = t_free,
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
    NameParamsObjects(Mod(ComputeCutOffInverse(V))/length(data), type = type)
}

ComputeCovarianceCgmm_NTS <- function(theta, Cmat = NULL, x, alphaReg,
                                      thetaHat, s_min, s_max, subdivisions = 50,
                                      IntegrationMethod = c("Uniform",
                                                            "Simpson"),
                                      randomIntegrationLaw = c("norm",
                                                               "unif"), ...) {
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    CovMat <- ComputeCgmmFcts_NTS(Fct = "Covariance", theta = theta,
                                  Cmat = Cmat, x = x, Weighting = "optimal",
                                  alphaReg = alphaReg, thetaHat = thetaHat,
                                  s_min = s_min, s_max = s_max,
                                  subdivisions = subdivisions,
                                  IntegrationMethod = IntegrationMethod,
                                  randomIntegrationLaw = randomIntegrationLaw,
                                  ...)
    CovMat/(n - 5)
}
