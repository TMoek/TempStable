##### ML#####
.asymptoticVarianceEstimML_CTS <- function(data, EstimObj,
                                           type = "CTS",
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
    asymptoticVarianceEstimML_CTS(thetaEst = EstimObj$Estim$par,
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

asymptoticVarianceEstimML_CTS <- function(thetaEst, n_sample,
                                          type = "CTS",
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
    NameParamsObjectsTemp(invFisherMatrix_CTS(as.numeric(thetaEst),
                                              subdivisions)/n_sample,
                          type = type)
}


invFisherMatrix_CTS <- function(theta, subdivisions = 100) {
    mat <- matrix(NA, 6, 6)
    integrand <- function(x, i, j) {
        invf <- 1/VectorialDensity_CTS(theta, x)
        df <- jacVectorialDensity_CTS(theta, x)
        y <- invf * df[, i] * df[, j]
    }
    for (i in 1:6) {
        for (j in 1:i) {
            mat[i, j] <- stats::integrate(f = integrand, lower = -Inf, upper = Inf,
                                   i = i, j = j,
                                   subdivisions = subdivisions)$value
            mat[j, i] <- mat[i, j]
        }
    }
    solve(mat)
}

VectorialDensity_CTS <- function(theta, xi) {
    dCTS(xi, theta[1], theta[2], theta[3], theta[4], theta[5], theta[6])
}

jacVectorialDensity_CTS <- function(theta, xi) {
    NumDeriv_jacobian_CTS(fctToDeriv = VectorialDensity_CTS,
                          WhereFctIsEvaluated = theta, xi = xi)
}

NumDeriv_jacobian_CTS <- function(fctToDeriv, WhereFctIsEvaluated, ...) {
    numDeriv::jacobian(fctToDeriv, WhereFctIsEvaluated, method = "Richardson",
             method.args = list(), ...)
}

##### GMM#####
.asymptoticVarianceEstimGMM_CTS <- function(data, EstimObj,
                                            type = "CTS", eps,
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
    V <- solve(GMMasymptoticVarianceEstim_CTS(theta = EstimObj$Estim$par,
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
.asymptoticVarianceEstimCgmm_CTS <- function(data, EstimObj,
                                             type = "CTS",
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
    V <- ComputeCovarianceCgmm_CTS(theta = EstimObj$Estim$par,
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

ComputeCovarianceCgmm_CTS <- function(theta, Cmat = NULL, x, alphaReg,
                                      thetaHat, s_min, s_max, subdivisions = 50,
                                      IntegrationMethod = c("Uniform",
                                                            "Simpson"),
                                      randomIntegrationLaw = c("norm",
                                                               "unif"), ...) {
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    CovMat <- ComputeCgmmFcts_CTS(Fct = "Covariance", theta = theta,
                                  Cmat = Cmat, x = x, Weighting = "optimal",
                                  alphaReg = alphaReg, thetaHat = thetaHat,
                                  s_min = s_min, s_max = s_max,
                                  subdivisions = subdivisions,
                                  IntegrationMethod = IntegrationMethod,
                                  randomIntegrationLaw = randomIntegrationLaw,
                                  ...)
    CovMat/(n - 6)
}

##### GMC#####
.asymptoticVarianceEstimGMC_CTS <- function(data, EstimObj,
                                            type = "CTS", eps,
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
    V <- solve(GMCasymptoticVarianceEstim_CTS(theta = EstimObj$Estim$par,
                                              ncond = EstimObj$ncond, x = data,
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
                                              IterationControl =
                                                IterationControl,
                                              ...))/length(data)
    NameParamsObjects(V, type = type)
}

GMCasymptoticVarianceEstim_CTS <- function(..., theta, x, ncond,
                                           WeightingMatrix, alphaReg = 0.01,
                                           regularization = "Tikhonov", eps) {
    K <- ComputeGMCWeightingMatrix_CTS(theta = theta, x = x, ncond = ncond,
                                       WeightingMatrix = WeightingMatrix, ...)
    B <- jacobianSampleRealCFMoment_CTS(t, theta)
    fct <- function(G) ComputeInvKbyG_CTS(K = K, G = G, alphaReg = alphaReg,
                                          regularization = regularization,
                                          eps = eps)
    invKcrossB <- apply(X = B, MARGIN = 2, FUN = fct)
    crossprod(B, invKcrossB)
}
