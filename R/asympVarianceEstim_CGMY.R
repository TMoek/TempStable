##### ML#####
.asymptoticVarianceEstimML_CGMY <- function(data, EstimObj,
                                            type = "Subordinator", ...) {
    asymptoticVarianceEstimML_CGMY(thetaEst = EstimObj$Estim$par,
                                   n_sample = length(data), type = type, ...)
}

asymptoticVarianceEstimML_CGMY <- function(thetaEst, n_sample,
                                           type = "Subordinator", subdivisions = 100, ...) {
    NameParamsObjectsTemp(invFisherMatrix_CGMY(as.numeric(thetaEst),
                                               subdivisions)/n_sample,
                          type = type)
}


invFisherMatrix_CGMY <- function(theta, subdivisions = 100) {
    mat <- matrix(NA, 4, 4)
    integrand <- function(x, i, j) {
        invf <- 1/VectorialDensity_CGMY(theta, x)
        df <- jacVectorialDensity_CGMY(theta, x)
        y <- invf * df[, i] * df[, j]
    }
    for (i in 1:4) {
        for (j in 1:i) {
            mat[i, j] <- integrate(f = integrand, lower = -Inf, upper = Inf,
                                   i = i, j = j,
                                   subdivisions = subdivisions)$value
            mat[j, i] <- mat[i, j]
        }
    }
    solve(mat)
}

VectorialDensity_CGMY <- function(theta, xi) {
    dCGMY(xi, theta[1], theta[2], theta[3], theta[4])
}

jacVectorialDensity_CGMY <- function(theta, xi) {
    NumDeriv_jacobian_CGMY(fctToDeriv = VectorialDensity_CGMY,
                           WhereFctIsEvaluated = theta, xi = xi)
}

NumDeriv_jacobian_CGMY <- function(fctToDeriv, WhereFctIsEvaluated, ...) {
    numDeriv::jacobian(fctToDeriv, WhereFctIsEvaluated, method = "Richardson",
             method.args = list(), ...)
}

##### GMM#####
.asymptoticVarianceEstimGMM_CGMY <- function(data, EstimObj,
                                             type = "Subordinator", eps, ...) {
    V <- solve(GMMasymptoticVarianceEstim_CGMY(theta = EstimObj$Estim$par,
                                               t = EstimObj$tEstim,
                                               x = data, eps = eps, ...))/
      length(data)
    NameParamsObjects(V, type = type)
}

##### CGMM#####
.asymptoticVarianceEstimCgmm_CGMY <- function(data, EstimObj,
                                              type = "Subordinator", ...) {
    V <- ComputeCovarianceCgmm_CGMY(theta = EstimObj$Estim$par,
                                    thetaHat = EstimObj$Estim$par,
                                    x = data, ...)
    NameParamsObjects(Mod(ComputeCutOffInverse(V))/length(data), type = type)
}

ComputeCovarianceCgmm_CGMY <- function(theta, Cmat = NULL, x, alphaReg,
                                       thetaHat, s_min, s_max,
                                       subdivisions = 50,
                                       IntegrationMethod = c("Uniform",
                                                             "Simpson"),
                                       randomIntegrationLaw = c("norm", "unif"),
                                       ...) {
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    CovMat <- ComputeCgmmFcts_CGMY(Fct = "Covariance", theta = theta,
                                   Cmat = Cmat, x = x, Weighting = "optimal",
                                   alphaReg = alphaReg, thetaHat = thetaHat,
                                   s_min = s_min, s_max = s_max,
                                   subdivisions = subdivisions,
                                   IntegrationMethod = IntegrationMethod,
                                   randomIntegrationLaw = randomIntegrationLaw,
                                   ...)
    CovMat/(n - 4)
}

##### GMC#####
.asymptoticVarianceEstimGMC_CGMY <- function(data, EstimObj,
                                             type = "Subordinator", eps, ...) {
    V <- solve(GMCasymptoticVarianceEstim_CGMY(theta = EstimObj$Estim$par,
                                               ncond = EstimObj$ncond, x = data,
                                               eps = eps, ...))/length(data)
    NameParamsObjects(V, type = type)
}

GMCasymptoticVarianceEstim_CGMY <- function(..., theta, x, ncond,
                                            WeightingMatrix, alphaReg = 0.01,
                                            regularization = "Tikhonov", eps) {
    K <- ComputeGMCWeightingMatrix_CGMY(theta = theta, x = x, ncond = ncond,
                                        WeightingMatrix = WeightingMatrix, ...)
    B <- jacobianSampleRealCFMoment_CGMY(t, theta)
    fct <- function(G) ComputeInvKbyG_CGMY(K = K, G = G, alphaReg = alphaReg,
                                           regularization = regularization,
                                           eps = eps)
    invKcrossB <- apply(X = B, MARGIN = 2, FUN = fct)
    crossprod(B, invKcrossB)
}
