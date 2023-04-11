##### ML#####
.asymptoticVarianceEstimML_TSS <- function(data, EstimObj,
                                           type = "TSS",
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
    asymptoticVarianceEstimML_TSS(thetaEst = EstimObj$Estim$par,
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

asymptoticVarianceEstimML_TSS <- function(thetaEst, n_sample,
                                          type = "TSS",
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
    NameParamsObjectsTemp(invFisherMatrix_TSS(as.numeric(thetaEst),
                                              subdivisions)/n_sample,
                          type = type)
}


invFisherMatrix_TSS <- function(theta, subdivisions = 100) {
    mat <- matrix(NA, 3, 3)
    integrand <- function(x, i, j) {
        invf <- 1/VectorialDensity_TSS(theta, x)
        df <- jacVectorialDensity_TSS(theta, x)
        y <- invf * df[, i] * df[, j]
    }
    for (i in 1:3) {
        for (j in 1:i) {
            mat[i, j] <- stats::integrate(f = integrand, lower = -Inf,
                                          upper = Inf, i = i, j = j,
                                          subdivisions = subdivisions)$value
            mat[j, i] <- mat[i, j]
        }
    }
    solve(mat)
}

VectorialDensity_TSS <- function(theta, xi) {
    dTSS(xi, theta[1], theta[2], theta[3])
}

jacVectorialDensity_TSS <- function(theta, xi) {
    NumDeriv_jacobian_TSS(fctToDeriv = VectorialDensity_TSS,
                          WhereFctIsEvaluated = theta, xi = xi)
}

NumDeriv_jacobian_TSS <- function(fctToDeriv, WhereFctIsEvaluated, ...) {
    numDeriv::jacobian(fctToDeriv, WhereFctIsEvaluated, method = "Richardson",
             method.args = list(), ...)
}

##### GMM#####
.asymptoticVarianceEstimGMM_TSS <- function(data, EstimObj,
                                            type = "TSS", eps,
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
    V <- solve(GMMasymptoticVarianceEstim_TSS(theta = EstimObj$Estim$par,
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
.asymptoticVarianceEstimCgmm_TSS <- function(data, EstimObj,
                                             type = "TSS",
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
    V <- ComputeCovarianceCgmm_TSS(theta = EstimObj$Estim$par,
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

ComputeCovarianceCgmm_TSS <- function(theta, Cmat = NULL, x, alphaReg,
                                      thetaHat, s_min, s_max, subdivisions = 50,
                                      IntegrationMethod = c("Uniform",
                                                            "Simpson"),
                                      randomIntegrationLaw = c("norm", "unif"),
                                      ...) {
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    CovMat <- ComputeCgmmFcts_TSS(Fct = "Covariance", theta = theta,
                                  Cmat = Cmat, x = x, Weighting = "optimal",
                                  alphaReg = alphaReg, thetaHat = thetaHat,
                                  s_min = s_min, s_max = s_max,
                                  subdivisions = subdivisions,
                                  IntegrationMethod = IntegrationMethod,
                                  randomIntegrationLaw = randomIntegrationLaw,
                                  ...)
    CovMat/(n - 3)
}

ComputeCutOffInverse <- function(X, alphaReg = 0.001) {
    s <- getSingularValueDecomposition(X)
    index <- (abs(s$lambda) < alphaReg)
    if (any(index)) {
        lambda <- s$lambda
        lambda[index] <- alphaReg
        D <- diag(lambda)
        Invmat <- s$ksi %*% D %*% t(s$phi)
    } else {
        qx <- qr(X)
        Invmat <- solve.qr(qx)
    }
    Invmat
}

# Added by Cedric 20220811
## if Kn is a symmetric matrix, we use eigenvalues analysis
getSingularValueDecomposition <- function(Kn){
  if (isSymmetric(Kn)){
    SingularValuesDecomposition <- eigen(x=Kn, symmetric=TRUE)
    phi <- SingularValuesDecomposition$vectors
    ksi <- SingularValuesDecomposition$vectors
    lambda <- SingularValuesDecomposition$values
  }
  else {
    SingularValuesDecomposition <- svd(Kn)
    phi <- SingularValuesDecomposition$v
    ksi <- SingularValuesDecomposition$u
    lambda <- SingularValuesDecomposition$d
  }
  return(list(lambda=lambda,phi=phi,ksi=ksi))
}


##### GMC#####
.asymptoticVarianceEstimGMC_TSS <- function(data, EstimObj,
                                            type = "TSS", eps,
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
    V <- solve(GMCasymptoticVarianceEstim_TSS(theta = EstimObj$Estim$par,
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

GMCasymptoticVarianceEstim_TSS <- function(..., theta, x, ncond,
                                           WeightingMatrix, alphaReg = 0.01,
                                           regularization = "Tikhonov",
    eps) {
    K <- ComputeGMCWeightingMatrix_TSS(theta = theta, x = x, ncond = ncond,
                                       WeightingMatrix = WeightingMatrix, ...)
    B <- jacobianSampleRealCFMoment_TSS(t, theta)
    fct <- function(G) ComputeInvKbyG_TSS(K = K, G = G, alphaReg = alphaReg,
                                          regularization = regularization,
                                          eps = eps)
    invKcrossB <- apply(X = B, MARGIN = 2, FUN = fct)
    crossprod(B, invKcrossB)
}
