##### ML#####
.asymptoticVarianceEstimML_STS <- function(data, EstimObj, type = "Subordinator", ...) {
    asymptoticVarianceEstimML_STS(thetaEst = EstimObj$Estim$par, n_sample = length(data), type = type, ...)
}

asymptoticVarianceEstimML_STS <- function(thetaEst, n_sample, type = "Subordinator", subdivisions = 100, ...) {
    NameParamsObjectsTemp(invFisherMatrix_STS(as.numeric(thetaEst), subdivisions)/n_sample, type = type)
}


invFisherMatrix_STS <- function(theta, subdivisions = 100) {
    mat <- matrix(NA, 3, 3)
    integrand <- function(x, i, j) {
        invf <- 1/VectorialDensity_STS(theta, x)
        df <- jacVectorialDensity_STS(theta, x)
        y <- invf * df[, i] * df[, j]
    }
    for (i in 1:3) {
        for (j in 1:i) {
            mat[i, j] <- integrate(f = integrand, lower = -Inf, upper = Inf, i = i, j = j, subdivisions = subdivisions)$value
            mat[j, i] <- mat[i, j]
        }
    }
    solve(mat)
}

VectorialDensity_STS <- function(theta, xi) {
    dSTS(xi, theta[1], theta[2], theta[3])
}

jacVectorialDensity_STS <- function(theta, xi) {
    NumDeriv_jacobian_STS(fctToDeriv = VectorialDensity_STS, WhereFctIsEvaluated = theta, xi = xi)
}

NumDeriv_jacobian_STS <- function(fctToDeriv, WhereFctIsEvaluated, ...) {
    jacobian(fctToDeriv, WhereFctIsEvaluated, method = "Richardson", method.args = list(), ...)
}

##### GMM#####
.asymptoticVarianceEstimGMM_STS <- function(data, EstimObj, type = "Subordinator", eps, ...) {
    V <- solve(GMMasymptoticVarianceEstim_STS(theta = EstimObj$Estim$par, t = EstimObj$tEstim, x = data, eps = eps, ...))/length(data)
    NameParamsObjects(V, type = type)
}

##### CGMM#####
.asymptoticVarianceEstimCgmm_STS <- function(data, EstimObj, type = "Subordinator", ...) {
    V <- ComputeCovarianceCgmm_STS(theta = EstimObj$Estim$par, thetaHat = EstimObj$Estim$par, x = data, ...)
    NameParamsObjects(Mod(ComputeCutOffInverse(V))/length(data), type = type)
}

ComputeCovarianceCgmm_STS <- function(theta, Cmat = NULL, x, alphaReg, thetaHat, s_min, s_max, subdivisions = 50, IntegrationMethod = c("Uniform",
    "Simpson"), randomIntegrationLaw = c("norm", "unif"), ...) {
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    CovMat <- ComputeCgmmFcts_STS(Fct = "Covariance", theta = theta, Cmat = Cmat, x = x, Weighting = "optimal", alphaReg = alphaReg,
        thetaHat = thetaHat, s_min = s_min, s_max = s_max, subdivisions = subdivisions, IntegrationMethod = IntegrationMethod,
        randomIntegrationLaw = randomIntegrationLaw, ...)
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


##### GMC#####
.asymptoticVarianceEstimGMC_STS <- function(data, EstimObj, type = "Subordinator", eps, ...) {
    V <- solve(GMCasymptoticVarianceEstim_STS(theta = EstimObj$Estim$par, ncond = EstimObj$ncond, x = data, eps = eps,
        ...))/length(data)
    NameParamsObjects(V, type = type)
}

GMCasymptoticVarianceEstim_STS <- function(..., theta, x, ncond, WeightingMatrix, alphaReg = 0.01, regularization = "Tikhonov",
    eps) {
    K <- ComputeGMCWeightingMatrix_STS(theta = theta, x = x, ncond = ncond, WeightingMatrix = WeightingMatrix, ...)
    B <- jacobianSampleRealCFMoment_STS(t, theta)
    fct <- function(G) ComputeInvKbyG_STS(K = K, G = G, alphaReg = alphaReg, regularization = regularization, eps = eps)
    invKcrossB <- apply(X = B, MARGIN = 2, FUN = fct)
    crossprod(B, invKcrossB)
}
