testFunc <- function(a,b){
  .Call("add",a,b, PACKAGE = "TempStable")
}




### fast rejection algorithm, R version

##' Sample a vector of random variates St ~ \tilde{S}(alpha, 1,
##' (cos(alpha*pi/2)*V_0)^{1/alpha}, V_0*I_{alpha = 1},
##' h*I_{alpha != 1}; 1) with LS transform
##' exp(-V_0((h+t)^alpha-h^alpha)) with the fast rejection
##' algorithm; see Nolan's book for the parametrization
##'
##' @title Sampling an exponentially tilted stable distribution
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @return vector of variates St
##' @author Marius Hofert, Martin Maechler
retstableR <- function(alpha, V0, h=1) {
  stopifnot(is.numeric(alpha), length(alpha) == 1,
            0 <= alpha, alpha <= 1) # alpha > 1 => cos(pi/2 *alpha) < 0
  n <- length(V0)
  ## case alpha == 1
  if(alpha == 1 || n == 0) return(V0) # alpha == 1 => point mass at V0
  ## else alpha != 1 => call fast rejection algorithm with optimal m
  m <- m.opt.retst(V0)
  mapply(retstablerej, m=m, V0=V0, alpha=alpha)
}


### state-of-the-art: fast rejection + Luc's algorithm, C version

##' Sample random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)), see Nolan's book for
##' the parametrization, with the fast rejection.
##' This procedure is more efficient than retstableR since it calls the C
##' function retstable_c and uses both the fast rejection and Luc Devroye's algorithm.
##'
##' @title Efficiently sampling an exponentially tilted stable distribution
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @param method which method to call ("Marius Hofert", "Luc Devroye")
##' @return vector of variates St
##' @author Martin Maechler
retstableC <- function(alpha, V0, h = 1, method = NULL) {
  n <- length(V0)
  stopifnot(is.numeric(alpha), length(alpha) == 1,
            0 < alpha, alpha <= 1,
            is.numeric(h), length(h) == 1, h > 0)
  if(alpha == 1 || n == 0) {
    ## alpha == 1 => St corresponds to a point mass at V0 with
    V0           # Laplace-Stieltjes transform exp(-V0*t)
  }
  else {
    if(is.null(method)) {
      if(any(diff(V.is.sml <- V0 * h^alpha < 4))) { ## use *both* methods
        r <- numeric(n)
        r[ V.is.sml] <- .Call(retstable_c, V0[ V.is.sml], h = h, alpha, "MH")
        r[!V.is.sml] <- .Call(retstable_c, V0[!V.is.sml], h = h, alpha, "LD")
        return(r)
      }
      else
        method <- if(V.is.sml[1]) "MH" else "LD"
    }
    else
      method <- match.arg(method, c("MH","LD"))
    .Call(retstable_c, V0, h = h, alpha, method)
  }
}

### fast rejection for fixed m, R version

##' Sample a random variate St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((1+t)^alpha-1)), see Nolan's book for
##' the parametrization, via an m-fold sum of random variates from
##' \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha}, (V_0/m)
##' *I_{alpha = 1}, I_{alpha != 1}; 1) with Laplace-Stieltjes transform
##' exp(-(V_0/m)*((1+t)^alpha-1)). This is a building block for the fast rejection.
##'
##' @title Sample an exponentially tilted stable distribution as an m-fold sum
##' @param m number of summands, any positive integer
##' @param V0 random variate
##' @param alpha parameter in (0,1]
##' @return St
##' @author Marius Hofert, Martin Maechler
retstablerej <- function(m, V0, alpha) {
  gamm. <- (cospi2(alpha)*V0/m)^(1/alpha)
  sum(unlist(lapply(integer(m),
                    function(.) {
                      ## apply standard rejection for sampling
                      ## \tilde{S}(alpha, 1, (cos(alpha*pi/2)
                      ##	*V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1},
                      ## h*I_{alpha != 1}; 1) with Laplace-Stieltjes
                      ## transform exp(-(V_0/m)*((h+t)^alpha-h^alpha))
                      repeat {
                        V__ <- rstable1(1, alpha, beta=1, gamma = gamm.)
                        if(runif(1) <= exp(-V__))
                          return(V__)
                      }})
             ## on acceptance, St_k ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
             ## *V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1}, h*I_{alpha != 1};
             ## 1) with Laplace-Stieltjes transform
             ## exp(-(V_0/m)*((h+t)^alpha-h^alpha))
  ))
}


##' @title tan(pi*x), exact for integer x
##' @param x numeric vector
##' @return numeric vector of values tan(pi*x)
##' @author Martin Maechler
##' tanpi <- function(x) tan(pi * (x %% 1))

##' @title cos(pi/2 * x), exact for integer x
##' @param x numeric vector
##' @return numeric vector of values cos(pi/2 *x)
##' @author Martin Maechler
cospi2 <- function(x) {
  x <- r <- (x %% 4)## cos(pi/2 x) == cos(pi/2(x + 4k))  \forall k \in \Z
  if(any(isI <- x == round(x))) {
    i <- which(isI)
    r[i	      ] <-  0 # for those where x is 1 or 3
    r[i[x[i] == 0]] <-  1
    r[i[x[i] == 2]] <- -1
  }
  io <- which(!isI)
  r[io] <- cospi(x[io]/2)
  r
}


##' Sample S ~ S(alpha, beta, gamma, delta; pm), see package
##' \pkg{stabledist} for the parameterization.
##'
##' @title Sampling stable distributions
##' @param n number of random variates to be generated
##' @param alpha, see code in fBasics
##' @param beta, see code in fBasics
##' @param gamma, see code in fBasics
##' @param delta, see code in fBasics
##' @param pm in {0,1} parameterization, see code in fBasics
##' @return vector of variates S
##' @author Martin Maechler, based on Diethelm Wuertz's code in fBasics
rstable1R <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
{
  stopifnot((la <- length(alpha)) >= 1, length(beta) >= 1,
            length(gamma) >= 1, length(delta) >= 1,
            0 < alpha, alpha <= 2, abs(beta) <= 1,
            length(pm) == 1, pm %in% 0:1)

  p2 <- pi/2
  ## Special case (a,b) = (1,0):
  if (all(alpha == 1) && all(beta == 0)) {
    Z <- rcauchy(n)
  }
  else {
    ## MM: Nolan(2009) "chapt1.pdf", p.21 has  "Theorem 1.19"
    ## -- and attributes that to  'Chambers et al. (1976)'

    ## Calculate uniform and exponential distributed random numbers:
    Theta <- pi * (runif(n)-1/2)
    W <- rexp(n)
    ##  ^^^^^^ was "-log(runif(n))"
    ## rexp() is faster, giving different numbers

    a.is.vec <- (la > 1)          # alpha is "vector" (not scalar)
    if(a.is.vec) {
      ## if alpha is not scalar, make sure that lengths of
      ##		alpha, beta, Theta, W  are all equal (== n)
      alpha <- rep(alpha, length.out = n)
      beta  <- rep(beta,  length.out = n)
    }

    norm <- alpha != 1 ## TODO:  abs(alpha - 1) > eps.alpha1
    ## FIXME(2): ditto for	  | alpha - 1 | << 1

    Z <- numeric(n)
    if(any(norm)) { ## alpha != 1
      alp <- alpha[norm]; Thet <- Theta[norm]
      b.tan.pa <- beta[norm]*tanpi(alp/2)
      th0 <- atan(b.tan.pa) / alp ## == \theta_0
      ## Now, from Nolan/Chambers' formula, we replace
      ##	1/(\cos(\alpha\theta_0) \cos\Theta)^{1/\alpha} with
      ## c / (\cos\Theta)^{1/\alpha} where
      ## c := (1 + (\beta*\tan(\pi\alpha/2))^2)^{1/{2\alpha}}
      ## need to show that c = 1/(\cos(\alpha\theta_0))^{1/\alpha}
      ## <==> 1 + (\beta*\tan(\pi\alpha/2))^2 = 1
      ## / (\cos(\alpha\theta_0))^2 and that's true, as
      ## 1 + (tan(y))^2 = 1 / (cos(y))^2 for
      ##   y = \alpha\theta_0 = \arc\tan(\beta*\tan(\pi\alpha/2))
      c. <- (1 + b.tan.pa^2)^(1/(2*alp))
      a.tht <- alp*(Thet+th0)
      Z[norm] <-
        sin(a.tht) * c. / cos(Thet)^(1/alp) *
        (cos(a.tht-Thet)/W[norm])^((1-alp)/alp)
    }
    ## {note that logicals vectorize, indices do *not* so easily}
    if(any(a1 <- !norm)) { ## alpha == 1
      bet <- beta[a1]; Thet <- Theta[a1]
      p2.bt <- p2 + bet*Thet
      Z[a1] <- (p2.bt*tan(Thet) - bet*log((p2*W[a1]*cos(Thet))/p2.bt))/p2
    }
  }

  if(pm == 0)
    ## delta_1 := delta_0 - .. [Nolan, chapt.1, (1.7)],
    ## since above 1-parametr.
    delta <- delta - gamma * {
      if(a.is.vec) {
        d.of <- numeric(n)
        d.of[norm] <- b.tan.pa
        d.of[a1	 ] <- beta * log(gamma)/p2
        d.of
      }
      else { ## alpha is scalar
        if(norm) b.tan.pa else beta * log(gamma)/p2
      }
    }

  ## Result: Location - Scale trafo -- only now using (gamma, delta):
  Z * gamma + delta
}


##' Sample S ~ S(alpha, beta, gamma, delta; pm), see rstable1R() above.
##' For beta == 1 and pm == 1, the fast C implementation is used.
##'
##' @title Efficiently sampling stable distributions
##' @param n number of random variates to be generated
##' @param alpha, see code in fBasics
##' @param beta, see code in fBasics
##' @param gamma, see code in fBasics
##' @param delta, see code in fBasics
##' @param pm in {0,1} parameterization, see code in fBasics
##' @return vector of variates S
##' @author Martin Maechler
rstable1C <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
{
  stopifnot(length(alpha) >= 1, length(beta) >= 1,
            length(gamma) >= 1, length(delta) >= 1,
            0 < alpha, alpha <= 2, abs(beta) <= 1, gamma >= 0,
            length(pm) == 1, pm %in% 0:1)

  if(beta == 1 && pm == 1)
    .Call(rstable_c, n, alpha) * gamma + delta
  else rstable1R(n, alpha=alpha, beta=beta,
                 gamma=gamma, delta=delta, pm=pm)
}

# je nach Bedarf, kann man das aendern
rstable1 <- rstable1R
# rstable1 <- rstable1C

### Clayton ####################################################################

##' Note: this is mainly to show that this function can be very well
##' approximated much more simply by just using m <- round(V0).
##'
##' @title Optimal constant for fast rejection
##' @param V0 numeric vector >= 0
##' @return optimal constant m for the fast rejection algorithm
##' @author Martin Maechler (based on Marius Hofert's code)
m.opt.retst <- function(V0) {
  n <- length(V0)
  fV <- floor(V0)
  cV <- ceiling(V0)
  v1 <- fV*exp(V0/fV)
  v2 <- cV*exp(V0/cV)

  m <- integer(n)
  l1 <- (V0 <= 1)
  m[which(l1)] <- 1L

  i2 <- which(!l1) ## those with V0 > 1
  l3 <- (v1[i2] <= v2[i2])
  i3 <- i2[l3]
  m[i3] <- fV[i3]
  i4 <- i2[!l3]
  m[i4] <- cV[i4]
  m
}
