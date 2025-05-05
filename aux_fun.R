library(INLA)

book.rMatern <- function(n, coords, sigma=1, range, 
                         kappa = sqrt(8*nu)/range, variance = sigma^2, nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}

varF <- function(m) {
  inla.tmarginal(function(x) 1 / sqrt(exp(x)), m)
}