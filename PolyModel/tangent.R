tangent <- function(fit, x){
  hermcoef <- hermite.he.polynomials(4)
  Dherm <- unlist(polynomial.values(polynomial.derivatives(hermcoef[-1]), x))
  slope <- t(Dherm) %*% fit$pbeta[-1]
  intercept <- c(hermite(x, 0), hermite(x, 1), hermite(x, 2), hermite(x, 3), hermite(x, 4)) %*% fit$pbeta- slope *x
  result <- list(slope = slope, intercept = as.numeric(intercept))
  return(result)
}