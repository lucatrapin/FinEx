arq.llk <- function(theta, r, p, x){

  n <- length(r)
  q <- stats::quantile(r, p)

  f <- rep(q, n)
  for (i in 2:n){
    f[i] <- theta[1] + theta[2]*f[i-1] + theta[3]*x[i-1]
  }

  obj <- sum((p - 1*(r < f))*(r - f))

  return(obj)

}
