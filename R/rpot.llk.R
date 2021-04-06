rpot.llk <- function(theta, y, x, mod, xp, xs, xx, lf){

  # Link functions
  fp.link <- link.fun(lf[1])
  fs.link <- link.fun(lf[2])
  fx.link <- link.fun(lf[3])

  # Number of covariates
  n.xs <- length(xs)
  n.xx <- length(xx)
  n.xp <- length(xp)

  # Compute parameters
  if (mod == "s"){
    phi   <- fp.link(as.matrix(x[, xp])%*%theta[1:n.xp])
    sigma <- fs.link(as.matrix(x[, xs])%*%theta[(n.xp+1):(n.xp+n.xs)])
    xi    <- fx.link(as.matrix(x[, xx])%*%theta[(n.xp+n.xs+1):(n.xp+n.xs+n.xx)])
  } else if (mod == "ar"){
    n <- length(y)
    if (n.xp>1){
      phi <- rep(theta[1]/(1-theta[n.xp+1]), n)
      for (i in 2:n) phi[i] <- x[i,xp]%*%theta[1:n.xp] + phi[i-1]*exp(theta[n.xp+1])/(1+exp(theta[n.xp+1]))
    } else{
      phi <- rep(theta[n.xp], n)
    }
    it <- n.xp + ifelse(n.xp>1, 1, 0)
    if (n.xs>1){
      sigma <- rep(theta[it+1]/(1-theta[it+n.xs+1]), n)
      for (i in 2:n) sigma[i] <- x[i,xs]%*%theta[(it+1):(it+n.xs)] + sigma[i-1]*exp(theta[(it+n.xs+1)])/(1+exp(theta[(it+n.xs+1)]))
    }else{
      sigma <- rep(theta[(it+1)], n)
    }
    it <- n.xp + ifelse(n.xp>1, 1, 0) + n.xs +  ifelse(n.xs>1, 1, 0)
    if (n.xx>1){
      xi <- rep(theta[it+1]/(1-theta[it+n.xx+1]), n)
      for (i in 2:n) xi[i] <- x[i,xx]%*%theta[(it+1):(it+n.xx)] + xi[i-1]*exp(theta[(it+n.xx+1)])/(1+exp(theta[(it+n.xx+1)]))
    }else{
      xi <- rep(theta[(it+1)], n)
    }
    phi <- fp.link(phi)
    sigma <- fs.link(sigma)
    xi <- fx.link(xi)
  }


  ix <- 1*(y>0)
  llk <- -sum((1-ix)*log(1-phi) + ix*log(phi/sigma) + ix*(-1/xi-1)*log(1+xi*y/sigma))

  return(llk)

}
