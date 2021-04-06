rpot.sv <- function(y, x, mod, xp, xs, xx, lf){

  ###########################
  # Unconditional estimates #
  ###########################
  sv    <- ismev::gpd.fit(y, 0, show = F)$mle
  sigma <- sv[1]
  xi    <- sv[2]
  phi   <- mean(y==0)

  if (lf[3]=="exp"){
    xi <- max(xi, 0.001)
  }

  ########################
  # Regression estimates #
  ########################

  # Number of covariates
  n.xp <- length(xp)
  n.xs <- length(xs)
  n.xx <- length(xx)

  # Inverse link functions
  fp.inv <- link.fun.inv(lf[1])
  fs.inv <- link.fun.inv(lf[2])
  fx.inv <- link.fun.inv(lf[3])

  # Compute starting values
  if (mod ==  "s"){
    theta.phi   <- c(fp.inv(phi),   rep(0, n.xp-1))
    theta.sigma <- c(fs.inv(sigma), rep(0, n.xs-1))
    theta.xi    <- c(fx.inv(xi),    rep(0, n.xx-1))
  } else if(mod == "ar"){
    if (n.xp>1)  theta.phi   <- c(fp.inv(phi)*(1-0.9), rep(0, n.xp-1), log(0.9/0.1)) else theta.phi <- fp.inv(phi)
    if (n.xs>1)  theta.sigma <- c(fs.inv(sigma)*(1-0.9), rep(0, n.xs-1), log(0.9/0.1)) else theta.sigma <- fs.inv(sigma)
    if (n.xx>1)  theta.xi    <- c(fx.inv(xi)*(1-0.9), rep(0, n.xx-1), log(0.9/0.1)) else theta.xi <- fx.inv(xi)
  }

  return(c(theta.phi, theta.sigma, theta.xi))

}
