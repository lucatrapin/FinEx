rpot.output <- function(...){

  arg <- list(...)

  # Link functions
  fp.link <- link.fun(arg$lf[1])
  fs.link <- link.fun(arg$lf[2])
  fx.link <- link.fun(arg$lf[3])

  # Number of covariates
  n.xp <- length(arg$xp)
  n.xs <- length(arg$xs)
  n.xx <- length(arg$xx)

  if (arg$mod=="s"){
    phi   <- fp.link(as.matrix(arg$x[,arg$xp])%*%arg$fit$par[1:n.xp])
    sigma <- fs.link(as.matrix(arg$x[,arg$xs])%*%arg$fit$par[(n.xp+1):(n.xp+n.xs)])
    xi    <- fx.link(as.matrix(arg$x[,arg$xx])%*%arg$fit$par[(n.xp+n.xs+1):(n.xp+n.xs+n.xx)])
  } else if (arg$mod=="ar"){
    n <- length(arg$y)
    if (n.xp>1){
      phi <- rep(arg$fit$par[1]/(1-arg$fit$par[n.xp+1]), n)
      for (i in 2:n) phi[i] <- arg$x[i,arg$xp]%*%arg$fit$par[1:n.xp] + phi[i-1]*exp(arg$fit$par[n.xp+1])/(1+exp(arg$fit$par[n.xp+1]))
    } else{
      phi <- rep(arg$fit$par[n.xp], n)
    }
    it <- n.xp + ifelse(n.xp>1, 1, 0)
    if (n.xs>1){
      sigma <- rep(arg$fit$par[it+1]/(1-arg$fit$par[it+n.xs+1]), n)
      for (i in 2:n) sigma[i] <- arg$x[i,arg$xs]%*%arg$fit$par[(it+1):(it+n.xs)] + sigma[i-1]*exp(arg$fit$par[(it+n.xs+1)])/(1+exp(arg$fit$par[(it+n.xs+1)]))
    }else{
      sigma <- rep(arg$fit$par[(it+1)], n)
    }
    it <- n.xp + ifelse(n.xp>1, 1, 0) + n.xs +  ifelse(n.xs>1, 1, 0)
    if (n.xx>1){
      xi <- rep(arg$fit$par[it+1]/(1-arg$fit$par[it+n.xx+1]), n)
      for (i in 2:n) xi[i] <- arg$x[i,arg$xx]%*%arg$fit$par[(it+1):(it+n.xx)] + xi[i-1]*exp(arg$fit$par[(it+n.xx+1)])/(1+exp(arg$fit$par[(it+n.xx+1)]))
    }else{
      xi <- rep(arg$fit$par[(it+1)], n)
    }
    phi <- fp.link(phi)
    sigma <- fs.link(sigma)
    xi <- fx.link(xi)
  }

  # Output
  out <- list()

  out$model <- arg$mod

  out$estimates <- arg$fit$par
  if (arg$mod=="s"){
    names(out$estimates) <- c(sprintf('psi%d',1:n.xp), sprintf('gamma%d',1:n.xs), sprintf('delta%d',1:n.xx))
  }else if (arg$mod=="ar"){
    names(out$estimates) <- c(sprintf('psi%d',1:n.xp),
                              if(n.xp>1) 'psi.ar' else NULL,
                              sprintf('gamma%d',1:n.xs),
                              if(n.xs>1) 'gamma.ar' else NULL,
                              sprintf('delta%d',1:n.xx),
                              if(n.xx>1) 'delta.ar' else NULL)
  }

  # Compute standard errors
  # V <- gpd.vcov(arg$fit$par, arg$y, arg$x, arg$xs, arg$xx, arg$lf, arg$na.ix)
  # H <- numDeriv::jacobian(gpd.gradient, arg$fit$par, , , , arg$y, arg$x, arg$xs, arg$xx, arg$lf, arg$na.ix)
  # out$se <- tryCatch(sqrt(diag(solve(H)%*%V%*%solve(H))), error=function(e){return(rep(NaN, n.xs + n.xx))})
  # names(out$se) <- c(sprintf('gamma%d',1:n.xs),
  #                    sprintf('delta%d',1:n.xx))

  out$value <- - arg$fit$value
  names(out$value) <- 'max.llk'

  out$parameters <- cbind(phi, sigma, xi)
  colnames(out$parameters) <- c('phi', 'sigma', 'xi')

  out$lf <- arg$lf

  out$threshold <- arg$q

  return(out)

}
