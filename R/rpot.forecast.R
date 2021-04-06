#' Realized Peaks over Threshold forecasts
#'
#' Computes one-day-ahead forecasts of Value at Risk (VaR) and Expected Shortfall (ES) for a fitted RPoT model.
#'
#' @param fit Output object of the \code{rpot()} function.
#' @param x Vector of realized measures used in the forecast. The intercept must be included.
#' @param alpha Probability level for the computation of the risk measures.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{parameters}}{Vector containing the forecasts of phi, sigma, and xi.}
#' \item{\code{var}}{Var forecast at level alpha.}
#' \item{\code{es}}{ES forecast at level alpha.}
#' }
#'
#' @details The risk measures are computed as in Bee et al. (2019):
#' \deqn{VaR^{\alpha}=q+\frac{\sigma}{\xi}\left(\frac{\phi}{1-\alpha}^{\xi}-1\right)}
#' \deqn{ES^{\alpha}=\frac{VaR^{\alpha}}{1-\xi}+\frac{\sigma-q\xi}{1-\xi}}
#'
#' @references
#' Bee, M., Dupuis, D. J., and Trapin, L. (2019). Realized Peaks over Threshold: A Time-Varying Extreme Value Approach with High-Frequency-Based Measures. \emph{Journal of Financial Econometrics}, 17(2), 254-283.
#'
#' @author Luca Trapin
#'
#' @export rpot.forecast
rpot.forecast <- function(fit, x, alpha){

  tau   <- fit$threshold
  phi   <- fit$parameters[,1]
  sigma <- fit$parameters[,2]
  xi    <- fit$parameters[,3]

  var <- tau + sigma/xi*((phi/(1-alpha))^xi-1)
  es  <- var/(1-xi) + (sigma-xi*tau)/(1-xi)

  if (fit$model=="ar"){
    n <- nrow(fit$parameters)
    fp <- link.fun(fit$lf[1])
    fs <- link.fun(fit$lf[2])
    fx <- link.fun(fit$lf[3])
    fp.inv <- link.fun.inv(fit$lf[1])
    fs.inv <- link.fun.inv(fit$lf[2])
    fx.inv <- link.fun.inv(fit$lf[3])
    psi <- fp.inv(fit$parameters[n,1])
    gamma <- fs.inv(fit$parameters[n,2])
    delta <- fx.inv(fit$parameters[n,3])

    psi   <- if (sum(substr(names(fit$estimates),1,3)=="psi")>1) fit$estimates[substr(names(fit$estimates),1,3)=="psi"]%*%c(x, psi) else fit$estimates[substr(names(fit$estimates),1,3)=="psi"]
    gamma <- if (sum(substr(names(fit$estimates),1,5)=="gamma")>1) fit$estimates[substr(names(fit$estimates),1,5)=="gamma"]%*%c(x, gamma) else fit$estimates[substr(names(fit$estimates),1,5)=="gamma"]
    delta <- if (sum(substr(names(fit$estimates),1,5)=="delta")>1) fit$estimates[substr(names(fit$estimates),1,5)=="delta"]%*%c(x, delta) else fit$estimates[substr(names(fit$estimates),1,5)=="delta"]

    parameters <- c(fp(psi), fs(gamma), fx(delta))
    names(parameters) <- c("phi", "sigma", "xi")

  } else {

    psi   <- if (sum(substr(names(fit$estimates),1,3)=="psi")>1) fit$estimates[substr(names(fit$estimates),1,3)=="psi"]%*% x else fit$estimates[substr(names(fit$estimates),1,3)=="psi"]
    gamma <- if (sum(substr(names(fit$estimates),1,5)=="gamma")>1) fit$estimates[substr(names(fit$estimates),1,5)=="gamma"]%*% x else fit$estimates[substr(names(fit$estimates),1,5)=="gamma"]
    delta <- if (sum(substr(names(fit$estimates),1,5)=="delta")>1) fit$estimates[substr(names(fit$estimates),1,5)=="delta"]%*% x else fit$estimates[substr(names(fit$estimates),1,5)=="delta"]

    parameters <- c(fp(psi), fs(gamma), fx(delta))
    names(parameters) <- c("phi", "sigma", "xi")

  }

  out <- list()
  out$parameters <- parameters
  out$var <- fit$threshold + parameters[2]/parameters[3]*((parameters[1]/(1-alpha))^parameters[3]-1)
  out$es  <- out$var/(1-parameters[3]) + (parameters[2]-parameters[3]*fit$threshold)/(1-parameters[3])

  return(out)

}
