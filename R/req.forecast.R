#' Realized extreme quantile forecasts
#'
#' Computes one-day-ahead forecasats of Value at Risk (VaR) and Expected Shortfall (ES) for a fitted REQ model.
#'
#' @param fit Output object of the \code{arq()} function.
#' @param x Scalar with the realized measures used in the forecast.
#' @param alpha Probability level for the computation of the risk measures.
#'
#' @return A matrix containing the VaR and ES forecasts at level alpha.
#'
#' @references
#' Bee, M., Dupuis, D. J., and Trapin, L. (2018). Realized extreme quantile: A joint model for conditional quantiles and measures of volatility with EVT refinements. \emph{Journal of Applied Econometrics}, 33(3), 398-415.
#'
#' @author Luca Trapin
#'
#' @export req.forecast
req.forecast <- function(fit, x, alpha){

  q <- fit$estimates[1] + fit$estimates[2]*utils::tail(fit$q,1) + fit$estimates[3]*x
  out <- q*fit$zk*(fit$p/alpha)^fit$estimates[4]

  return(out)

}

