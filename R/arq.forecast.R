#' Dynamic Quantile forecasts
#'
#' Computes one-day-ahead forecasts of the dynamic quantile.
#'
#' @param fit Output object of the \code{arq()} function.
#' @param x Scalar with the realized measures used in the forecast.
#' @param alpha Probability level for the computation of the risk measures.
#'
#' @return The quantile forecast at level alpha.
#'
#' @author Luca Trapin
#'
#' @export arq.forecast
arq.forecast <- function(fit, x, alpha){

  q <- fit$estimates[1] + fit$estimates[2]*utils::tail(fit$q,1) + fit$estimates[3]*x

  return(q)

}
