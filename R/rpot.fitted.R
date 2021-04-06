#' Realized Peaks over Threshold fitted values
#'
#' Computes fitted Value at Risk (VaR) and Expected Shortfall (ES) for a fitted RPoT model.
#'
#' @param fit Output object of the \code{rpot()} function.
#' @param alpha Probability level for the computation of the risk measures.
#'
#' @return A matrix containing the VaR and ES fitted values at level alpha.
#'
#' @details The risk measures are computed as in Bee et al. (2019):
#' \deqn{VaR^{\alpha}=q+\frac{\sigma}{\xi}\left(\frac{\phi}{1-\alpha}^{\xi}-1\right)}
#' \deqn{ES^{\alpha}=\frac{VaR^{\alpha}}{1-\xi}+\frac{\sigma-q\xi}{1-\xi}}
#'
#' @references
#' Bee, M., Dupuis, D. J., and Trapin, L. (2019). Realized Peaks over Threshold: A Time-Varying Extreme Value Approach with High-Frequency-Based Measures. \emph{Journal of Financial Econometrics}, 17(2), 254-283.
#'
#' @seealso
#' \code{\link{rpot.forecast}}
#'
#' @author Luca Trapin
#'
#' @export rpot.fitted
rpot.fitted <- function(fit, alpha){

  tau   <- fit$threshold
  phi   <- fit$parameters[,1]
  sigma <- fit$parameters[,2]
  xi    <- fit$parameters[,3]

  var <- tau + sigma/xi*((phi/(1-alpha))^xi-1)
  es  <- var/(1-xi) + (sigma-xi*tau)/(1-xi)

  out <- cbind(var, es)

  return(out)

}
