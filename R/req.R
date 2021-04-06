#' Realized Extreme Quantile
#'
#' Computes Value at Risk (VaR) using a dynamic threshold estimated with the \code{arq()} function
#' and extreme value refinements as in Bee et al. (2018).
#'
#' @param r Vector of returns on an asset.
#' @param p Probability level of the dynamic threshold.
#' @param x Vector of observations for a realized measure.
#' @param alpha Probability level of the VaR.
#' @param pa Probability level used to set the threshold in the Hill estimator.
#' @param model Type of model: "s" is default. See details in the \code{arq()} function.
#' @param sv Vector of starting values for the \code{arq()} model. This is an optional argument.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{estimates}}{Parameters of the \code{arq()} model and the tail index.}
#' \item{\code{value}}{Minimized value of the loss function from the \code{arq()} model.}
#' \item{\code{p}}{Probability level of the dynamic threshold.}
#' \item{\code{zk}}{Threshold value of the quantile residuals used to estimate the tail index.}
#' \item{\code{q}}{Vector containing the dynamic threshold computed with the \code{arq()} model.}
#' \item{\code{var}}{Vector containing the VaR estimates.}
#' }
#'
#' @details
#' Let \eqn{q^{p}} be the estimated dynamic threshold with the \code{arq()} model, the tail index (\eqn{\xi}) is computed on the quantile residuals, \eqn{z^{p}=\frac{r}{q^{p}},} using the the Hill estimator,
#' \deqn{\widehat{\xi}=\frac{1}{k} \sum^k_{j=1}\log\left(\frac{z^{p}_{(j)}}{z^{p}_{(k)}}\right),} where \eqn{z^{p}_{(k)}} is the order statistics associated to the probability level \eqn{pa}.
#' The VaR is computed as \deqn{VaR^{\alpha}_t=q^{p}_t z^{p}_k \left(\frac{p}{\alpha}\right)^{\xi}}
#'
#' @examples
#' # Fit the REQ model
#' p <- 0.95
#' pa <- 0.975
#' alpha <- 0.99
#' fit <- req(sp500$r, p, sp500$rv, alpha, pa)
#'
#' @references
#' Bee, M., Dupuis, D. J., and Trapin, L. (2018). Realized extreme quantile: A joint model for conditional quantiles and measures of volatility with EVT refinements. \emph{Journal of Applied Econometrics}, 33(3), 398-415.
#'
#' @seealso
#' \code{\link{arq}}, \code{\link{req.forecast}}
#'
#' @author Luca Trapin
#'
#' @export req
req <- function(r, p, x, alpha, pa, model="s", sv=NULL){

  fit <- arq(r, p, x, model, sv)

  z  <- r/fit$q
  ev <- hill(z, pa)
  xi <- ev$xi
  zk <- ev$zk

  out <- list()

  out$estimates <- c(fit$estimates, xi)
  names(out$estimates)[4] <- "xi"
  out$value <- fit$value
  out$p <- p
  out$zk <- zk
  out$q <- fit$q
  out$var <- fit$q*zk*(p/alpha)^xi


  return(out)

}

