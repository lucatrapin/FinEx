#' Autoregressive quantile model
#'
#' Fits an autoregressive quantile model with realized measures.
#'
#' @param r Vector returns on an asset.
#' @param p Probability level of the quantile (scalar).
#' @param x Vector of observations for a realized measure.
#' @param model Type of model. "s" is default. See the details.
#' @param sv Vector of starting values. This is an optional argument.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{estimates}}{ Vector of estimated parameters.}
#' \item{\code{value}}{ Minimized value of the loss function.}
#' \item{\code{q}}{ Vector of fitted quantiles.}
#' }
#'
#' @details
#' Parameters are estimated minimizing the quantile loss function. Optimization is performed with the \code{optim()} function
#' alternating "BFGS" and "Nelder-Mead" algorithms until convergence as in Engle and Manganelli (2003).
#'
#' \code{model} currently allows to select one quantile models: "s". For a quantile \eqn{q^{p}} at level \eqn{p} the model is
#' \deqn{q^{p}_t = b_0 + b_1 q^{p}_{t-1} + b_2 x_{t-1}} where \eqn{x_{t-1}} is the lagged value of the realized measure.
#'
#' @examples
#' # Fit the ARQ model
#' p   <- 0.95
#' fit <- arq(sp500$r, p, sp500$rv)
#'
#' @references
#' Engle, R. F., and Manganelli, S. (2004). CAViaR: Conditional autoregressive value at risk by regression quantiles. \emph{Journal of Business & Economic Statistics}, 22(4), 367-381.
#'
#' @seealso
#' \code{\link{arq.forecast}}
#' @author Luca Trapin
#'
#' @export arq
arq <- function(r, p, x, model="s", sv = NULL){

  # Check arguments
  if (p >=1 || p <=0) {stop("Error: p must be betweeen (0,1)")}

  ### Estimate parameters ###
  n   <- length(r)
  q   <- stats::quantile(r, p)
  if (is.null(sv)) sv <- arq.sv(r, p, x)
  est <- stats::optim(sv, arq.llk, , r, p, x)

  old     <- est$par
  old.val <- est$value
  old.est <- est

  conv <- 1

  while (conv==1){

    est1 <- stats::optim(old, arq.llk, , r, p, x, method="BFGS")

    est <- stats::optim(est1$par, arq.llk, , r, p, x)

    if (est$value < old.val){
      old <- est$par
      old.val <- est$value
    } else {
      break
    }
    old.est <- est
    conv <- est$convergence
  }

  # output
  par <- est$par
  q <- rep(0, n)
  for (i in 2:n){
    q[i] <- par[1] + par[2]*q[i-1] + par[3]*x[i-1]
  }

  out <- list()

  out$estimates <- par
  names(out$estimates) <- c("b0", "b1", "b2")

  out$value <- est$value

  out$q <- q

  return(out)

}
