#' Realized Peaks over Threshold
#'
#' Fits a Realized Peaks over Threshold (RPoT) model as in Bee et al. (2019). Realized measures can be used to model the dynamics
#' of the exceedance probability (\eqn{\phi}), and the scale (\eqn{\sigma}) and tail (\eqn{\xi}) parameters of the Generalized Pareto distribution.
#'
#' @param r Vector of returns on an asset.
#' @param q A constant threshold (scalar).
#' @param x Matrix of realized measures. The intercept must be included, i.e. set the first column to a vector of ones.
#' @param model Type of model. "s" is default. See the details.
#' @param xp Select which realized measures (columns of x) appear in phi.
#' @param xs Select which realized measures (columns of x) appear in sigma.
#' @param xx Select which realized measures (columns of x) appear in xi.
#' @param lf Link function for phi, sigma and xi. See details.
#' @param sv Vector of starting values. This is an optional argument.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{model}}{Type of model.}
#' \item{\code{estimates}}{Vector of regression parameters: psi, gamma, delta.}
#' \item{\code{value}}{Maximized log-likelihood.}
#' \item{\code{parameters}}{Matrix containing the fitted values of phi, sigma, and xi.}
#' \item{\code{lf}}{Link function for phi, sigma and xi.}
#' \item{\code{threshold}}{Constant threshold (q).}
#' }
#'
#'@details
#' Parameters are estimated by Maximum Likelihood. Optimization is performedusing the "Nelder-Mead" algorithm of the optim() function. The joint likelihood function of the RPoT model is
#' \deqn{\prod^{T}_{t=1}(1-\phi_t)^{1-I_t}\left(\frac{\phi_t}{\sigma_t}\left(1+\frac{\xi_t}{\sigma_t}\right)^{-\frac{1}{\xi_t}-1}\right)^{I_t}}
#' where \eqn{I_t} is the indicator function of the exceedancees.
#'
#' \code{lf} is used to constraint the parameters to the appropriate parameter space. The link function can be either \code{"logit"} (logistic function), \code{"exp"} (exponential function) or
#' \code{"identity"} (identity function). The logistic is the only function that can be used for phi.
#' The exponential and identity functions can be used for sigma and xi. Default is \code{lf=c("logit","exp","identity")}.
#'
#' \code{model} allows to select two different models:
#' \itemize{
#' \item "s" is the simple RPoT where the dynamic parameters are fully characterized by the realized measures,
#' \deqn{\phi_t=f_p(\psi x_t), \quad \sigma_t=f_s(\gamma x_t), \quad \xi_t=f_x(\delta x_t)} with \eqn{f(\cdot)} the corresponding link function.
#' \item "ar" includes also an autoregressive component in the dynamic parameters. Let \eqn{\phi=f_p(\tilde{\phi})}, \eqn{\sigma=f_s(\tilde{\sigma})}, \eqn{\xi=f_x(\tilde{\xi})} then
#' \deqn{\tilde{\phi}_t=\psi_1 + \psi_{ar} \tilde{\phi}_{t-1} + \psi_{-} x_t,}
#' \deqn{\tilde{\sigma}_t=\gamma_1 + \gamma_{ar} \tilde{\sigma}_{t-1} + \gamma_{-} x_t,}
#' \deqn{\tilde{\xi}_t=\delta_1 + \delta_{ar} \tilde{\xi}_{t-1} + \delta_{-} x_t}
#' where \eqn{\psi_{-}}, \eqn{\gamma_{-}}, \eqn{\delta_{-}} are the remaining parameters in the \code{estimates} vector.
#' }
#'
#' @examples
#' # Define a lag structure for the return and the realized variane
#' n   <- nrow(sp500)
#' rx  <- sp500$r[-1]
#' rvx <- sp500$rv[-n]
#'
#' # Define the threshold
#' q <- quantile(rx, 0.9)
#'
#' # Define the covariate structure
#' x  <- cbind(1,log(rvx))
#' xp <- c(1,2)
#' xs <- c(1,2)
#' xx <- 1
#'
#' # Fit the RPOT model
#' fit <- rpot(rx, q, x, model="s", xp, xs, xx)
#'
#' @references
#' Bee, M., Dupuis, D. J., and Trapin, L. (2019). Realized Peaks over Threshold: A Time-Varying Extreme Value Approach with High-Frequency-Based Measures. \emph{Journal of Financial Econometrics}, 17(2), 254-283.
#'
#' @seealso
#' \code{\link{rpot.fitted}}, \code{\link{rpot.forecast}}
#'
#' @author Luca Trapin
#'
#' @export rpot
rpot <- function(r, q, x, model="s", xp, xs, xx, lf = c("logit", "exp", "identity"), sv=NULL){

  y <- ifelse(r>q, r-q, 0)

  # Check correct link function
  if(!is.matrix(x)){stop("Error: x must be a matrix")}
  if (lf[1]!="logit"){stop("Error: phi link function must be logit.")}
  if (lf[2]=="logit"){stop("Error: sigma link function must be either exp or identity")}
  if (lf[3]=="logit"){stop("Error: xi link function must be either exp or identity")}

  if (is.null(sv)) sv <- rpot.sv(y, x, model, xp, xs, xx, lf)

  fit <- stats::optim(sv, rpot.llk, gr = NULL, y,  x, model, xp, xs, xx, lf, hessian = FALSE, control = list(maxit=10000))

  out <- rpot.output(fit=fit, y=y, x=x, mod=model, xp=xp, xs=xs, xx=xx, lf=lf, q=q)

  return(out)

}
