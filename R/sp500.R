#' Daily returns and 5-minute realized variance for the SP500.
#'
#' Daily observations from 2000 to 2014 from the Oxford-Man Realized library (Heber et al., 2009). These time series were used in the
#' empirical analysis by Bee and Trapin (2018).
#'
#' @docType data
#'
#' @usage data(sp500)
#'
#' @references
#' Heber, G., Lunde, A., Shephard, N., and Sheppard, K. (2009). \emph{Oxford-Man Institute’s realized library}, version 0.1.
#'
#' Bee, M., and Trapin, L. (2018). Estimating and forecasting conditional risk measures with extreme value theory: A review. \emph{Risks}, 6(2), 45.
#'
#' @examples
#' data(sp500)
#' returns <- sp500$r
#' realized_variance <- sp500$rv
"sp500"
