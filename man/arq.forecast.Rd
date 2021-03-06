% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arq.forecast.R
\name{arq.forecast}
\alias{arq.forecast}
\title{Dynamic Quantile forecasts}
\usage{
arq.forecast(fit, x, alpha)
}
\arguments{
\item{fit}{Output object of the \code{arq()} function.}

\item{x}{Scalar with the realized measures used in the forecast.}

\item{alpha}{Probability level for the computation of the risk measures.}
}
\value{
The quantile forecast at level alpha.
}
\description{
Computes one-day-ahead forecasts of the dynamic quantile.
}
\author{
Luca Trapin
}
