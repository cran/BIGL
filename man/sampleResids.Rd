% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaleResids.R
\name{sampleResids}
\alias{sampleResids}
\title{Sample residuals according to a new model}
\usage{
sampleResids(means, sampling_errors, method, rescaleResids, ...)
}
\arguments{
\item{means}{a vector of means}

\item{sampling_errors}{Sampling vector to resample errors from. Used only if
\code{error} is 4 and is passed as argument to \code{\link{generateData}}.
If \code{sampling_errors = NULL} (default), mean residuals at off-axis
points between observed and predicted response are taken.}

\item{method}{What assumption should be used for the variance of on- and 
off-axis points. This argument can take one of the values from 
\code{c("equal", "model", "unequal")}. With the value \code{"equal"} as the
default. \code{"equal"} assumes that both on- and off-axis points have the 
same variance, \code{"unequal"} estimates a different parameter for on- and 
off-axis points and \code{"model"} predicts variance based on the average 
effect of an off-axis point. If no transformations are used the 
\code{"model"} method is recommended. If transformations are used, only the
\code{"equal"} method can be chosen.}

\item{rescaleResids}{a boolean indicating whether to rescale residuals,
or else normality of the residuals is assumed.}

\item{...}{passed on to predictVar}
}
\value{
sampled residuals
}
\description{
Sample residuals according to a new model
}
