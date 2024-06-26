% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harbronLoewe.R
\name{harbronLoewe}
\alias{harbronLoewe}
\title{Alternative Loewe generalization}
\usage{
harbronLoewe(
  doseInput,
  parmInput,
  asymptotes = 2,
  startvalues = NULL,
  newtonRaphson = FALSE,
  ...
)
}
\arguments{
\item{doseInput}{Dose-response dataframe containing \code{"d1"} and
\code{"d2"} columns}

\item{parmInput}{Numeric vector or list with appropriately named
parameter inputs. Typically, it will be coefficients from a
\code{MarginalFit} object.}

\item{asymptotes}{Number of asymptotes. It can be either \code{1}
as in standard Loewe model or \code{2} as in generalized Loewe model.}

\item{startvalues}{Starting values for the non-linear equation,
from the observed data}

\item{newtonRaphson}{a boolean, is Newton raphson used for finding the
response surface? May be faster but also less stable}

\item{...}{Further arguments that are currently unused}
}
\description{
Alternative Loewe generalization
}
