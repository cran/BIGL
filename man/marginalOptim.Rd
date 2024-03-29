% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitMarginals.R
\name{marginalOptim}
\alias{marginalOptim}
\title{Fit two 4-parameter log-logistic functions with common baseline}
\usage{
marginalOptim(data, transforms = NULL, start, model, ...)
}
\arguments{
\item{data}{Dose-response dataframe. Marginal data will be extracted from
it automatically.}

\item{transforms}{Transformation functions. If non-null, \code{transforms} is
a list containing 5 elements, namely biological and power transformations
along with their inverse functions and \code{compositeArgs} which is a list
with argument values shared across the 4 functions. See vignette for more
information.}

\item{start}{Starting parameter values. If not specified, they will be
obtained from \code{\link{initialMarginal}}.}

\item{model}{List with model parameters. Typically, this is an output from
\code{\link{constructFormula}}.}

\item{...}{Further parameters passed to \code{\link[stats]{optim}} function}
}
\value{
Variance-covariance matrix which is returned by \code{\link{optim}}
  is based on the fact that minimization of sum-of-squared residuals leads
  essentially to a maximum likelihood estimator and so variance-covariance
  matrix can be estimated using inverse Hessian evaluated at the optimal
  parameters. In some cases, so obtained variance-covariance matrix might not
  be positive-definite which probably means that estimates are unstable
  because of either a poor choice of initial values or poor properties of the
  data itself.
}
\description{
This function is an alternative to non-linear least squares and
provides optimization framework with \code{\link{optim}} function.
It is however noticeably slower than NLS methods and can be especially
time consuming in large datasets, in particular if bootstrap statistics
are calculated.
}
