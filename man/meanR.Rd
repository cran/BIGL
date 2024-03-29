% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanR.R
\name{meanR}
\alias{meanR}
\title{Compute meanR statistic for the estimated model}
\usage{
meanR(
  data_off,
  fitResult,
  transforms = fitResult$transforms,
  null_model = c("loewe", "hsa", "bliss", "loewe2"),
  R,
  CP,
  reps,
  nested_bootstrap = FALSE,
  B.B = NULL,
  B.CP = NULL,
  cl = NULL,
  method = c("equal", "model", "unequal"),
  bootStraps,
  paramsBootstrap,
  idUnique,
  n1,
  transFun,
  invTransFun,
  ...
)
}
\arguments{
\item{data_off}{data frame with off -axis information}

\item{fitResult}{Monotherapy (on-axis) model fit, e.g. produced by
\code{\link{fitMarginals}}. It has to be a \code{"MarginalFit"} object or a
list containing \code{df}, \code{sigma}, \code{coef},
\code{shared_asymptote} and \code{method} elements for, respectively,
marginal model degrees of freedom, residual standard deviation, named
vector of coefficient estimates, logical value of whether shared asymptote
is imposed and method for estimating marginal models during bootstrapping
(see \code{\link{fitMarginals}}). If biological and power transformations
were used in marginal model estimation, \code{fitResult} should contain
\code{transforms} elements with these transformations. Alternatively, these
can also be specified via \code{transforms} argument.}

\item{transforms}{Transformation functions. If non-null, \code{transforms} is
a list containing 5 elements, namely biological and power transformations
along with their inverse functions and \code{compositeArgs} which is a list
with argument values shared across the 4 functions. See vignette for more
information.}

\item{null_model}{Specified null model for the expected response surface.
Currently, allowed options are \code{"loewe"} for generalized Loewe model,
\code{"hsa"} for Highest Single Agent model, \code{"bliss"} for Bliss additivity,
and \code{"loewe2"} for the alternative Loewe generalization.}

\item{R}{Numeric vector containing mean deviation of predicted response
surface from the observed one at each of the off-axis points. If missing,
it will be calculated automatically from output of
\code{\link{predictOffAxis}} function.}

\item{CP}{Matrix which is part of covariance matrix for the \code{R} argument}

\item{reps}{Numeric vector containing number of replicates for each off-axis
dose combination. If missing, it will be calculated automatically from output
of \code{\link{predictOffAxis}} function.}

\item{nested_bootstrap}{When statistics are calculated, if
\code{nested_bootstrap = TRUE}, \code{CP} matrix is recalculated at each
bootstrap iteration of \code{B.B} using \code{B.CP} iterations. Using such
nested bootstrap may however significantly increase computational time. If
\code{nested_bootstrap = FALSE}, \code{CP} bootstrapped data reuses
\code{CP} matrix calculated from the original data.}

\item{B.B}{Number of iterations to use in bootstrapping null distribution for
either meanR or maxR statistics.}

\item{B.CP}{Number of bootstrap iterations to use for CP matrix estimation}

\item{cl}{If parallel computations are desired, \code{cl} should be a cluster
object created by \code{\link[parallel]{makeCluster}}. If parallel
computing is active, progress reporting messages are not necessarily
ordered as it should be expected.}

\item{method}{What assumption should be used for the variance of on- and 
off-axis points. This argument can take one of the values from 
\code{c("equal", "model", "unequal")}. With the value \code{"equal"} as the
default. \code{"equal"} assumes that both on- and off-axis points have the 
same variance, \code{"unequal"} estimates a different parameter for on- and 
off-axis points and \code{"model"} predicts variance based on the average 
effect of an off-axis point. If no transformations are used the 
\code{"model"} method is recommended. If transformations are used, only the
\code{"equal"} method can be chosen.}

\item{bootStraps}{precomputed bootstrap objects}

\item{paramsBootstrap}{parameters for the nested bootstrap}

\item{idUnique}{unique combinations of on-axis points, a character vector}

\item{n1}{the number of off-axis points}

\item{transFun, invTransFun}{the transformation and inverse transformation functions for the variance}

\item{...}{Further arguments that will be later passed to
\code{\link{generateData}} function during bootstrapping}
}
\value{
This function returns a \code{meanR} object with estimates for the
  meanR statistical test. \code{meanR} object is essentially a list with
  appropriately named elements.

  \code{meanR} object list includes notably the calculated F-statistic,
  p-value and degrees of freedom (\code{"n1"} and \code{"df0"} respectively)
  used to find the critical value of the F-distribution under the null.

  If \code{\link{meanR}} test is run with bootstrapping, then p-value
  estimate is based on bootstrapped null distribution of test statistic and an
  additional element \code{"FDist"} (of class \code{"ecdf"}) is returned.
}
\description{
\code{\link{meanR}} computes the meanR statistic for the provided model
and returns the computed F-statistic and the estimated p-value. p-value
can be calculated either by assuming an exact distribution or using
bootstrapping procedure. In the latter case, null distribution of
bootstrapped F-statistics is also returned.
}
