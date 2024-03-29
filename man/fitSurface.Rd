% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitSurface.R
\name{fitSurface}
\alias{fitSurface}
\title{Fit response surface model and compute meanR and maxR statistics}
\usage{
fitSurface(
  data,
  fitResult,
  transforms = fitResult$transforms,
  null_model = c("loewe", "hsa", "bliss", "loewe2"),
  effect = "effect",
  d1 = "d1",
  d2 = "d2",
  statistic = c("none", "meanR", "maxR", "both"),
  CP = NULL,
  B.CP = 50,
  B.B = NULL,
  nested_bootstrap = FALSE,
  error = 4,
  sampling_errors = NULL,
  wild_bootstrap = FALSE,
  wild_bootType = "normal",
  control = "FWER",
  cutoff = 0.95,
  parallel = FALSE,
  progressBar = TRUE,
  method = c("equal", "model", "unequal"),
  confInt = TRUE,
  digits = 9,
  bootRS = TRUE,
  trans = "identity",
  rescaleResids = FALSE,
  invtrans = switch(trans, identity = "identity", log = "exp"),
  newtonRaphson = FALSE,
  asymptotes = 2,
  bootmethod = method
)
}
\arguments{
\item{data}{Dose-response dataframe.}

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

\item{effect}{Name of the response column in the data ("effect")}

\item{d1}{Name of the column with doses of the first compound ("d1")}

\item{d2}{Name of the column with doses of the second compound ("d2")}

\item{statistic}{Which statistics should be computed. This argument can take
one of the values from \code{c("none", "meanR", "maxR", "both")}.}

\item{CP}{Prediction covariance matrix. If not specified, it will be estimated
by bootstrap using \code{B.CP} iterations.}

\item{B.CP}{Number of bootstrap iterations to use for CP matrix estimation}

\item{B.B}{Number of iterations to use in bootstrapping null distribution for
either meanR or maxR statistics.}

\item{nested_bootstrap}{When statistics are calculated, if
\code{nested_bootstrap = TRUE}, \code{CP} matrix is recalculated at each
bootstrap iteration of \code{B.B} using \code{B.CP} iterations. Using such
nested bootstrap may however significantly increase computational time. If
\code{nested_bootstrap = FALSE}, \code{CP} bootstrapped data reuses
\code{CP} matrix calculated from the original data.}

\item{error}{Type of error for resampling in the bootstrapping procedure.
This argument will be passed to \code{\link{generateData}}. If \code{error
= 4} (default), the error terms for generating distribution of the null
will be resampled from the vector specified in \code{sampling_errors}. If
\code{error = 1}, normal errors are added. If \code{error = 2}, errors are
sampled from a mixture of two normal distributions. If \code{error = 3},
errors are generated from a rescaled chi-square distribution.}

\item{sampling_errors}{Sampling vector to resample errors from. Used only if
\code{error} is 4 and is passed as argument to \code{\link{generateData}}.
If \code{sampling_errors = NULL} (default), mean residuals at off-axis
points between observed and predicted response are taken.}

\item{wild_bootstrap}{Whether special bootstrap to correct for
heteroskedasticity should be used. If \code{wild_bootstrap = TRUE}, errors
are generated from \code{sampling_errors} multiplied by a random variable
following Rademacher distribution. Argument is used only if \code{error = 4}.}

\item{wild_bootType}{Type of distribution to be used for wild bootstrap. If \code{wild_bootstrap = TRUE}, 
errors are generated from "rademacher", "gamma", "normal" or "two-point" distribution.}

\item{control}{If \code{control = "FCR"} then algorithm controls false coverage rate, if \code{control = "dFCR"} then 
algorithm controls directional false coverage rate, if \code{control = "FWER"} then 
algorithm controls family wise error rate}

\item{cutoff}{Cut-off to use in maxR procedure for declaring non-additivity
(default is 0.95).}

\item{parallel}{Whether parallel computing should be used for bootstrap. This
parameter can take either integer value to specify the number of threads to
be used or logical \code{TRUE/FALSE}. If \code{parallel = TRUE}, then
\code{max(1, detectCores()-1)} is set to be the number of threads. If
\code{parallel = FALSE}, then a single thread is used and cluster object
is not created.}

\item{progressBar}{A boolean, should progress of bootstraps be shown?}

\item{method}{What assumption should be used for the variance of on- and 
off-axis points. This argument can take one of the values from 
\code{c("equal", "model", "unequal")}. With the value \code{"equal"} as the
default. \code{"equal"} assumes that both on- and off-axis points have the 
same variance, \code{"unequal"} estimates a different parameter for on- and 
off-axis points and \code{"model"} predicts variance based on the average 
effect of an off-axis point. If no transformations are used the 
\code{"model"} method is recommended. If transformations are used, only the
\code{"equal"} method can be chosen.}

\item{confInt}{a boolean, should confidence intervals be returned?}

\item{digits}{Numeric value indicating the number of digits used for numeric values in confidence intervals}

\item{bootRS}{a boolean, should bootstrapped response surfaces be used in the
calculation of the confidence intervals?}

\item{trans, invtrans}{the transformation function for the variance and its
inverse, possibly as strings}

\item{rescaleResids}{a boolean indicating whether to rescale residuals,
or else normality of the residuals is assumed.}

\item{newtonRaphson}{A boolean, should Newton-Raphson be used to find Loewe
response surfaces? May be faster but also less stable to switch on}

\item{asymptotes}{Number of asymptotes. It can be either \code{1}
as in standard Loewe model or \code{2} as in generalized Loewe model.}

\item{bootmethod}{The resampling method to be used in the bootstraps. Defaults to the same as method}
}
\value{
This function returns a \code{ResponseSurface} object with estimates
  of the predicted surface. \code{ResponseSurface} object is essentially a
  list with appropriately named elements.

  Elements of the list include input data, monotherapy model coefficients and
  transformation functions, null model used to construct the surface as well
  as estimated CP matrix, occupancy level at
  each dose combination according to the generalized Loewe model and
  \code{"offAxisTable"} element which contains observed and predicted effects
  as well as estimated z-scores for each dose combination.

  If statistical testing was done, returned object contains \code{"meanR"}
  and \code{"maxR"} elements with output from \code{\link{meanR}} and
  \code{\link{maxR}} respectively.
}
\description{
This function computes predictions for off-axis dose combinations according
to the BIGL or HSA null model and, if required, computes appropriate meanR
and maxR statistics. Function requires as input dose-response dataframe and
output of \code{\link{fitMarginals}} containing estimates for the monotherapy
model. If transformation functions were used in monotherapy estimation, these
should also be provided.
}
\details{
Please see the example vignette \code{vignette("analysis", package = "BIGL")}
and the report "Lack of fit test for detecting synergy" included in the
\code{papers} folder for further details on the test statistics used:
\code{system.file("papers", "newStatistics.pdf", package = "BIGL")}
}
\examples{
\dontrun{
  data <- subset(directAntivirals, experiment == 4)
  ## Data should contain d1, d2 and effect columns
  transforms <- list("PowerT" = function(x, args) with(args, log(x)),
                     "InvPowerT" = function(y, args) with(args, exp(y)),
                     "BiolT" = function(x, args) with(args, N0 * exp(x * time.hours)),
                     "InvBiolT" = function(y, args) with(args, 1/time.hours * log(y/N0)),
                     "compositeArgs" = list(N0 = 1, time.hours = 72))
  fitResult <- fitMarginals(data, transforms)
  surf <- fitSurface(data, fitResult, statistic = "meanR")
  summary(surf)
}
}
