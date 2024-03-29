% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{L4}
\alias{L4}
\title{4-parameter logistic dose-response function}
\usage{
L4(dose, b, L, U, logEC50)
}
\arguments{
\item{dose}{Dose level}

\item{b}{Hill's coefficient (slope of the curve)}

\item{L}{Baseline effect (at zero dose)}

\item{U}{Asymptote effect (at infinite dose)}

\item{logEC50}{Point of inflection (in logarithmic terms)}
}
\description{
4-parameter logistic dose-response function
}
