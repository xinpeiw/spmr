% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two-stage-prediction.R
\name{sp_para}
\alias{sp_para}
\title{Conduct the two-stage prediction estimation}
\usage{
sp_para(d)
}
\arguments{
\item{d}{output of the function format_sp()}
}
\value{
the coefficient and corresponding standard error of every causal function f_j(X).

the F-value and corresponding P-value for hypothesis testing H_0: f = 0 (the causal effect of the exposure on the outcome is zero).
}
\description{
conduct the two-stage prediction estimation.
}
