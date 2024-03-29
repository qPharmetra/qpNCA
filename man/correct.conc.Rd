% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct.conc.r
\name{correct.conc}
\alias{correct.conc}
\title{Correct Missing Concentration}
\usage{
correct.conc(
  x,
  by = NULL,
  nomtimevar = "ntad",
  tau = NA,
  tstart = NA,
  tend = NA,
  teval = NA,
  reg = "SD",
  ss = "N",
  route = "EV",
  method = 1
)
}
\arguments{
\item{x}{input dataset name (after Time Deviation Correction Rules have been applied by \code{\link{correct.time}})}

\item{by}{character: column names in x indicating grouping variables; default is as.character(dplyr::groups(x))}

\item{nomtimevar}{variable name containing the nominal sampling time after dose}

\item{tau}{dosing interval (for multiple dosing); NA (default) for if single dose; x$tau overrides}

\item{tstart}{start time of partial AUC (start>0); NA (default) if not requested; x$tstart overrides}

\item{tend}{end time of partial AUC; NA (default) if not requested; x$tend overrides}

\item{teval}{user selected AUC interval; NA (default) if not requested; x$teval overrides}

\item{reg}{regimen, "sd" or "md"; x$reg overrides}

\item{ss}{is steady state reached (y/n); x$ss overrides}

\item{route}{route of drug administration ("EV","IVB","IVI"); x$route overrides}

\item{method}{method for trapezoidal rule;  x$method overrides
\itemize{
\item 1: linear up - linear down
\item 2: linear up - logarithmic down
\item 3: linear before first Tmax, logarithmic after first Tmax
}}
}
\value{
a dataset with missing concentrations imputed. The following variables are added:\tabular{ll}{
   \strong{Variable} \tab \strong{Description} \cr
   crule.nr \tab correction rule number \cr
   crule.txt \tab text explaining what was altered \cr
   applies.to.conc \tab lists all critical time points to which the concentration correction rule applies \cr
}
}
\description{
Corrects missing concentration at critical time points
(e.g, predose, TAU, start and end of user selected AUC interval).
}
\details{
\itemize{
\item Use interpolation if there is a measurable concentration BEFORE and AFTER the missing concentration
\item Use extrapolation if there is NO measurable concentration AFTER the missing concentration
\item Set missing concentration at predose to 0 (SD, non-endogenous) or value at t=TAU (steady state only)
\item Set missing concentration at t=TAU to value at t=0 (steady state only)
}

The following Concentration Deviation Correction Rules will be applied to critical time points (t=0, tau, tstart, tend, teval), if needed:\tabular{llll}{
   \strong{Rule} \tab \strong{Regimen} \tab \strong{Description} \tab \strong{Applied to} \cr
   SDC-1 \tab sd \tab Set concentration to 0 (only non-endogenous compounds) \tab t=0 \cr
   SDC-2 \tab sd \tab impute missing concentration by interpolation \tab t=tau,tstart,tend,teval \cr
   SDC-3 \tab sd \tab impute missing concentration by extrapolation \tab t=tau,tend,teval \cr
   SDC-4 \tab sd (IVB) \tab impute missing concentration by back-extrapolation \tab t=0 \cr
   MDC-1 \tab md \tab impute missing concentration by existing conc at t=0 or t=tau (only if steady state has been reached) \tab t=0,tau \cr
   MDC-2 \tab md \tab impute missing concentration by interpolation \tab t=tau,tstart,tend,teval \cr
   MDC-3 \tab md \tab impute missing concentration by extrapolation \tab t=tau,tend,teval \cr
   MDC-4 \tab md (IVB) \tab impute missing concentration by back-extrapolation \tab t=0 \cr
}
}
\examples{
\donttest{
library(magrittr)
library(dplyr)
data(ncx)
x <- ncx
x \%<>\% group_by(subject)
x \%<>\% correct.loq
x \%<>\% correct.time
x \%>\% correct.conc \%>\% head
}
}
