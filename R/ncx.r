#' NCA Example Dataset
#'
#' Sample data formatted to support examples in the qpNCA package.
#' Derived from 'Theoph'.
#'
#' @format A data frame with 132 rows and 10 variables:
#' \describe{
#'   \item{subject}{subject identifier}
#'   \item{wt}{Weight (kg)}
#'   \item{dose}{dose (mg)}
#'   \item{time}{time since drug administration}
#'   \item{tad}{time after dose (h)}
#'   \item{ntad}{nominal time after dose (h)}
#'   \item{dv}{drug concentration (mg/L)}
#'   \item{bloq}{below limit of quantitation (1) or not (0)}
#'   \item{loq}{limit of quantitation (mg/L)}
#'   \item{loqrule}{loq rule}
#' }
#' @usage data(ncx)
#' @source
#' library(magrittr)
#'
#' library(dplyr)
#'
#' library(qpNCA)
#'
#' x <- Theoph
#'
#' ntad <- c(0,0.25,0.5,1,2,4,5,7,9,12,24)
#'
#' for(i in 1:nrow(x))\{
#'
#'   time  <- x$Time\[\[i]]
#'
#'   delta <- abs(ntad - time)
#'
#'   best  <- min(delta)
#'
#'   index <- match(best, delta)
#'
#'   nom   <- ntad\[\[index]]
#'
#'   x$ntad\[\[i]] <- nom
#'
#' \}
#'
#' rm(list = c('time','delta','best','index','nom', 'i','ntad'))
#'
#' x %<>% rename(time = Time, dv = conc)
#'
#' x %<>% mutate(bloq = ifelse(dv==0,1,0), loq = 0.01, tad = time, loqrule = 1,
#'               subject=as.numeric(Subject), ntad=as.numeric(ntad))
#'
#' x %<>% select(-subject)
#'
#' x %<>% mutate(Dose = signif(digits = 3, Dose * Wt))
#'
#' names(x) %<>% tolower
#'
#' x %<>% select(subject, wt, dose, time, tad, ntad, dv, everything())
#'
"ncx"
