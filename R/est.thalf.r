#' Calculate Lambda_z and Elimination Half-life
#'
#' Calculates lambda_z and thalf for each PK curve identified using \code{by}. \cr
#'
#' The function starts with the last three sample points and performs
#' log-linear regression on it. It then adds one sampling point at a time
#' (including and ending at tmax) and performs the regression again.
#' The results of the regression with the highest adjusted R-squared are returned. \cr
#' \cr
#' Visual outliers can be excluded from the regression analysis.
#'
#' @param x a dataset
#' @param by character: column names in x indicating grouping variables; default is as.character(dplyr::groups(x))
#' @param timevar variable name containing the actual sampling time after dose
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param includeCmax include results of regression including Cmax in selection? (y/n); x$includeCmax overrides if provided
#' @param exclvar a variable name containing information about points to be excluded (these should have exclvar = 1)
#' @param startth user defined start time for regression
#' @param endth user defined end time for regression
#' @importFrom stats lm
#' @importFrom dplyr first last
#' @return
#' a dataset with estimates for each regression analysis in one observation.
#' The following parameters are available.
#' * **no.points** number of data points used in the regression analysis
#' * **intercept** estimated intercept
#' * **lambda_z** -1*estimated slope
#' * **r.squared** square of the correlation coefficient
#' * **adj.r.squared** adjusted square of the correlation coefficient
#' * **thalf** elimination half-life
#' * **start_th** time of first sample included in the thalf estimation
#' * **end_th time** of last sample included in the thalf estimation
#' * **includeCmax** include results of regression including Cmax in selection? (y/n)
#' * **points_excluded** are time points excluded from the half-life estimation? (y/n)
#' * **usrmod** are start and/or end point of regression interval user defined? (y/n)
#' @export
#' @examples
#' \donttest{
#' library(dplyr)
#' data(ncx)
#' x <- ncx
#' x <- x |> group_by(subject)
#' x <- x |> correct.loq
#' x |> est.thalf |> head
#' }

est.thalf <- function(
  x,
  by = NULL,
  timevar = "time",
  depvar = "dv",
  includeCmax = "Y",
  exclvar = NA,
  startth = NA,
  endth = NA
) {
  if (is.null(by)) by <- as.character(groups(x))
  supplied <- character(0)
  if (!missing(includeCmax)) supplied <- c(supplied, 'includeCmax')
  if (!missing(startth)) supplied <- c(supplied, 'startth')
  if (!missing(endth)) supplied <- c(supplied, 'endth')
  x <- group_by_at(x, vars(by))
  x <- do(
    .data = x,
    .est.thalf(
      .,
      timevar = timevar,
      depvar = depvar,
      includeCmax = includeCmax,
      exclvar = exclvar,
      startth = startth,
      endth = endth,
      supplied = supplied
    )
  )
  # x <- ungroup(x)
  x
}
.est.thalf <- function(
  x,
  timevar,
  depvar,
  includeCmax,
  exclvar,
  startth,
  endth,
  supplied
) {
  if ('includeCmax' %in% names(x)) {
    if ('includeCmax' %in% supplied) {
      warning('includeCmax supplied as column overrides like-named argument')
    }
    includeCmax <- unique(x$includeCmax)
    x$includeCmax <- NULL
  }
  if (length(includeCmax) > 1) {
    warning('includeCmax has length > 1; only first value will be used')
    includeCmax <- includeCmax[[1]]
  }

  if ('startth' %in% names(x)) {
    if (!is.na(unique(x$startth))) {
      if ('startth' %in% supplied) {
        warning('startth supplied as column overrides like-named argument')
      }
      startth <- unique(x$startth)
      x$startth <- NULL
    }
    x$startth <- NULL
  }
  if (length(startth) > 1) {
    warning('startth has length > 1; only first value will be used')
    startth <- startth[[1]]
  }

  if ('endth' %in% names(x)) {
    if (!is.na(unique(x$endth))) {
      if ('endth' %in% supplied) {
        warning('endth supplied as column overrides like-named argument')
      }
      endth <- unique(x$endth)
      x$endth <- NULL
    }
    x$endth <- NULL
  }
  if (length(endth) > 1) {
    warning('endth has length > 1; only first value will be used')
    endth <- endth[[1]]
  }

  if (!'timevar' %in% names(x)) x <- x |> rename(timevar = !!timevar)
  if (!'depvar' %in% names(x)) x <- x |> rename(depvar = !!depvar)

  if (!is.na(exclvar) & !(exclvar %in% names(x)))
    stop(paste("Exclusion variable", exclvar, "does not exist"), call. = F)

  if (!(is.na(exclvar)) & exclvar %in% names(x)) {
    x <- x |> rename(exclvar = !!exclvar)
  }

  if (!is.na(exclvar)) {
    usrmod = 0 # flags whether user set start and/or stop time of range
    anyexcl = 0
    if (any(x$exclvar == 1)) {
      anyexcl = 1
    }

    x <- x |> #      mutate(tmax=first(timevar[depvar==max(depvar,na.rm=T)])) |>
      filter(exclvar != 1 | is.na(exclvar)) # remove samples to be excluded from the regression
  }

  x <- x |> filter(!is.na(depvar) & depvar > 0) |>
    mutate(tmax = first(timevar[depvar == max(depvar)])) |>
    filter(timevar >= tmax) # greater or equal than ORIGINAL tmax (before possible exclusions)

  if (tolower(includeCmax) == "n") {
    x = x |> filter(timevar > tmax) # greater than ORIGINAL tmax (before possible exclusions)
  }

  # 1 both start and end: use these two, no optimalisation

  if (!is.na(startth) & !is.na(endth)) {
    x <- x |> filter(timevar >= startth & timevar <= endth) # keep only samples between startth and endth

    usrmod = 1
    est = 1
    result = data.frame(matrix(ncol = 7, nrow = 1))

    lmcall = lm(log(x$depvar) ~ x$timevar)
    lmcall.estimates = as.numeric(summary(lmcall)$coef[, "Estimate"])
    ## save results of loglin lm fit
    result[1] = length(x$timevar)
    result[2] = exp(lmcall.estimates[1]) # exponentiate to see the actual intercept
    result[3] = lmcall.estimates[2] * -1
    result[4] = summary(lmcall)$r.squared
    result[5] = summary(lmcall)$adj.r.squared
    result[6] = startth
    result[7] = endth
  } # 2 only start: remove all before that start time, then start optimalisation
  else {
    if (!is.na(startth) & is.na(endth)) {
      usrmod = 1
      x <- x |> filter(timevar >= startth) # keep only samples from startth to the end
    }

    # 3 only end: remove all after that end time, then start optimalisation

    if (is.na(startth) & !is.na(endth)) {
      usrmod = 1
      x <- x |> filter(timevar <= endth) # keep only samples from start to endth
    }

    est = 0
    i = length(x$timevar) - 2
    result = data.frame(matrix(ncol = 7, nrow = 1))

    if (i >= 1) {
      est = 1
      result = data.frame(matrix(ncol = 7, nrow = i))
    }

    while (i >= 1) {
      ## make subset of data frame for each number of data points
      xx = x[i:nrow(x), ]
      ## execute loglinear model fit
      lmcall = lm(log(xx$depvar) ~ xx$timevar)
      lmcall.estimates = as.numeric(summary(lmcall)$coef[, "Estimate"])
      ## save results of loglin lm fit
      result[i, 1] = length(x$timevar) - i + 1
      result[i, 2] = exp(lmcall.estimates[1]) # exponentiate to see the actual intercept
      result[i, 3] = lmcall.estimates[2] * -1
      result[i, 4] = summary(lmcall)$r.squared
      result[i, 5] = summary(lmcall)$adj.r.squared
      result[i, 6] = x$timevar[i]
      result[i, 7] = last(x$timevar)
      i = i - 1
    }
  }
  names(result) = c(
    'no.points',
    'intercept',
    'lambda_z',
    'r.squared',
    'adj.r.squared',
    'start_th',
    'end_th'
  )
  if (est == 1) {
    result = result |>
      mutate(
        sel = no.points[adj.r.squared == max(adj.r.squared)],
        thalf = log(2) / lambda_z,
        includeCmax = includeCmax
      ) |>
      filter(sel == no.points) |>
      select(-sel)
  } else {
    result = result |>
      mutate(
        no.points = NA,
        intercept = NA,
        lambda_z = NA,
        r.squared = NA,
        adj.r.squared = NA,
        thalf = NA,
        start_th = NA,
        end_th = NA
      )
  }

  if (!is.na(exclvar)) {
    if (anyexcl == 1) {
      result = result |>
        mutate(points_excluded = "Y")
    } else {
      result = result |>
        mutate(points_excluded = "N")
    }
  } else {
    result = result |>
      mutate(points_excluded = "N")
  }

  result = result |>
    mutate(usrmod = ifelse(usrmod == 1, "Y", "N"))

  return(result)
}
