#' Correct Missing Concentration
#'
#' Corrects missing concentration at critical time points
#' (e.g, predose, TAU, start and end of user selected AUC interval).
#'
#' * Use interpolation if there is a measurable concentration BEFORE and AFTER the missing concentration
#' * Use extrapolation if there is NO measurable concentration AFTER the missing concentration
#' * Set missing concentration at predose to 0 (SD, non-endogenous) or value at t=TAU (steady state only)
#' * Set missing concentration at t=TAU to value at t=0 (steady state only)
#'
#' The following Concentration Deviation Correction Rules will be applied to critical time points (t=0, tau, tstart, tend, teval), if needed:
#'
#' **Rule** | **Regimen** | **Description** | **Applied to**
#' ---- | ------- | ----------- | ----------
#' SDC-1  |   sd   |       Set concentration to 0 (only non-endogenous compounds)      |      t=0
#' SDC-2   |    sd  |      impute missing concentration by interpolation   |                  t=tau,tstart,tend,teval
#' SDC-3   |    sd    |    impute missing concentration by extrapolation   |                  t=tau,tend,teval
#' SDC-4   |    sd (IVB)  | impute missing concentration by back-extrapolation   |             t=0
#' MDC-1   |    md  |      impute missing concentration by existing conc at t=0 or t=tau (only if steady state has been reached)  |  t=0,tau
#' MDC-2   |    md    |    impute missing concentration by interpolation   |                   t=tau,tstart,tend,teval
#' MDC-3   |    md    |    impute missing concentration by extrapolation   |                  t=tau,tend,teval
#' MDC-4   |    md (IVB)  |  impute missing concentration by back-extrapolation  |              t=0
#'
#' @importFrom dplyr left_join lead
#' @param x input dataset name (after Time Deviation Correction Rules have been applied by \code{\link{correct.time}})
#' @param by character: column names in x indicating grouping variables; default is as.character(dplyr::groups(x))
#' @param nomtimevar variable name containing the nominal sampling time after dose
#' @param tau dosing interval (for multiple dosing); NA (default) for if single dose; x$tau overrides
#' @param tstart start time of partial AUC (start>0); NA (default) if not requested; x$tstart overrides
#' @param tend end time of partial AUC; NA (default) if not requested; x$tend overrides
#' @param teval user selected AUC interval; NA (default) if not requested; x$teval overrides
# @param th lamdba_z information for each curve; pass NA to suppress; default is est.thalf(x, by = by, timevar = timevar, depvar = depvar )
#' @param reg regimen, "sd" or "md"; x$reg overrides
#' @param ss is steady state reached (y/n); x$ss overrides
#' @param route route of drug administration ("EV","IVB","IVI"); x$route overrides
#' @param method method for trapezoidal rule;  x$method overrides

#' * 1: linear up - linear down
#' * 2: linear up - logarithmic down
#' * 3: linear before first Tmax, logarithmic after first Tmax
#'
#' @return  a dataset with missing concentrations imputed. The following variables are added:
#'
#' **Variable** | **Description**
#' -------- | -----------
#'  crule.nr     |   correction rule number
#'  crule.txt     |  text explaining what was altered
#'  applies.to.conc  | lists all critical time points to which the concentration correction rule applies
#' @export
#' @examples
#' \donttest{
#' library(magrittr)
#' library(dplyr)
#' data(ncx)
#' x <- ncx
#' x %<>% group_by(subject)
#' x %<>% correct.loq
#' x %<>% correct.time
#' x %>% correct.conc %>% head
#' }
correct.conc <- function(
  x,
  by = NULL,
  nomtimevar = "ntad",
  tau = NA,
  tstart = NA,
  tend = NA,
  teval = NA,
  #  th=NA,
  reg = "SD",
  ss = "N",
  route = "EV",
  method = 1
) {
  if (is.null(by)) by <- as.character(groups(x))
  supplied <- character(0)
  if (!missing(tau)) supplied <- c(supplied, 'tau')
  if (!missing(tstart)) supplied <- c(supplied, 'tstart')
  if (!missing(tend)) supplied <- c(supplied, 'tend')
  if (!missing(teval)) supplied <- c(supplied, 'teval')
  if (!missing(reg)) supplied <- c(supplied, 'reg')
  if (!missing(method)) supplied <- c(supplied, 'method')
  if (!missing(route)) supplied <- c(supplied, 'route')
  if (!missing(ss)) supplied <- c(supplied, 'ss')
  if (!'conc.tau' %in% names(x))
    stop('conc.tau not found; probably you need to run correct.time()')
  x <- group_by_at(x, vars(by))
  x <- do(
    .data = x,
    .correct.conc(
      .,
      by = by,
      nomtimevar = nomtimevar,
      tau = tau,
      tstart = tstart,
      tend = tend,
      teval = teval,
      th = th,
      reg = reg,
      ss = ss,
      route = route,
      method = method,
      supplied = supplied
    )
  )
  # x <- ungroup(x)
  x
}

.correct.conc <- function(
  x,
  by,
  nomtimevar,
  tau,
  tstart,
  tend,
  teval,
  th,
  reg,
  ss,
  route,
  method,
  supplied
) {
  #tau=tau,tstart=tstart,tend=tend,teval=teval,reg=reg,ss=ss,route=route,method=method

  for (arg in c(
    'tau',
    'tstart',
    'tend',
    'teval',
    'reg',
    'method',
    'route',
    'ss'
  )) {
    if (arg %in% names(x)) {
      if (arg %in% supplied) {
        warning(arg, ' supplied as column overrides like-named argument')
      }
      assign(arg, unique(x[[arg]]))
      # x[[arg]] <- NULL
    }
    if (length(get(arg)) > 1) {
      warning(arg, ' has length > 1; only first value will be used')
      assign(arg, get(arg)[[1]])
    }
    if (arg %in% names(x)) {
      # assure x preserves value actually used
      x[[arg]] <- get(arg)
    }
  }

  # @ 1.1.8, following is unreachable because correct.conc requires
  # output of correct.time which ALWAYS makes column lambda_z, even if NA.
  # if (
  #   !identical(NA,th) &
  #   !("lambda_z" %in% names(x))
  # ) {
  #   x %<>% left_join(
  #     th %>% select(-no.points,-intercept,-r.squared,-adj.r.squared,-thalf),
  #     by = by
  #   )
  # }

  x %<>%
    mutate(
      ptime = x[[nomtimevar]], # nominal time                            (internal)
      crule.nr = "", # correction rule number
      crule.txt = "", # explanation of concentration substitution
      applies.to.conc = "" #,                     # lists all AUCS to which the concentration correction rule applies
      #lambda_z=ifelse("lambda_z"%in%names(.),lambda_z,NA)
    )

  # create lead and lag variables for each AUC

  #TAU
  if (!is.na(tau)) {
    x %<>%
      lag_lead(
        nomtimevar1 = "ptime",
        depvar1 = "conc.tau",
        timevar1 = "time.tau",
        lagc = "lag.ctau",
        lagt = "lag.ttau",
        leadc = "lead.ctau",
        leadt = "lead.ttau"
      )
  }
  #PARTIAL
  if (!is.na(tstart) & !is.na(tend)) {
    x %<>%
      lag_lead(
        nomtimevar1 = "ptime",
        depvar1 = "conc.part",
        timevar1 = "time.part",
        lagc = "lag.cpart",
        lagt = "lag.tpart",
        leadc = "lead.cpart",
        leadt = "lead.tpart"
      )
  }
  #TEVAL
  if (!is.na(teval)) {
    x %<>%
      lag_lead(
        nomtimevar1 = "ptime",
        depvar1 = "conc.teval",
        timevar1 = "time.teval",
        lagc = "lag.cteval",
        lagt = "lag.tteval",
        leadc = "lead.cteval",
        leadt = "lead.tteval"
      )
  }

  # Start concentration substitutions

  # Create some variables

  x %<>%
    mutate(
      t0val = ifelse(
        is.element(0, ptime),
        conc.tau[ptime == 0],
        NA
      ),
      tauval = ifelse(
        is.element(tau, ptime),
        conc.tau[ptime == tau],
        NA
      )
    )

  #TAU
  if (!is.na(tau)) {
    x %<>%
      mutate_cond(
        ptime == tau &
          is.na(conc.tau) &
          !is.na(t0val) &
          tolower(reg) == "md" &
          tolower(ss) == "y" &
          tolower(route) != "ivb",
        conc.tau = t0val, # take value at t = 0
        time.tau = tau,
        tau.flag = 1,
        crule.nr = "MDC-1",
        crule.txt = paste(
          "Missing concentration at t=",
          tau,
          "
        corrected to ",
          round(conc.tau, 3),
          " by imputation with existing concentration at t=0",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TAU ")
      )

    x %<>%
      mutate_cond(
        ptime == tau &
          is.na(conc.tau) &
          !is.na(lag.ctau) &
          !is.na(lead.ctau),
        conc.tau = interpol(
          c1 = lag.ctau,
          c2 = lead.ctau,
          t1 = tau,
          t2 = lag.ttau,
          t3 = lead.ttau,
          method = method
        ), #interpolate
        #c1+((c2-c1)*(t1-t2)/(t3-t2))
        #conc.tau=lag.ctau +
        #         (lead.ctau-lag.ctau)*(tau-lag.ttau)/(lead.ttau-lag.ttau),  # interpolate
        time.tau = tau,
        tau.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-2", "MDC-2"),
        crule.txt = paste(
          "Missing concentration at t=",
          tau,
          " corrected to ",
          round(conc.tau, 3),
          " by interpolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TAU ")
      )
    x %<>%
      mutate_cond(
        ptime == tau &
          is.na(conc.tau) &
          !is.na(lag.ctau) &
          is.na(lead.ctau) &
          !is.na(lambda_z),
        conc.tau = lag.ctau * exp(-1 * lambda_z * (tau - lag.ttau)), # extrapolate
        time.tau = tau,
        tau.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-3", "MDC-3"),
        crule.txt = paste(
          "Missing concentration at t=",
          tau,
          " corrected to ",
          round(conc.tau, 3),
          " by extrapolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TAU ")
      )
    x %<>% mutate(tauval = conc.tau[ptime == tau])
  }

  #TSTART and TEND
  #TSTART
  if (!is.na(tstart) & !is.na(tend)) {
    x %<>%
      mutate_cond(
        ptime == tstart &
          !is.na(tau) &
          tstart == tau &
          is.na(conc.part) &
          !is.na(t0val) &
          tolower(reg) == "md" &
          tolower(ss) == "y" &
          tolower(route) != "ivb",
        conc.part = t0val, # take value at t=0 if TSTART=TAU
        time.part = tstart,
        tstart.flag = 1,
        crule.nr = "MDC-1",
        crule.txt = paste(
          "Missing concentration at t=",
          tstart,
          " corrected to ",
          round(conc.part, 3),
          " by imputation with existing concentration at t=0",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TSTART ")
      )
    x %<>%
      mutate_cond(
        ptime == tstart &
          is.na(conc.part) &
          !is.na(lag.cpart) &
          !is.na(lead.cpart),
        conc.part = interpol(
          c1 = lag.cpart,
          c2 = lead.cpart,
          t1 = tstart,
          t2 = lag.tpart,
          t3 = lead.tpart,
          method = method
        ), #interpolate
        #c1+((c2-c1)*(t1-t2)/(t3-t2))
        #conc.part=lag.cpart +                                         # interpolate
        #          (lead.cpart-lag.cpart)*(tstart-lag.tpart)/(lead.tpart-lag.tpart),
        time.part = tstart,
        tstart.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-2", "MDC-2"),
        crule.txt = paste(
          "Missing concentration at t=",
          tstart,
          " corrected to ",
          round(conc.part, 3),
          " by interpolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TSTART ")
      )
    x %<>%
      mutate_cond(
        ptime == tstart &
          is.na(conc.part) &
          !is.na(lag.cpart) &
          is.na(lead.cpart) &
          !is.na(lambda_z),
        conc.part = lag.cpart * exp(-1 * lambda_z * (tstart - lag.tpart)), # extrapolate
        time.part = tstart,
        tstart.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-3", "MDC-3"),
        crule.txt = paste(
          "Missing concentration at t=",
          tstart,
          " corrected to ",
          round(conc.part, 3),
          " by extrapolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TSTART ")
      )
  }

  #TEND
  if (!is.na(tstart) & !is.na(tend)) {
    x %<>%
      mutate_cond(
        ptime == tend &
          !is.na(tau) &
          tend == tau &
          is.na(conc.part) &
          !is.na(t0val) &
          tolower(reg) == "md" &
          tolower(ss) == "y" &
          tolower(route) != "ivb",
        conc.part = t0val, # take value at t=0 if TEND=TAU
        time.part = tend,
        tend.flag = 1,
        crule.nr = "MDC-1",
        crule.txt = paste(
          "Missing concentration at t=",
          tend,
          " corrected to ",
          round(conc.part, 3),
          " by imputation with existing concentration at t=0",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TEND ")
      )
    x %<>%
      mutate_cond(
        ptime == tend &
          is.na(conc.part) &
          !is.na(lag.cpart) &
          !is.na(lead.cpart),
        conc.part = interpol(
          c1 = lag.cpart,
          c2 = lead.cpart,
          t1 = tend,
          t2 = lag.tpart,
          t3 = lead.tpart,
          method = method
        ), #interpolate
        #c1+((c2-c1)*(t1-t2)/(t3-t2))
        #conc.part=lag.cpart +                                    # interpolate
        #          (lead.cpart-lag.cpart)*(tend-lag.tpart)/(lead.tpart-lag.tpart),
        time.part = tend,
        tend.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-2", "MDC-2"),
        crule.txt = paste(
          "Missing concentration at t=",
          tend,
          " corrected to ",
          round(conc.part, 3),
          " by interpolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TEND ")
      )
    x %<>%
      mutate_cond(
        ptime == tend &
          is.na(conc.part) &
          !is.na(lag.cpart) &
          is.na(lead.cpart) &
          !is.na(lambda_z),
        conc.part = lag.cpart * exp(-1 * lambda_z * (tend - lag.tpart)), # extrapolate
        time.part = tend,
        tend.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-3", "MDC-3"),
        crule.txt = paste(
          "Missing concentration at t=",
          tend,
          " corrected to ",
          round(conc.part, 3),
          " by extrapolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TEND ")
      )
  }
  #TEVAL
  if (!is.na(teval)) {
    x %<>%
      mutate_cond(
        ptime == teval &
          !is.na(tau) &
          teval == tau &
          is.na(conc.teval) &
          !is.na(t0val) &
          tolower(reg) == "md" &
          tolower(ss) == "y" &
          tolower(route) != "ivb",
        conc.teval = t0val, # take value at t=0 if TEVAL=TAU
        time.teval = teval,
        teval.flag = 1,
        crule.nr = "MDC-1",
        crule.txt = paste(
          "Missing concentration at t=",
          teval,
          " corrected to ",
          round(conc.part, 3),
          " by imputation with existing concentration at t=0",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TEVAL ")
      )
    x %<>%
      mutate_cond(
        ptime == teval &
          is.na(conc.teval) &
          !is.na(lag.cteval) &
          !is.na(lead.cteval),
        conc.teval = interpol(
          c1 = lag.cteval,
          c2 = lead.cteval,
          t1 = teval,
          t2 = lag.tteval,
          t3 = lead.tteval,
          method = method
        ), #interpolate
        #c1+((c2-c1)*(t1-t2)/(t3-t2))
        #conc.teval=lag.cteval +                             # interpolate
        #(lead.cteval-lag.cteval)*(teval-lag.tteval)/(lead.tteval-lag.tteval),
        time.teval = teval,
        teval.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-2", "MDC-2"),
        crule.txt = paste(
          "Missing concentration at t=",
          teval,
          " corrected to ",
          round(conc.teval, 3),
          " by interpolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TEVAL ")
      )
    x %<>%
      mutate_cond(
        ptime == teval &
          is.na(conc.teval) &
          !is.na(lag.cteval) &
          is.na(lead.cteval) &
          !is.na(lambda_z),
        conc.teval = lag.cteval * exp(-1 * lambda_z * (teval - lag.tteval)), # extrapolate
        time.teval = teval,
        teval.flag = 1,
        crule.nr = ifelse(tolower(reg) == "sd", "SDC-3", "MDC-3"),
        crule.txt = paste(
          "Missing concentration at t=",
          teval,
          " corrected to ",
          round(conc.teval, 3),
          " by extrapolation",
          sep = ""
        ),
        applies.to.conc = paste(applies.to.conc, "TEVAL ")
      )
  }

  # correct t=0 conc for all aucs where t=0 is needed (for PARTIAL this is NOT needed as it does not start at t=0)

  # back-extrapolate t=0 concentration for all AUCs if route is IVB
  x %<>%
    mutate(
      back_extrap = 0,
      firstmeasc = conc.lastall[which(conc.lastall > 0)][1],
      firstmeast = ptime[which(conc.lastall > 0)][1]
    )

  # if there are LOQs between t=0 and first measurable conc, set these equal to first measurable conc
  x %<>%
    mutate_cond(
      tolower(route) == "ivb" &
        ptime > 0 &
        ptime < firstmeast &
        conc.lastall == 0,
      conc.lastall = firstmeasc,
      conc.tau = firstmeasc,
      conc.teval = firstmeasc
    )
  x %<>%
    mutate(
      lc1 = conc.lastall[which(conc.lastall > 0 & ptime > 0)][1],
      lc2 = conc.lastall[which(conc.lastall > 0 & ptime > 0)][2],
      lt1 = ptime[which(conc.lastall > 0 & ptime > 0)][1],
      lt2 = ptime[which(conc.lastall > 0 & ptime > 0)][2],
    ) %>%
    mutate_cond(
      tolower(route) == "ivb" &
        ptime == 0,
      back_extrap = ifelse(
        !is.na(lc1) &
          !is.na(lc2) &
          lc1 > 0 &
          lc2 > 0 &
          lc1 > lc2,
        1,
        0
      ),
      conc.lastall = ifelse(
        back_extrap == 1,
        exp(log(lc1) + (0 - lt1) / (lt2 - lt1) * (log(lc2) - log(lc1))),
        firstmeasc
      ),
      t0.flag = 1,
      crule.nr = ifelse(tolower(reg) == "sd", "SDC-4", "MDC-4"),
      crule.txt = paste(
        "Concentration at t=0 back-extrapolated to ",
        round(conc.lastall, 3),
        " (IV bolus)",
        sep = ""
      ),
      time.tau = 0,
      time.teval = 0,
      conc.tau = conc.lastall,
      conc.teval = conc.lastall
    )

  # ALL and LAST
  x %<>%
    mutate_cond(
      ptime == 0 &
        (is.na(conc.lastall) | conc.lastall > 0) &
        tolower(reg) == "sd" &
        tolower(route) != "ivb",
      conc.lastall = 0, # set conc to 0
      time.lastall = 0,
      t0.flag = 1,
      crule.nr = "SDC-1",
      crule.txt = paste(
        "Missing or measurable concentration at (SD) PREDOSE set to 0",
        sep = ""
      )
    )
  x %<>%
    mutate_cond(
      ptime == 0 &
        is.na(conc.lastall) &
        !is.na(tauval) &
        tolower(reg) == "md" &
        tolower(ss) == "y" &
        tolower(route) != "ivb",
      conc.lastall = tauval, # take value at t=tau
      time.lastall = 0,
      t0.flag = 1,
      crule.nr = "MDC-1",
      crule.txt = paste(
        "Missing concentration at PREDOSE corrected to ",
        round(conc.tau, 3),
        " by imputation with existing concentration at t=TAU",
        sep = ""
      ),
      applies.to.conc = paste("PREDOSE")
    )

  #TAU
  x %<>%
    mutate_cond(
      ptime == 0 &
        (is.na(conc.tau) | conc.tau > 0) &
        tolower(reg) == "sd" &
        tolower(route) != "ivb",
      conc.tau = 0, # set conc to 0
      time.tau = 0,
      t0.flag = 1,
      crule.nr = "SDC-1",
      crule.txt = paste(
        "Missing or measurable concentration at (SD) PREDOSE set to 0",
        sep = ""
      ),
      applies.to.conc = paste("PREDOSE")
    )
  x %<>%
    mutate_cond(
      ptime == 0 &
        is.na(conc.tau) &
        !is.na(tauval) &
        tolower(reg) == "md" &
        tolower(ss) == "y" &
        tolower(route) != "ivb",
      conc.tau = tauval, # take value at t=tau
      time.tau = 0,
      t0.flag = 1,
      crule.nr = "MDC-1",
      crule.txt = paste(
        "Missing concentration at PREDOSE corrected to ",
        round(conc.tau, 3),
        " by imputation with existing concentration at t=TAU",
        sep = ""
      ),
      applies.to.conc = paste("PREDOSE")
    )

  #TEVAL
  x %<>%
    mutate_cond(
      ptime == 0 &
        (is.na(conc.teval) | conc.teval > 0) &
        tolower(reg) == "sd" &
        tolower(route) != "ivb",
      conc.teval = 0, # set conc to 0
      time.teval = 0,
      t0.flag = 1,
      crule.nr = "SDC-1",
      crule.txt = paste(
        "Missing or measurable concentration at (SD) PREDOSE set to 0",
        sep = ""
      ),
      applies.to.conc = paste("PREDOSE")
    )
  x %<>%
    mutate_cond(
      ptime == 0 &
        is.na(conc.teval) &
        !is.na(tauval) &
        tolower(reg) == "md" &
        tolower(ss) == "y" &
        tolower(route) != "ivb",
      conc.teval = tauval, # take value at t=tau
      time.teval = 0,
      t0.flag = 1,
      crule.nr = "MDC-1",
      crule.txt = paste(
        "Missing concentration at PREDOSE corrected to ",
        round(conc.teval, 3),
        " by imputation with existing concentration at t=TAU",
        sep = ""
      ),
      applies.to.conc = paste("PREDOSE")
    )
  x
}
