globalVariables('par')
#' Calculate Lambda_z Parameters
#'
#' Calculates PK parameters, including those that need lambda_z.
#'
#' @param x output from \code{\link{correct.conc}}
#' @param by character: column names in x indicating grouping variables; default is as.character(dplyr::groups(x))
#' @param tau dosing interval (for multiple dosing); NA (default) for if single dose; x$tau overrides
#' @param tstart start time of partial AUC (start>0); NA (default) if not requested; x$tstart overrides
#' @param tend end time of partial AUC; NA (default) if not requested; x$tend overrides
#' @param teval user selected AUC interval; NA (default) if not requested; x$teval overrides
#' @param route route of drug administration ("EV","IVB","IVI"); x$route overrides
#' @param method method for trapezoidal rule
#' @param th lamdba_z information for each curve; default is est.thalf(x, by = by, timevar = timevar, depvar = depvar )
#' @param covariates covariates dataset (containing at least dose for CL calculation); defaults to unique combinations of \code{by} and \code{dose} evaluated on \code{x}; can be character name of csv file or local object
#' @param dose variable containing the dose amount; default 'dose' set to 1 if not in \code{names(x)}
#' @param infdur variable containing the infusion duration; default 'infdur' set to 1 if not in \code{names(x)}
#' @param factor conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000); x$factor overrides
#' @param reg regimen, "sd" or "md"; x$reg overrides
#' @param ss is steady state reached (y/n); x$ss overrides
#' @param timevar variable name containing the actual sampling time after dose
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @importFrom dplyr left_join
#'
#'
#' @return A dataset containing all parameters calculated in \code{\link{est.thalf}} and \code{\link{calc.par}} \cr
#' with estimates for the following parameters added, one observation per profile:
#'
#' **Parameter** | **Description**
#' ------------ | -----------
#' clast.pred   | predicted concentration at tlast
#' aucinf.obs   | aucinf based on observed concentration at tlast
#' aucinf.pred  | aucinf based on predicted concentration at tlast
#' aumcinf.obs  | area under the first moment curve extrapolated to infinity, based on observed concentration at tlast
#' aumcinf.pred | area under the first moment curve extrapolated to infinity, based on predicted concentration at tlast
#' cl.obs, cl.f.obs     | clearance based on aucinf.obs, at steady state based on auctau
#' cl.pred, cl.f.pred    | clearance based on aucinf.pred
#' cl.ss, cl.f.ss     | clearance at steady state, based on auctau
#' mrt.obs      | mean residence time based on aumcinf.obs and aucinf.obs
#' mrt.pred     | mean residence time based on aumcinf.pred and aucinf.pred
#' vz.obs, vz.f.obs     | distribution volume based on cl.f.obs, at steady state based on auctau
#' vz.pred, vz.f.pred    | distribution based on cl.pred/cl.f.pred
#' vss.obs      | steady-state volume based on cl.obs and mrt.obs
#' vss.pred     | steady-state volume based on cl.pred and mrt.pred
#' pctextr.pred	| percentage of AUC extrapolated to infinity, based on aucinf.pred
#' pctextr.obs	| percentage of AUC extrapolated to infinity, based on aucinf.obs
#' pctback.pred	| percentage of AUC extrapolated back to 0, based on aucinf.pred
#' pctback.obs	| percentage of AUC extrapolated back to 0, based on aucinf.obs
#'
#' Note: ctmax must be merged separately as those were calculated from uncorrected data.
#' @export
#' @examples
#' \donttest{
#' library(magrittr)
#' library(dplyr)
#' data(ncx)
#' x <- ncx
#' x %<>% group_by(subject)
#' x %<>% correct.loq
#' cm <- calc.ctmax(x)              # using uncorrected data
#' th <- est.thalf(x)               # possibly expensive!
#' x %<>% correct.time              # or pass (th = th)
#' x %<>% correct.conc
#' x %<>% calc.par.th               # or pass (th = th)
#' x %<>% left_join(cm)
#' x %>% data.frame %>% head(2)
#' }
calc.par.th <- function(
  x,
  by = NULL,
  tau = NA,
  tstart = NA,
  tend = NA,
  teval = NA,
  route = "EV",
  method = 1,
  th = NULL,
  covariates = NA,
  dose = "dose",
  infdur = "infdur",
  factor = 1,
  reg = "SD",
  ss = "N",
  timevar = "time",
  depvar = "dv"
) {
  if (is.null(by)) by <- as.character(groups(x))

  targs <- list(x = x, by = by)
  if (!missing(timevar)) targs <- c(targs, list(timevar = timevar))
  if (!missing(depvar)) targs <- c(targs, list(depvar = depvar))
  if (is.null(th)) th <- do.call(est.thalf, targs)

  pargs <- list(x = x, by = by)
  if (!missing(tau)) pargs <- c(pargs, list(tau = tau))
  if (!missing(tstart)) pargs <- c(pargs, list(tstart = tstart))
  if (!missing(tend)) pargs <- c(pargs, list(tend = tend))
  if (!missing(teval)) pargs <- c(pargs, list(tval = tval))
  if (!missing(route)) pargs <- c(pargs, list(route = route))
  if (!missing(method)) pargs <- c(pargs, list(method = method))

  # direct call causes nothing to be missing, with spurious warnings
  #   x %<>% calc.par(
  #   by = by,
  #   tau = tau,
  #   tstart = tstart,
  #   tend = tend,
  #   teval = teval,
  #   route = route,
  #   method = method
  # )
  p <- do.call(calc.par, pargs)

  # retain if available ...
  nms <- c('reg', 'ss', 'factor', 'route', 'loqrule', 'tauval')
  suppressMessages(
    x <- left_join(p, unique(select(x, !!by, !!intersect(nms, names(x)))))
  )
  stopifnot(
    is.character(dose),
    length(dose) == 1
  )
  if (!dose %in% names(x)) {
    x[[dose]] <- 1
  }
  stopifnot(
    is.character(infdur),
    length(infdur) == 1
  )
  if (!infdur %in% names(x)) {
    x[[infdur]] <- 1
  }

  if (identical(NA, covariates)) {
    covariates <- unique(x[, c(by, dose, infdur), drop = FALSE])
  }
  # At this point, utility of having assigned dose on x is ended.
  # User must have supplied dose in covariates using column named in arg 'dose'.
  # To avoid merge conflicts of covariates and x, we squelch x[[dose]].
  x[[dose]] <- NULL
  x[[infdur]] <- NULL

  for (arg in c('reg', 'ss', 'factor', 'route')) {
    if (arg %in% names(x)) {
      if (!eval(substitute(missing(arg)))) {
        warning(arg, ' supplied as column overrides like-named argument')
      }
      # assign(arg,x[[arg]])
      # x[[arg]] <- NULL
    } else {
      x %<>% mutate(!!arg := get(arg)) # deliberate switch to data priority vs global env
    }

    # if(length(get(arg)) > 1) {
    #   warning(arg, ' has length > 1; only first value will be used')
    #   assign(arg, get(arg)[[1]])
    # }
  }

  if (!missing(covariates)) {
    if (is.character(covariates)) {
      if (file.exists(covariates)) {
        covariates <- read.csv(covariates)
      } else {
        covariates = get(covariates) # to convert the covariates string to the real data frame
      }
    }
  }
  x %<>%
    left_join(th, by = by) %>%
    left_join(covariates, by = intersect(names(x), names(covariates)))
  x$dosevar <- x[[dose]]
  x$infdurvar <- x[[infdur]]
  x %<>%
    mutate(
      factor = ifelse(is.na(factor), 1, factor),
      reg = tolower(reg),
      ss = tolower(ss)
    )

  x %<>% mutate(route = tolower(route))
  x %<>%
    mutate(
      clast.pred = exp(log(intercept) - lambda_z * tlast), # intercept was already exponentiated in calc_thalf
      aucinf.obs = auclast + clast.obs / lambda_z,
      aucinf.pred = auclast + clast.pred / lambda_z,
      aumcinf.obs = aumclast +
        tlast * clast.obs / lambda_z +
        clast.obs / lambda_z^2,
      aumcinf.pred = aumclast +
        tlast * clast.pred / lambda_z +
        clast.pred / lambda_z^2,
      cl.f.obs = ifelse(
        tolower(ss) == "n",
        dosevar * factor / aucinf.obs,
        dosevar * factor / auctau
      ),
      cl.f.pred = ifelse(
        tolower(ss) == "n",
        dosevar * factor / aucinf.pred,
        dosevar * factor / auctau
      ),
      inffactor = ifelse(tolower(route) == "ivi", infdurvar / 2, 0),
      mrtinf.obs = ifelse(
        tolower(ss) == "n",
        aumcinf.obs / aucinf.obs - inffactor,
        (aumctau + tau * (aucinf.obs - auctau)) / auctau - inffactor
      ),
      mrtinf.pred = ifelse(
        tolower(ss) == "n",
        aumcinf.pred / aucinf.pred - inffactor,
        (aumctau + tau * (aucinf.pred - auctau)) / auctau - inffactor
      ),
      mrtall = mrtall - inffactor, # mrtall is calculated in calc.par
      mrtlast = mrtlast - inffactor, # mrtlast is calculated in calc.par
      vz.f.obs = ifelse(
        tolower(ss) == "n",
        dosevar * factor / (lambda_z * aucinf.obs),
        dosevar * factor / (lambda_z * auctau)
      ),
      vz.f.pred = ifelse(
        tolower(ss) == "n",
        dosevar * factor / (lambda_z * aucinf.pred),
        dosevar * factor / (lambda_z * auctau)
      ),
      #   vss.obs= mrtinf.obs*cl.f.obs,
      #   vss.pred= mrtinf.pred*cl.f.pred
    )
  x %<>%
    mutate(
      vss.obs = ifelse(
        route == "ivb" | route == "ivi",
        mrtinf.obs * cl.f.obs,
        NA
      )
    ) #cannot be calculated after EV
  x %<>%
    mutate(
      vss.pred = ifelse(
        route == "ivb" | route == "ivi",
        mrtinf.pred * cl.f.pred,
        NA
      )
    ) #cannot be calculated after EV
  x %<>%
    mutate(
      pctextr.obs = (clast.obs / lambda_z) / aucinf.obs * 100,
      pctextr.pred = (clast.pred / lambda_z) / aucinf.pred * 100,
      pctback.obs = area.back.extr / aucinf.obs * 100,
      pctback.pred = area.back.extr / aucinf.pred * 100
    )
  x %<>% select(-dosevar)

  # rename CL and V after IVB/IVI

  x %<>% mutate(qpiv = tolower(route) == "ivb" | tolower(route) == "ivi")
  x %<>% mutate(cl.obs = ifelse(qpiv, cl.f.obs, NA))
  x %<>% mutate(cl.pred = ifelse(qpiv, cl.f.pred, NA))
  x %<>% mutate(vz.obs = ifelse(qpiv, vz.f.obs, NA))
  x %<>% mutate(vz.pred = ifelse(qpiv, vz.f.pred, NA))
  x %<>% mutate(cl.f.obs = ifelse(!qpiv, cl.f.obs, NA))
  x %<>% mutate(cl.f.pred = ifelse(!qpiv, cl.f.pred, NA))
  x %<>% mutate(vz.f.obs = ifelse(!qpiv, vz.f.obs, NA))
  x %<>% mutate(vz.f.pred = ifelse(!qpiv, vz.f.pred, NA))

  # rename CL at steady state
  x %<>% mutate(ssiv = tolower(ss) == 'y' & qpiv)
  x %<>% mutate(ssev = tolower(ss) == 'y' & !qpiv)

  x %<>% mutate(cl.f.ss = NA, cl.ss = NA)

  x %<>% mutate(cl.ss = ifelse(ssiv, cl.obs, cl.ss)) #cl.obs and cl.pred are equal at steady state (AUCtau used)
  x %<>% mutate(cl.obs = ifelse(!ssiv, cl.obs, NA))
  x %<>% mutate(cl.pred = ifelse(!ssiv, cl.pred, NA))

  x %<>% mutate(cl.f.ss = ifelse(ssev, cl.f.obs, cl.f.ss)) #cl.f.obs and cl.f.pred are equal at steady state (AUCtau used)
  x %<>% mutate(cl.f.obs = ifelse(!ssev, cl.f.obs, NA))
  x %<>% mutate(cl.f.pred = ifelse(!ssev, cl.f.pred, NA))

  # done with these:
  x %<>% select(-qpiv, -ssiv, -ssev)

  # if (tolower(ss)=="y") {
  #   if (tolower(route)=="ivb"|tolower(route)=="ivi") {
  #     x = x %>% mutate( cl.ss   = cl.obs,   cl.obs   = NA, cl.pred   = NA, cl.f.ss = NA )
  #   }else{
  #     x = x %>% mutate( cl.f.ss = cl.f.obs, cl.f.obs = NA, cl.f.pred = NA, cl.ss = NA )
  #   }
  # }else{
  #   x = x %>% mutate(cl.f.ss=NA, cl.ss=NA)
  # }
  # x %<>% select((!!names(.)[[1]]):loqrule, dose, everything()) # to preserve traditional column order
  x
}
