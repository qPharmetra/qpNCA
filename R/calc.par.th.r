globalVariables('par')
#' Calculate Lambda_z Parameters
#'
#' Calculates PK parameters that need lambda_z.
#'
#' @param x result parameter dataset from \code{\link{calc.par}}
#' @param by column names in x indicating grouping variables
#' @param th result dataset from \code{\link{est.thalf}}
#' @param covariates covariates dataset (containing at least dose for CL calculation); defaults to unique combinations of \code{by} and \code{dose} evaluated on \code{x}; can be character name of csv file or local object
#' @param dose variable containing the dose amount; default 'dose' set to 1 if not in \code{names(x)}
#' @param factor conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000); x$factor overrides
#' @param reg regimen, "sd" or "md"; x$reg overrides
#' @param ss is steady state reached (y/n); x$ss overrides
#' @param route of drug administration ("EV","IVB","IVI"); x$route overrides
#' @importFrom dplyr left_join
#'
#'
#' @return A dataset containing all parameters calculated in \code{\link{est.thalf}} and \code{\link{calc.par}} \cr
#' with estimates for the following parameters added, one observation per subject:
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
#' example(calc.par) # creates par
#' # notice x includes (optional) loqrule, includeCmax, reg, method, route, ss
#' covs <- Theoph %>%
#'   select(subject = Subject, Wt, dose = Dose) %>%
#'   unique %>%
#'   mutate(dose = dose * Wt, subject=as.numeric(as.character(subject))) # see ?Theoph
#' y <- x %>% select(subject, reg, ss, loqrule) %>% unique
#' y %<>% mutate(factor = 1)
#' par %<>% left_join(y, by = 'subject')
#' par %<>% calc.par.th(by = 'subject', th = th, covariates = covs)
#' par %<>% left_join(ctmax, ., by = 'subject')
#' par %>% head
#' par %>% data.frame %>% head(2)
#' }
calc.par.th <- function(
  x,
  by=character(0),
  th=th,
  covariates=NA,
  dose="dose",
  factor=1,
  reg="SD",
  ss="N",
  route="EV"
  ){
  stopifnot(
    is.character(dose),
    length(dose) == 1
  )
  if(!dose %in% names(x)){
    x[[dose]] <- 1
  }
  if(identical(NA, covariates)){
    covariates <- unique(x[,c(by, dose),drop=FALSE])
  }
  # At this point, utility of having assigned dose on x is ended.
  # User must have supply dose in covariates using column named in arg 'dose'.
  # To avoid merge conflicts of covariates and x, we squelch x[[dose]].
  x[[dose]] <- NULL

  for(arg in c('reg','ss','factor','route')){
    if(arg %in% names(x)){
      if(!eval(substitute(missing(arg)))){
        warning(arg,' supplied as column overrides like-named argument')
      }
      assign(arg,x[[arg]])
      x[[arg]] <- NULL
    }
    # if(length(get(arg)) > 1) {
    #   warning(arg, ' has length > 1; only first value will be used')
    #   assign(arg, get(arg)[[1]])
    # }
  }

  if(!missing(covariates)){
    if(is.character(covariates)){
      if(file.exists(covariates)){
        covariates <- read.csv(covariates)
      }else{
        covariates=get(covariates) # to convert the covariates string to the real data frame
      }
    }
  }
  result=left_join(x,th,by=by) %>%
  left_join(covariates,by=intersect(names(x), names(covariates)))
  result$dosevar <- result[[dose]]
  result %<>%
  mutate(
    factor=ifelse(is.na(factor),1,factor),
    reg=tolower(reg),
    ss=tolower(ss),
    route=tolower(route),
    clast.pred=exp(log(intercept)-lambda_z*tlast), # intercept was already exponentiated in calc_thalf
    aucinf.obs=auclast+clast.obs/lambda_z,
    aucinf.pred=auclast+clast.pred/lambda_z,
    aumcinf.obs=aumclast+tlast*clast.obs/lambda_z + clast.obs/lambda_z^2,
    aumcinf.pred=aumclast+tlast*clast.pred/lambda_z + clast.pred/lambda_z^2,
    cl.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/aucinf.obs, dosevar*factor/auctau),
    cl.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/aucinf.pred, dosevar*factor/auctau),
    mrtinf.obs= ifelse(
      tolower(ss)=="n",
      aumcinf.obs/aucinf.obs,
      (aumctau + tau*(aucinf.obs-auctau))/auctau
    ),
    mrtinf.pred= ifelse(
      tolower(ss)=="n",
      aumcinf.pred/aucinf.pred,
      (aumctau + tau*(aucinf.pred-auctau))/auctau
    ),
    vz.f.obs=  ifelse(
      tolower(ss)=="n",
      dosevar*factor/(lambda_z*aucinf.obs),
      dosevar*factor/(lambda_z*auctau)
    ),
    vz.f.pred= ifelse(
      tolower(ss)=="n",
      dosevar*factor/(lambda_z*aucinf.pred),
      dosevar*factor/(lambda_z*auctau)
    ),
    #   vss.obs= mrtinf.obs*cl.f.obs,
    #   vss.pred= mrtinf.pred*cl.f.pred,
    vss.obs= ifelse(route=="ivb"|route=="ivi",mrtinf.obs*cl.f.obs,NA),
    vss.pred= ifelse(route=="ivb"|route=="ivi",mrtinf.pred*cl.f.pred,NA),
    pctextr.obs=(clast.obs/lambda_z)/aucinf.obs*100,
    pctextr.pred=(clast.pred/lambda_z)/aucinf.pred*100,
    pctback.obs=area.back.extr/aucinf.obs*100,
    pctback.pred=area.back.extr/aucinf.pred*100
  ) %>%
  select(-dosevar)

  # rename CL and V after IVB/IVI

  qpiv <- tolower(route)=="ivb"|tolower(route)=="ivi"
  result %<>% mutate( cl.obs =    ifelse( qpiv, cl.f.obs,  NA))
  result %<>% mutate( cl.pred =   ifelse( qpiv, cl.f.pred, NA))
  result %<>% mutate( vz.obs =    ifelse( qpiv, vz.f.obs,  NA))
  result %<>% mutate( vz.pred =   ifelse( qpiv, vz.f.pred, NA))
  result %<>% mutate( cl.f.obs =  ifelse(!qpiv, cl.f.obs,  NA))
  result %<>% mutate( cl.f.pred = ifelse(!qpiv, cl.f.pred, NA))
  result %<>% mutate( vz.f.obs =  ifelse(!qpiv, vz.f.obs,  NA))
  result %<>% mutate( vz.f.pred = ifelse(!qpiv, vz.f.pred, NA))

# rename CL at steady state
  ssiv <- tolower(ss) == 'y' &  qpiv
  ssev <- tolower(ss) == 'y' & !qpiv

  result %<>% mutate(cl.f.ss = NA, cl.ss = NA)

  result %<>% mutate(cl.ss     = ifelse( ssiv, cl.obs,   cl.ss   ))
  result %<>% mutate(cl.obs    = ifelse(!ssiv, cl.obs,   NA      ))
  result %<>% mutate(cl.pred   = ifelse(!ssiv, cl.pred,  NA      ))

  result %<>% mutate(cl.f.ss   = ifelse( ssev, cl.f.obs, cl.f.ss ))
  result %<>% mutate(cl.f.obs  = ifelse(!ssev, cl.f.obs, NA      ))
  result %<>% mutate(cl.f.pred = ifelse(!ssev, cl.f.pred,  NA    ))

  # if (tolower(ss)=="y") {
  #   if (tolower(route)=="ivb"|tolower(route)=="ivi") {
  #     result = result %>% mutate( cl.ss   = cl.obs,   cl.obs   = NA, cl.pred   = NA, cl.f.ss = NA )
  #   }else{
  #     result = result %>% mutate( cl.f.ss = cl.f.obs, cl.f.obs = NA, cl.f.pred = NA, cl.ss = NA )
  #   }
  # }else{
  #   result = result %>% mutate(cl.f.ss=NA, cl.ss=NA)
  # }

  return(result)
}
