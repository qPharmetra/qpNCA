#' A wrapper for the individual qpNCA functions
#' This wrapper function consecutively executes the following qpNCA functions:
#' 1. correct.loq \cr
#' 2. est.thalf\cr
#' 2a. plot.reg\cr
#' 3. calc.ctmax\cr
#' 4. correct.time\cr
#' 5. correct.conc\cr
#' 6. tab.corr\cr
#' 7. calc.par\cr
#' 8. calc.par.th\cr
#'
#' USAGE:
#'
#' qpNCA(x, by=c("subject"), nomtimevar="ntad", timevar="time",depvar="dv",
#'       includeCmax="Y",excl="excl",savedir=NA,
#'       bloqvar="bloq",loqvar="loq",loqrule=1,
#'       tau=NA,tstart=NA,tend=NA,teval=NA,cov=NA,dose=NA,factor=NA,reg="SD",ss="N",route="PO",method=1) {
#'
#' @param x input dataset name (if called within dplyr: .)
#' @param by by-variable(s), e.g. c("subject","day")
#' @param nomtimevar variable name containing the nominal sampling time
#' @param timevar variable name containing the sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param X-axis label (default: "timevar") in plots
#' @param Y-axis label (default: "depvar") in plots
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param loqvar variable name containing the LOQ value
#' @param loqrule rule number to be applied to the LOQ values in the curve
#' @param includeCmax include results of regression including Cmax in selection? (y/n)
#' @param excl variable name containing information about points to be excluded (these should have <excl>=1)
#' @param savedir folder where regression plots (.PNG) will be saved; leave empty for standard output
#' @param tau dosing interval (for multiple dosing), if single dose, leave empty
#' @param tstart start time of partial AUC (start>0), if not requested, leave empty
#' @param tend end time of partial AUC, if not requested, leave empty
#' @param teval user selected AUC interval, if not requested, leave empty
#' @param cov covariates dataset
#' @param dose variable containing the dose amount
#' @param factor onversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000)
#' @param reg regimen, "sd" or "md"
#' @param ss is steady state reached (y/n)
#' @param route route of drug administration ("po","iv")
#' @param method method for trapezoidal rule:
#'              1: linear up - linear down
#'              2: linear up - logarithmic down
#'              3: linear before first Tmax, logarithmic after first Tmax
#'
#' @return
#' covariates : contains covariates selected with the 'cov' argument
#' half_life  : contains linear regression parameters
#' ctmax      : contains cmax and tmax estimated from uncorrected data
#' ct_corr    : contains the time and concentration corrected dataset
#' corrections: contains descriptions of the corrections applied
#' pkpar      : contains all estimated PK parameters
#'
#' @export
qpNCA <- function(x, by=c("subject"), nomtimevar="ntad", timevar="time",depvar="dv",
                  bloqvar="bloq",loqvar="loq",loqrule=1,
                  includeCmax="Y",excl="excl",savedir=NA,timelab="timevar",deplab="depvar",
                  tau=NA,tstart=NA,tend=NA,teval=NA,cov=NA,dose=NA,factor=NA,reg="SD",ss="N",route="PO",method=1) {

  # 1. Apply LOQ rules

  loqed = x %>%
    group_by_(by) %>%
    do(correct.loq(.,nomtimevar=nomtimevar,timevar=timevar,depvar=depvar,bloqvar=bloqvar,loqvar=loqvar,loqrule=loqrule)) %>%
    ungroup

  # 2. estimate thalf ON UNCORRECTED DATA

  th = loqed %>%
    group_by_(by) %>%
    do(est.thalf(.,timevar=timevar,depvar=depvar,includeCmax=includeCmax,excl=excl)) %>%
    ungroup

  # 2a.

  plot.reg(loqed,by=by,th=th,bloqvar=bloqvar,timevar=timevar,depvar=depvar,excl=excl,savedir=savedir,timelab=timelab,deplab=deplab)

  # 3. find Cmax and tmax ON UNCORRECTED DATA

  ctmax = loqed %>%
    group_by_(by) %>%
    do(calc.ctmax(.,timevar=timevar,depvar=depvar)) %>%
    ungroup

  # 4. create dataset with corrected time deviations

  tc = loqed %>%
    group_by_(by) %>%
    do(correct.time(.,nomtimevar=nomtimevar,timevar=timevar,depvar=depvar,
                    tau=tau,tstart=tstart,tend=tend,teval=teval,th=th,reg=reg,method=method)) %>%
    do(correct.conc(.,nomtimevar=nomtimevar,
                    tau=tau,tstart=tstart,tend=tend,teval=teval,th=th,reg=reg,ss=ss,route=route,method=method)) %>%
    ungroup

  corrtab = tc %>% tab.corr(.,nomtimevar=nomtimevar,by=by)

  # 5. Calculate PK parameters NOT based on lambda_z ON CORRECTED DATA

  par = tc %>%
    group_by_(by) %>%
    do(calc.par(.,tau=tau,tstart=tstart,tend=tend,teval=teval,route=route,method=method)) %>%
    ungroup

  # 6.  Calculate PK parameters that need lambda_z

  par_all = calc.par.th(x=par,th=th,cov=cov,dose=dose,reg=reg,ss=ss,factor=factor)

  par_all = left_join(ctmax,par_all)

  result=list(covariates=cov,
              half_life=th,
              ctmax=ctmax,
              ct_corr=tc,
              corrections=corrtab,
              pkpar=par_all)

  return(result)

}
