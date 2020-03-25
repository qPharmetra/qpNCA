#' A wrapper for the individual qpNCA functions
#'
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
#' @param plotdir folder where regression plots (.PNG) will be saved; leave empty for standard output
#' @param tau dosing interval (for multiple dosing), if single dose, leave empty
#' @param tstart start time of partial AUC (start>0), if not requested, leave empty
#' @param tend end time of partial AUC, if not requested, leave empty
#' @param teval user selected AUC interval, if not requested, leave empty
#' @param cov covariates dataset
#' @param dose variable containing the dose amount
#' @param factor onversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000)
#' @param reg regimen, "sd" or "md"
#' @param ss is steady state reached (y/n)
#' @param route route of drug administration ("EV","IVB", "IVI")
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
#' @export
qpNCA <- function(x, by=c("subject"), nomtimevar="ntad", timevar="time",depvar="dv",
                  bloqvar="bloq",loqvar="loq",loqrule=1,
                  includeCmax="Y",exclvar=NA,plotdir=NA,pdfdir=NA,timelab="timevar",deplab="depvar",
                  tau=NA,tstart=NA,tend=NA,teval=NA,covfile=NA,dose=NA,factor=NA,reg="SD",ss="N",route="EV",method=1) {
  
  # 0. Check input
  
  cat("\n")
  cat("Checking function arguments...")
  
  check.input(x, by=by, nomtimevar=nomtimevar, timevar=timevar, depvar=depvar,bloqvar=bloqvar, loqvar=loqvar, loqrule=loqrule,
             includeCmax=includeCmax, exclvar=exclvar, plotdir=plotdir, pdfdir=pdfdir, timelab=timelab, deplab=deplab,
             tau=tau, tstart=tstart, tend=tend, teval=teval, covfile=covfile, dose=dose, factor=factor, reg=reg, ss=ss,
             route=route, method=method)
  
  # 1. Apply LOQ rules
  
  cat("Applying LOQ rules...\n")
  
  loqed = x %>%
    group_by_at(by) %>%
    do(correct.loq(.,nomtimevar=nomtimevar,timevar=timevar,depvar=depvar,bloqvar=bloqvar,loqvar=loqvar,loqrule=loqrule)) %>%
    ungroup
  
    # 2. estimate thalf ON UNCORRECTED DATA
  
  cat("Performing Thalf estimation...\n")
  
  th = loqed %>% 
    group_by_at(by) %>%
    do(est.thalf(.,timevar=timevar,depvar=depvar,includeCmax=includeCmax,excl=excl)) %>%
    ungroup

  # 2a.
  
  if (is.na(plotdir)) cat("Creating regression plots in standard output...\n")
  else cat(paste("Writing regression plots to folder",plotdir,"...\n"))
  
  plot.reg(loqed,by=by,th=th,bloqvar=bloqvar,timevar=timevar,depvar=depvar,exclvar=exclvar,plotdir=plotdir,timelab=timelab,deplab=deplab)
  
  cat("\n")
  
  # 3. find Cmax and tmax ON UNCORRECTED DATA
  
  cat("Calculating Cmax/Tmax...\n")
  
  ctmax = loqed %>% 
    group_by_at(by) %>%
    do(calc.ctmax(.,timevar=timevar,depvar=depvar)) %>%
    ungroup

  # 4. and 5. create dataset with corrected time deviations
  
  cat("Applying time deviation corrections and missing concentration imputations...\n")
  
  tc = loqed %>%
    group_by_at(by) %>%
    do(correct.time(.,by=by,nomtimevar=nomtimevar,timevar=timevar,depvar=depvar,
                    tau=tau,tstart=tstart,tend=tend,teval=teval,th=th,reg=reg,method=method)) %>%
    do(correct.conc(.,by=by,nomtimevar=nomtimevar,
                    tau=tau,tstart=tstart,tend=tend,teval=teval,th=,reg=reg,ss=ss,route=route,method=method)) %>%
    ungroup
  
  cat("\n")
  
  # 6. create table with corrections
  
  cat("Creating correction tables...\n")
  
  corrtab = tc %>% tab.corr(.,nomtimevar=nomtimevar,by=by) 
  
  # 7. Calculate PK parameters NOT based on lambda_z ON CORRECTED DATA
  
  cat("Calculating parameters that do not need lambda_z...\n")
  
  par = tc %>% 
    group_by_at(by) %>%
    do(calc.par(.,tau=tau,tstart=tstart,tend=tend,teval=teval,route=route,method=method)) %>%
    ungroup

  # 8. Calculate PK parameters that need lambda_z
  
  cat("Calculating parameters that DO need lambda_z...\n")
  
  par_all = calc.par.th(x=par,by=by,th=th,covfile=covfile,dose=dose,reg=reg,ss=ss,factor=factor,route=route)
  
  cat("Combining all parameters...\n")
  
  par_all = left_join(ctmax,par_all,by=by)
  
  # 9. create summary PDFs
  
  if (is.na(pdfdir)) cat("No PDF summaries created\n")
  else cat(paste("Writing summary PDF documents to folder",pdfdir,"...\n"))
  
  nca.sum(par_all,corrfile=corrtab,by=by,plotdir=plotdir,pdfdir=pdfdir)
  
  cat("\nWriting results...\n")
  
  covfile=get(covfile) # to convert the cov string to the real data frame
  
  result=list(covariates=covfile,
              half_life=th,
              ct_corr=tc,
              corrections=corrtab,
              pkpar=par_all)

  cat("\nDone!\n")
  
  return(result)

}
