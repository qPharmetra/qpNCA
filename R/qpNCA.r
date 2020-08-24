#' Perform Non-compartmental Analysis
#'
#' Consecutively executes the following NCA steps:
#' \describe{
#'  \item{\code{link{correct.loq}}}{impute LOQ values}
#'  \item{\code{link{est.thalf}}}{calculate lambda-z and halflife}
#'  \item{\code{link{plot.reg}}}{plot each curve}
#'  \item{\code{link{calc.ctmax}}}{calculate Cmax and Tmax}
#'  \item{\code{link{correct.time}}}{supply records at critical times with correct concentration}
#'  \item{\code{link{correct.conc}}}{supply correct concentrations at critical times}
#'  \item{\code{link{tab.corr}}}{tabulate data alterations}
#'  \item{\code{link{calc.par}}}{calculates profile-specific summary statistics}
#'  \item{\code{link{calc.par.th}}}{calculates parameters depednent on lambda-z}
#'  }
#'
#' USAGE:
#'
#' qpNCA(x, by=c("subject"), nomtimevar="ntad", timevar="time",depvar="dv",
#'       includeCmax="Y",exclvar="exclvar",savedir=NA,
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
#' @param loqrule rule number to be applied to the LOQ values in the curve; x$loqrule overrides if provided
#' @param includeCmax include results of regression including Cmax in selection? (y/n) x$includeCmax overrides if provided
#' @param exclvar variable name containing information about points to be excluded (these should have <exclvar>=1)
#' @param plotdir folder where regression plots (.PNG) will be saved; leave empty for standard output
#' @param tau dosing interval (for multiple dosing), if single dose, leave empty; x$tau overrides if provided
#' @param tstart start time of partial AUC (start>0), if not requested, leave empty; x$tstart overrides if provided
#' @param tend end time of partial AUC, if not requested, leave empty; x$tend overrides if provided
#' @param teval user selected AUC interval, if not requested, leave empty; x$tend overrides if provided
#' @param cov covariates dataset
#' @param dose variable containing the dose amount
#' @param factor conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000); x$factor overrides if provided
#' @param reg regimen, "sd" or "md"; can be character column name in x; x$regimen overrides if provided
#' @param ss is steady state reached (y/n); x$ss overrides if provided
#' @param route route of drug administration ("EV","IVB", "IVI"); x$overrides if provided
#' @param method method for trapezoidal rule;  x$method overrides if provided
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

qpNCA <- function(
  x,
  by=c("subject"),
  nomtimevar="ntad",
  timevar="time",
  depvar="dv",
  bloqvar="bloq",
  loqvar="loq",
  loqrule=1,
  includeCmax="Y",
  exclvar=NA,
  plotdir=NA,
  pdfdir=NA,
  timelab="timevar",
  deplab="depvar",
  tau=NA,
  tstart=NA,
  tend=NA,
  teval=NA,
  covfile=NA,
  dose=NA,
  factor=NA,
  reg="SD",
  ss="N",
  route="EV",
  method=1
){

  # 0. Check input

  cat("\n")
  cat("Checking function arguments...")

  enforce <- c(
    'ss','route','reg','method','loqrule','includeCmax',
    'tau','tstart','tend','teval','factor'
  )

  enforced <- setdiff(enforce, names(x))

  for(arg in enforce){
    if(!arg %in% names(x)) x[[arg]] <- get(arg)
  }

 rm(list = enforce)

check.input(
    x, by=by, nomtimevar=nomtimevar, timevar=timevar, depvar=depvar,
    bloqvar=bloqvar, loqvar=loqvar,
    #loqrule=loqrule, includeCmax=includeCmax,
    exclvar=exclvar, plotdir=plotdir, pdfdir=pdfdir, timelab=timelab,
    deplab=deplab,
    #tau=tau, tstart=tstart, tend=tend, teval=teval,
    covfile=covfile, dose=dose#,
    #factor=factor, reg=reg, ss=ss,route=route, method=method
  )

  # 1. Apply LOQ rules

  cat("Applying LOQ rules...\n")

  loqed = x %>%
    group_by_at(by) %>%
    do(
      correct.loq(
        .,
        nomtimevar=nomtimevar,
        timevar=timevar,
        depvar=depvar,
        bloqvar=bloqvar,
        loqvar=loqvar#,
       # loqrule=loqrule
      )
    ) %>%
    ungroup

   loqed$loqrule <- x$loqrule

    # 2. estimate thalf ON UNCORRECTED DATA

  cat("Performing Thalf estimation...\n")

  th = loqed %>%
    group_by_at(by) %>%
    do(est.thalf(
      .,timevar=timevar,depvar=depvar,
      #includeCmax=includeCmax,
      exclvar=exclvar
    )) %>%
    ungroup

  # 2a.

  if(is.null(plotdir)){
    cat('Skipping regression plots')
  }else{
    if (is.na(plotdir)){
      cat("Creating regression plots in standard output...\n")
    }else{
      cat(paste("Writing regression plots to folder",plotdir,"...\n"))
    }
    plot.reg(loqed,by=by,th=th,bloqvar=bloqvar,timevar=timevar,depvar=depvar,exclvar=exclvar,plotdir=plotdir,timelab=timelab,deplab=deplab)
  }
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
    do(correct.time(
      .,by=by,nomtimevar=nomtimevar,timevar=timevar,depvar=depvar,th=th#,
      #tau=tau,tstart=tstart,tend=tend,teval=teval,reg=reg,method=method
    ))

  for(arg in c(
    'tau','tstart','tend','teval',
    'reg','ss','route','method'
  )){
    if(arg %in% names(loqed))if(!arg %in% names(tc)){
      tc[[arg]] <- loqed[[arg]]
    }
  }

  tc %<>%
    do(correct.conc(
      .,by=by,nomtimevar=nomtimevar,th=#,
      #tau=tau,tstart=tstart,tend=tend,teval=teval,reg=reg,ss=ss,route=route,method=method
    )) %>%
    ungroup

  cat("\n")

  for(arg in enforce){
    if(arg %in% names(loqed))if(!arg %in% names(tc))tc[[arg]] <- loqed[[arg]]
  }


  # 6. create table with corrections

  cat("Creating correction tables...\n")

  corrtab = tc %>% tab.corr(.,nomtimevar=nomtimevar,by=by)

  # 7. Calculate PK parameters NOT based on lambda_z ON CORRECTED DATA

  cat("Calculating parameters that do not need lambda_z...\n")

  par = tc %>%
    group_by_at(by) %>%
    do(calc.par(
      .#,tau=tau,tstart=tstart,tend=tend,teval=teval,route=route,method=method
    )) %>%
    ungroup


  # 8. Calculate PK parameters that need lambda_z

  cat("Calculating parameters that DO need lambda_z...\n")

  par_all = par %>%
    left_join(
      tc %>% select(!!by, reg, ss, factor, route, loqrule) %>% unique
    ) %>%
    calc.par.th(
    #x=par,
    by=by,th=th,covfile=covfile,dose=dose#,
    #reg=reg,ss=ss,factor=factor,route=route
  )

  cat("Combining all parameters...\n")

  par_all = left_join(ctmax,par_all,by=by)

  # 9. create summary PDFs

  if (is.na(pdfdir)){
    cat("No PDF summaries created\n")
  }else{
    cat(paste("Writing summary PDF documents to folder",pdfdir,"...\n"))
  }

  nca.sum(par_all,corrfile=corrtab,by=by,pdfdir=pdfdir)

  cat("\nWriting results...\n")

  if(!missing(covfile)){
    if(is.character(covfile)){
      if(file.exists(covfile)){
        covfile <- read.csv(covfile)
      }else{
        covfile=get(covfile) # to convert the covfile string to the real data frame
      }
    }
  }

  result=list(
    covariates=covfile,
    half_life=th,
    ct_corr=tc,
    corrections=corrtab,
    pkpar=par_all
  )

  cat("\nDone!\n")

  return(result)

}

check.input <- function(
  x, by=NA, nomtimevar=NA, timevar=NA, depvar=NA,
  bloqvar=NA, loqvar=NA, loqrule=NA,
  includeCmax=NA, exclvar=NA, plotdir=NA, pdfdir=NA, timelab=NA, deplab=NA,
  tau=NA, tstart=NA, tend=NA, teval=NA, covfile=NA, dose=NA, factor=NA, reg=NA, ss=NA,
  route=NA, method=NA
) {

  if('ss' %in% names(x)) ss <- x$ss
  if('route' %in% names(x)) route <- x$route
  if('reg' %in% names(x)) reg <- x$reg
  if('method' %in% names(x)) method <- x$method
  if('loqrule' %in% names(x)) loqrule <- x$loqrule
  if('includeCmax' %in% names(x)) includeCmax <- x$includeCmax
  if('tau' %in% names(x)) tau <- x$tau
  if('tstart' %in% names(x)) tstart <- x$tstart
  if('tend' %in% names(x)) tend <- x$tend
  if('teval' %in% names(x)) teval <- x$teval
  if('factor' %in% names(x)) factor <- x$factor
  # 'ss','route','method','loqrule','includeCmax','tau','tstart','tend','teval','factor'


  chkfile <- data.frame(Errors_Warnings="delete",
                        stringsAsFactors = F
  )

  # 1 by variables

  if ( is.na(by[1]) )
    chkfile=rbind(chkfile,"Error: By variable cannot be empty")
  if ( !is.na(by[1]) & length(setdiff(by,names(x)))>0 )
    chkfile=rbind(chkfile,paste("Error: By variable(s) not in input data:",paste(setdiff(by,names(x)), collapse=' ')))

  # 2 other mandatory input dataset columns (nomtimevar,timevar,depvar,bloqvar,loqvar)

  checknames=c("nomtimevar","timevar","depvar","bloqvar","loqvar")
  checkvalues=c(nomtimevar,timevar,depvar,bloqvar,loqvar)

  misstext=paste(checknames[is.na(checkvalues)],collapse=", ")

  if (length(checknames[is.na(checkvalues)])==1)
    chkfile=rbind(chkfile,paste("Error: Argument",misstext,"is mandatory"))
  if (length(checknames[is.na(checkvalues)])>1)
    chkfile=rbind(chkfile,paste("Error: Arguments",misstext,"are mandatory"))

  if ( !(all(checkvalues[!is.na(checkvalues)]%in%names(x))) )
    chkfile=rbind(chkfile,paste("Error: These column(s) are not in input data:",paste(setdiff(checkvalues[!is.na(checkvalues)],names(x)), collapse=' ')))

  # 3 check optional exclusion variable

  if ( !is.na(exclvar) & !(exclvar%in%names(x)))
    chkfile=rbind(chkfile,paste("Error: Exclusion variable (", exclvar, ") not in input data"))

  # 4 valid LOQ rule number?

  if ( !(all(loqrule%in%c(1,2,3,4))))
    chkfile=rbind(chkfile,"Error: Loqrule argument should be 1, 2, 3 or 4")

  # 5 check optional includeCmax argument

  if ( any(!is.na(includeCmax) & !(includeCmax%in%c("Y","y","N","n"))))
    chkfile=rbind(chkfile,"Error: IncludeCmax argument can only be Y (default) or N")

  # 6 check covariate variable

  if (!(missing(covfile))) {
    if(is.character(covfile)){
      if (!exists(covfile)) {
        chkfile=rbind(chkfile,paste("Error: Covariate dataframe",covfile,"does not exist"))
      } else {
        covfile=get(covfile)  # to convert the cov string to the real data frame
      }
    }
    if (!(dose %in% names(covfile))){
      chkfile=rbind(chkfile,paste("Error: Dose variable (",dose,") not in covariate dataframe"))
    }
  } else {
    chkfile=rbind(chkfile,paste("Error: Argument cov is mandatory"))
  }

  # 7 check factor argument

  if (any(is.na(factor)))
    chkfile=rbind(chkfile,"Warning: Factor argument is empty, assuming a value of 1")

  # 8 check regimen argument

  if ( !(all(reg%in%c("sd","SD","md","MD"))))
    chkfile=rbind(chkfile,"Error: Regimen argument should be SD or MD")

  # 9 check steady state argument

  if ( !(all(ss%in%c("y","Y","n","N"))))
    chkfile=rbind(chkfile,"Error: Steady state argument should be Y or N")

  # 9a check if tau is defined when at steady state (ss="Y"). Warn user for this.

  if ( any(ss %in% c("y","Y") & is.na(tau)))
    chkfile=rbind(chkfile,"Warning: Tau not defined while at steady state, no clearances or volumes will be calculated")

  # 10 check route argument

  if ( !all((route%in%c("ev","EV","ivb","IVB","ivi","IVI"))))
    chkfile=rbind(chkfile,"Error: Route argument should be EV, IVB or IVI")

  # 11 check method argument

  if ( !all((method%in%c(1,2,3))))
    chkfile=rbind(chkfile,"Error: Method argument should be 1, 2 or 3")

  chkfile = chkfile %>% filter(Errors_Warnings!="delete")

  if (nrow(chkfile)>0) {

    print(kable(chkfile))
    cat("\n")
    if (any(grepl("Error",chkfile$Errors_Warnings))) stop("Execution of qPNCA aborted due to errors", call.= F)

  }

  else cat("all OK!\n")

}





