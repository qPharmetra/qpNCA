globalVariables('Errors_Warnings')
#' Perform Non-compartmental Analysis
#'
#' Consecutively executes the following NCA steps:
#' * [correct.loq] impute LOQ values
#' * [est.thalf] calculate lambda-z and halflife
#' * [plot_reg] plot each curve
#' * [calc.ctmax] calculate Cmax and Tmax
#' * [correct.time] supply records at critical times with correct concentration
#' * [correct.conc] supply correct concentrations at critical times
#' * [tab.corr] tabulate data alterations
#' * [calc.par] calculates profile-specific summary statistics
#' * [calc.par.th] calculates parameters dependent on lambda-z
#'
#' @param x input dataset name
#' @param by column names in x indicating grouping variables
#' @param nomtimevar variable name containing the nominal sampling time
#' @param timevar variable name containing the sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param loqvar variable name containing the LOQ value
#' @param loqrule rule number to be applied to the LOQ values in the curve; x$loqrule overrides if provided
#' @param includeCmax include results of regression including Cmax in selection? (y/n) x$includeCmax overrides if provided
#' @param exclvar variable name containing information about points to be excluded (these should have exclvar = 1)
#' @param plotdir directory where regression plots (.PNG) will be saved; NA gives default location, NULL skips regression plots
#' @param pdfdir directory where pdf summaries will be saved; NA gives default location, NULL skips summary
#' @param tau dosing interval (for multiple dosing); NA (default) for if single dose; x$tau overrides
#' @param tstart start time of partial AUC (start>0); NA (default) if not requested; x$tstart overrides
#' @param tend end time of partial AUC; NA (default) if not requested; x$tend overrides
#' @param teval user selected AUC interval; NA (default) if not requested; x$teval overrides
#' @param covariates covariates dataset
#' @param dose variable containing the dose amount
#' @param factor conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000); x$factor overrides if provided
#' @param reg regimen, "sd" or "md"; can be character column name in x; x$regimen overrides if provided
#' @param ss is steady state reached (y/n); x$ss overrides if provided
#' @param route route of drug administration ("EV","IVB", "IVI"); x$overrides if provided
#' @param timelab label for time axis
#' @param deplab label for dependent variable axis
#' @param method method for trapezoidal rule;  x$method overrides if provided
#' * 1: linear up - linear down
#' * 2: linear up - logarithmic down
#' * 3: linear before first Tmax, logarithmic after first Tmax
#'
#' @return (list)
#' * **covariates** covariates selected with the \code{cov} argument
#' * **half_life** linear regression parameters
#' * **ctmax** cmax and tmax estimated from uncorrected data
#' * **ct_corr** the time and concentration corrected dataset
#' * **corrections** descriptions of the corrections applied
#' * **pkpar** all estimated PK parameters
#' * **plotdir** directory for plots
#' * **pdfdir** directory for pdf summaries
#'
#' @export
#' @importFrom utils read.csv
#' @importFrom knitr kable
#' @importFrom dplyr group_by_at ungroup left_join
#' @examples
#' library(magrittr)
#' library(dplyr)
#' library(qpNCA)
#' x <- Theoph
#' ntad <- c(0,0.25,0.5,1,2,4,5,7,9,12,24)
#' for(i in 1:nrow(x)){
#'   time  <- x$Time[[i]]
#'   delta <- abs(ntad - time)
#'   best  <- min(delta)
#'   index <- match(best, delta)
#'   nom   <- ntad[[index]]
#'   x$ntad[[i]] <- nom
#' }
#' rm(list = c('time','delta','best','index','nom', 'i','ntad'))
#' x %<>% rename(time = Time, dv = conc, subject = Subject)
#' x %<>% mutate(bloq = 0, loq = 0.01, tad = time)
#' x %<>% filter(dv > 0)
#' covs <- Theoph %>%
#'   select(subject = Subject, Wt, dose = Dose) %>%
#'   unique %>%
#'   mutate(dose = dose * Wt) # see ?Theoph
#' y <- qpNCA(x, by = 'subject', covariates = covs)

qpNCA <- function(
  x,
  by=character(0),
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
  covariates=NA,
  dose='dose',
  factor=1,
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
    covariates=covariates, dose=dose#,
    #factor=factor, reg=reg, ss=ss,route=route, method=method
  )

  # 1. Apply LOQ rules

  cat("Applying LOQ rules...\n")

  # loqed = x %>%
  #   group_by_at(by) %>%
  #   do(
  #     correct.loq(
  #       .,
  #       nomtimevar=nomtimevar,
  #       timevar=timevar,
  #       depvar=depvar,
  #       bloqvar=bloqvar,
  #       loqvar=loqvar#,
  #      # loqrule=loqrule: has been embedded in the dataset at this point
  #     )
  #   ) %>%
  #   ungroup

  loqed <- x %>% correct.loq(
    by = by,
    nomtimevar = nomtimevar,
    timevar = timevar,
    depvar = depvar,
    bloqvar = bloqvar,
    loqvar = loqvar
  )

   # loqed$loqrule <- x$loqrule

    # 2. estimate thalf ON UNCORRECTED DATA

  cat("Performing Thalf estimation...\n")

  th <- loqed %>% est.thalf(
    by = by,
    timevar = timevar,
    depvar = depvar,
    exclvar = exclvar
  )
    # group_by_at(by) %>%
    # do(est.thalf(
    #   .,timevar=timevar,depvar=depvar,
    #   #includeCmax=includeCmax,
    #   exclvar=exclvar
    # )) %>%
    # ungroup

  # 2a.

  if(is.null(plotdir)){
    cat('Skipping regression plots')
  }else{
    if (is.na(plotdir)){
      cat("Creating regression plots in standard output...\n")
    }else{
      cat(paste("Writing regression plots to directory",plotdir,"...\n"))
    }
    plotdir <- plot_reg(
      loqed, by = by, th = th, bloqvar = bloqvar, timevar = timevar,
      depvar = depvar, exclvar = exclvar, plotdir = plotdir,
      timelab = timelab, deplab = deplab
    )
  }
  cat("\n")

  # 3. find Cmax and tmax ON UNCORRECTED DATA

  cat("Calculating Cmax/Tmax...\n")

  # ctmax = loqed %>%
  #   group_by_at(by) %>%
  #   do(calc.ctmax(., timevar=timevar, depvar=depvar)) %>%
  #   ungroup

  ctmax <- loqed %>% calc.ctmax(by = by, timevar = timevar, depvar = depvar)

  # 4. and 5. create dataset with corrected time deviations

  cat("Applying time deviation corrections and missing concentration imputations...\n")

  # tc = loqed %>%
  #   group_by_at(by) %>%
  #   do(correct.time(
  #     .,by=by,nomtimevar=nomtimevar,timevar=timevar,depvar=depvar,th=th#,
  #     #tau=tau,tstart=tstart,tend=tend,teval=teval,reg=reg,method=method
  #   ))
  #
  tc <- loqed %>% correct.time(
    by = by,
    nomtimevar = nomtimevar,
    timevar = timevar,
    depvar = depvar,
    th = th
  )

  # now in correct.time()
  # for(arg in c(
  #   'tau','tstart','tend','teval',
  #   'reg','ss','route','method'
  # )){
  #   if(arg %in% names(loqed))if(!arg %in% names(tc)){
  #     tc[[arg]] <- loqed[[arg]]
  #   }
  # }

  # tc %<>%
  #   do(correct.conc(
  #     .,
  #     by=by,
  #     nomtimevar=nomtimevar
  #     #,th=,
  #     #tau=tau,tstart=tstart,tend=tend,teval=teval,reg=reg,ss=ss,route=route,method=method
  #   )) %>%
  #   ungroup

  tc %<>% correct.conc(by = by, nomtimevar = nomtimevar, th = th)

  cat("\n")

  # for(arg in enforce){
  #   if(arg %in% names(loqed))if(!arg %in% names(tc))tc[[arg]] <- loqed[[arg]]
  # }


  # 6. create table with corrections

  cat("Creating correction tables...\n")

  corrtab = tc %>% tab.corr(., by=by, nomtimevar=nomtimevar)

  # 7. Calculate PK parameters NOT based on lambda_z ON CORRECTED DATA

  cat("Calculating parameters that do not need lambda_z...\n")

  # par = tc %>%
  #   group_by_at(by) %>%
  #   do(calc.par(
  #     .,tau=tau,tstart=tstart,tend=tend,teval=teval,route=route,method=method
  #   )) %>%
  #   ungroup

  par <- tc %>% calc.par(by = by)

  # 8. Calculate PK parameters that need lambda_z

  cat("Calculating parameters that DO need lambda_z...\n")

  par_all = par %>%
    left_join(
      tc %>% select(!!by, reg, ss, factor, route, loqrule) %>% unique
    ) %>%
    calc.par.th(
    #x=par,
    by=by,th=th,covariates=covariates,dose=dose
    #,
    #reg=reg,ss=ss,factor=factor,route=route
  )

  cat("Combining all parameters...\n")

  par_all = left_join(ctmax,par_all,by=by)

  # 9. create summary PDFs

  if (is.null(pdfdir)){
    cat("No PDF summaries created\n")
  }else{
    if(!is.na(pdfdir)){
      cat(paste("Writing summary PDF documents to directory",pdfdir,"...\n"))
    }else{
      cat('Writing summary PDF documents to standard location ...\n')
    }
    pdfdir <- nca.sum(par_all,corrfile=corrtab,by=by,pdfdir=pdfdir)
  }



  cat("\nWriting results...\n")

  if(!missing(covariates)){
    if(is.character(covariates)){
      if(file.exists(covariates)){
        covariates <- read.csv(covariates)
      }else{
        covariates=get(covariates) # to convert the covariates string to the real data frame
      }
    }
  }

  result = list(
    covariates=covariates,
    half_life=th,
    ct_corr=tc,
    corrections=corrtab,
    pkpar=par_all,
    plots = plotdir,
    pdfdir = pdfdir
  )

  cat("\nDone!\n")

  return(result)

}

check.input <- function(
  x, by=NA, nomtimevar=NA, timevar=NA, depvar=NA,
  bloqvar=NA, loqvar=NA, loqrule=NA,
  includeCmax=NA, exclvar=NA, plotdir=NA, pdfdir=NA, timelab=NA, deplab=NA,
  tau=NA, tstart=NA, tend=NA, teval=NA, covariates=NA, dose=NA, factor=NA, reg=NA, ss=NA,
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

  if (!(missing(covariates))) {
    if(is.character(covariates)){
      if (!exists(covariates)) {
        chkfile=rbind(chkfile,paste("Error: Covariate dataframe",covariates,"does not exist"))
      } else {
        covariates=get(covariates)  # to convert the cov string to the real data frame
      }
    }
    if (!(dose %in% names(covariates))){
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





