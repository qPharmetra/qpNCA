globalVariables('Errors_Warnings')
#' Perform Non-compartmental Analysis
#'
#' Consecutively executes the following NCA steps:
#' * [correct.loq] impute LOQ values
#' * [est.thalf] calculate lambda_z and half-life
#' * [plot_reg] plot each regression curve
#' * [calc.ctmax] calculate Cmax and Tmax
#' * [correct.time] correct time deviations at critical time points
#' * [correct.conc] impute missing concentrations at critical time points
#' * [tab.corr] tabulate data alterations
#' * [calc.par] calculates parameters not dependent on lambda_z
#' * [calc.par.th] calculates parameters dependent on lambda_z
#'
#' @param x input dataset name
#' @param by character: column names in x indicating grouping variables; default is as.character(dplyr::groups(x))
#' @param nomtimevar variable name containing the nominal sampling time after dose
#' @param timevar variable name containing the actual sampling time after dose
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param loqvar variable name containing the LOQ value
#' @param loqrule rule number to be applied to the LOQ values in the curve; x$loqrule overrides if provided
#' @param includeCmax include Cmax in half-life estimation? (y/n); x$includeCmax overrides if provided
#' @param exclvar variable name indicating points to be excluded in half-life estimation (these should have exclvar = 1)
#' @param plotdir directory where regression plots (.PNG) will be saved; NA gives default location, NULL skips regression plots
#' @param timelab label for time axis in regression plots
#' @param deplab label for dependent variable axis in regression plots
#' @param tau dosing interval (for multiple dosing); NA (default) if single dose; x$tau overrides
#' @param tstart start time of partial AUC (start>0); NA (default) if not requested; x$tstart overrides
#' @param tend end time of partial AUC; NA (default) if not requested; x$tend overrides
#' @param teval user selected AUC interval (starting at t=0); NA (default) if not requested; x$teval overrides
#' @param covariates covariates dataset (containing at least dose for CL calculation); defaults to unique combinations of \code{by} and \code{dose} evaluated on \code{x}; can be character path of csv file or local object
#' @param dose variable containing the dose amount; default 'dose' set to 1 if not in \code{names(x)}
#' @param factor conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000); x$factor overrides if provided
#' @param reg regimen, "SD" or "MD"; x$reg overrides if provided
#' @param ss is steady state reached (y/n); x$ss overrides if provided
#' @param route route of drug administration ("EV","IVB", "IVI"); x$route overrides if provided
#' @param method method for trapezoidal rule; x$method overrides if provided
#' * 1: linear up - linear down
#' * 2: linear up - logarithmic down
#' * 3: linear before first Tmax, logarithmic after first Tmax
#'
#' @return (list)
#' * **covariates** covariates selected with the \code{covariates} argument
#' * **half_life** linear regression parameters
#' * **ct_corr** the time and concentration corrected dataset
#' * **corrections** descriptions of the corrections applied
#' * **pkpar** all estimated PK parameters
#' * **plots** generated plots
#'
#' @export
#' @importFrom utils read.csv
#' @importFrom knitr kable
#' @importFrom dplyr group_by_at ungroup left_join
#' @examples
#' \donttest{
#' library(magrittr)
#' library(dplyr)
#' data(ncx)
#' x <- ncx
#' x %<>% group_by(subject)
#' x %<>% mutate(excl_th=0)
#' covs <- x %>% select(subject, wt, dose) %>% unique
#' z <- qpNCA(x, by = 'subject', covariates = covs, exclvar='excl_th')
#' z %>% names
#' }

qpNCA <- function(
  x,
  by=NULL,
  nomtimevar="ntad",
  timevar="time",
  depvar="dv",
  bloqvar="bloq",
  loqvar="loq",
  loqrule=1,
  includeCmax="Y",
  exclvar=NA,
  plotdir=NA,
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

  message("\n")
  message("Checking function arguments...")

  enforce <- c(
    'ss','route','reg','method','loqrule','includeCmax',
    'tau','tstart','tend','teval','factor'
  )

  enforced <- setdiff(enforce, names(x))

  for(arg in enforce){
    if(!arg %in% names(x)) x[[arg]] <- get(arg)  # overrule column values with attribute
  }

 rm(list = enforce)

 if(is.null(by)) by <- as.character(groups(x))

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

check.input(
    x, by=by, nomtimevar=nomtimevar, timevar=timevar, depvar=depvar,
    bloqvar=bloqvar, loqvar=loqvar,
    #loqrule=loqrule, includeCmax=includeCmax,
    exclvar=exclvar, plotdir=plotdir, timelab=timelab,
    deplab=deplab,
    #tau=tau, tstart=tstart, tend=tend, teval=teval,
    covariates=covariates, dose=dose#,
    #factor=factor, reg=reg, ss=ss,route=route, method=method
  )

  # 1. Apply LOQ rules

  message("Applying LOQ rules...\n")

    loqed <- x %>% correct.loq(
    by = by,
    nomtimevar = nomtimevar,
    timevar = timevar,
    depvar = depvar,
    bloqvar = bloqvar,
    loqvar = loqvar
  )

    # 2. estimate thalf ON UNCORRECTED DATA

  message("Performing Thalf estimation...\n")

  th <- loqed %>% est.thalf(
    by = by,
    timevar = timevar,
    depvar = depvar,
    exclvar = exclvar
  )

  # 2a.

  if(is.null(plotdir)){
    message('Skipping regression plots')
  }else{
    if (is.na(plotdir)){
      message("Creating regression plots in standard output...\n")
    }else{
      message(paste("Writing regression plots to directory",plotdir,"...\n"))
    }
    plotdir <- plot_reg(
      loqed, by = by, th = th, bloqvar = bloqvar, timevar = timevar,
      depvar = depvar, exclvar = exclvar, plotdir = plotdir,
      timelab = timelab, deplab = deplab
    )
  }
  message("\n")

  # 3. find Cmax and tmax ON UNCORRECTED DATA

  message("Calculating Cmax/Tmax...\n")

  ctmax <- loqed %>% calc.ctmax(by = by, timevar = timevar, depvar = depvar)

  # 4. and 5. create dataset with corrected time deviations

  message("Applying time deviation corrections and missing concentration imputations...\n")

  tc <- loqed %>% correct.time(
    by = by,
    nomtimevar = nomtimevar,
    timevar = timevar,
    depvar = depvar,
    th = th
  )

  tc %<>% correct.conc(by = by, nomtimevar = nomtimevar)

  # 6. create table with corrections

  message("Creating correction tables...\n")

  corrtab = tc %>% tab.corr(., by=by, nomtimevar=nomtimevar)

  # 7. Calculate PK parameters NOT based on lambda_z ON CORRECTED DATA

  message("Calculating parameters that do not need lambda_z...\n")

  # par <- tc %>% calc.par(by = by)

  # 8. Calculate PK parameters that need lambda_z
  # @ 1.1.8, tc passed directly to calc.par.th

  message("Calculating parameters that DO need lambda_z...\n")

  # par_all = par %>%
  #   left_join(
  #     tc %>% select(!!by, reg, ss, factor, route, loqrule) %>% unique
  #   ) %>%

  par_all <- tc %>%
    calc.par.th(
    by=by,th=th,covariates=covariates,dose=dose
    #,
    #reg=reg,ss=ss,factor=factor,route=route
  )

  message("Combining all parameters...")

  par_all = left_join(ctmax,par_all,by=by)

   message("\nWriting results...\n")

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
    plots = plotdir
  )

  message("\nDone!\n")

  return(result)

}

#' Check qpNCA function arguments for validity
#'
#' Checks whether all function arguments are valid and entered column names are present in input data \cr
#' See \code{\link{qpNCA}} for description of the arguments
#'
#' @param x data.frame
#' @param by column names in x indicating grouping variables
#' @param nomtimevar variable name containing the nominal sampling time after dose
#' @param timevar variable name containing the actual sampling time after dose
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param loqvar variable name containing the LOQ value
#' @param loqrule rule number to be applied to the LOQ values in the curve; x$loqrule overrides if provided
#' @param includeCmax include Cmax in half-life estimation? (y/n); x$includeCmax overrides if provided
#' @param exclvar variable name indicating points to be excluded in half-life estimation (these should have exclvar = 1)
#' @param plotdir directory where regression plots (.PNG) will be saved; NA gives default location, NULL skips regression plots
#' @param timelab label for time axis in regression plots
#' @param deplab label for dependent variable axis in regression plots
#' @param tau dosing interval (for multiple dosing); NA (default) if single dose; x$tau overrides
#' @param tstart start time of partial AUC (start>0); NA (default) if not requested; x$tstart overrides
#' @param tend end time of partial AUC; NA (default) if not requested; x$tend overrides
#' @param teval user selected AUC interval (starting at t=0); NA (default) if not requested; x$teval overrides
#' @param covariates covariates dataset; Must contain the dose variable
#' @param dose variable containing the dose amount
#' @param factor conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000); x$factor overrides if provided
#' @param reg regimen, "SD" or "MD"; x$reg overrides if provided
#' @param ss is steady state reached (y/n); x$ss overrides if provided
#' @param route route of drug administration ("EV","IVB", "IVI"); x$route overrides if provided
#' @param method method for trapezoidal rule; x$method overrides if provided
#' @return Check results
check.input <- function(
  x, by=NA, nomtimevar=NA, timevar=NA, depvar=NA,
  bloqvar=NA, loqvar=NA, loqrule=NA,
  includeCmax=NA, exclvar=NA, plotdir=NA, timelab=NA, deplab=NA,
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

  # 10 check if tstart, tend and teval differ from 0 (use teval instead of start/tend for interval 0-x)

  if ( any(tstart==0 & !is.na(tstart)))
    chkfile=rbind(chkfile,"Error: Tstart cannot be 0, use Teval to calculate auc(0-Teval)")

  if ( any(tend==0 & !is.na(tend)))
    chkfile=rbind(chkfile,"Error: Tend cannot be 0")

  if ( any(teval==0 & !is.na(teval)))
    chkfile=rbind(chkfile,"Error: Teval cannot be 0")

  # 11 check route argument

  if ( !all((route%in%c("ev","EV","ivb","IVB","ivi","IVI"))))
    chkfile=rbind(chkfile,"Error: Route argument should be EV, IVB or IVI")

  # 12 check method argument

  if ( !all((method%in%c(1,2,3))))
    chkfile=rbind(chkfile,"Error: Method argument should be 1, 2 or 3")

  chkfile = chkfile %>% filter(Errors_Warnings!="delete")

  if (nrow(chkfile)>0) {

    print(kable(chkfile))
    message("\n")
    if (any(grepl("Error",chkfile$Errors_Warnings))) stop("Execution of qPNCA aborted due to errors", call.= F)

  }

  else message("all OK!\n")

}





