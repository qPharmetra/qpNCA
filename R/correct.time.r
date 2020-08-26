#' Correct Concentrations for Time Deviations
#'
#' Corrects concentrations at critical, but deviating time points
#' (e.g, predose, TAU, start and end of user selected AUC interval),
#' and adds missing records at these critical time points.
#' * Records with missing NOMINAL time will be removed and this must be corrected before the function is called
#' * If a record at the critical time point is missing and add it and set time to nominal time and set dv conc to NA
#' * Use interpolation if there is a measurable concentration AFTER the nominal time point (i.e. sample is taken too late)
#' * Use extrapilation if there is NO measurable concentration AFTER the nominal time point (i.e. sample is taken too early)
#' * Set deviating time at predose to 0
#' * Original time and conc will be kept in original variables.
#'
#' The following Time Deviation Correction Rules will be applied to critical time points (t=0, tau, tstart, tend, teval), if needed:
#'
#'   **Rule** | **Regimen** | **Description** | **Applied to**
#'   ---- | ------- | ----------- | ----------
#'   SDT-1 | sd | Set actual time to 0 | t=0
#'   SDT-2 | sd | Correct concentration at deviating time by interpolation | t=tau,tstart,tend,teval
#'   SDT-3 | sd | Correct concentration at deviating time by extrapolation | t=tau,tend,teval
#'   MDT-1 | md | If predose sample taken after dosing, set actual time to 0 and conc to NA | t=0
#'   MDT-2 | md | Correct concentration at deviating time by interpolation (too late) |  t=tau,tstart,tend,teval
#'   MDT-3 | md | Correct concentration at deviating time by extrapolation (too early) | t=0,tau,tend,teval
#'   MDT-3a | md | Set actual time to zero if concentration is BLOQ (too early) | t=0
#'
#' @param x input dataset name (contains all data, including LOQ (set conc to zero for these))
#' @param nomtimevar variable name containing the nominal sampling time
#' @param timevar variable name containing the actual sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param tau dosing interval (for multiple dosing), if single dose, leave empty
#' @param teval user selected AUC interval, if not requested, leave empty
#' @param th file name of file with lamdba_z information for each curve (can be derived from est.thalf)
#' @param reg regimen, "sd" or "md"
#' @param by column names in x indicating grouping variables
#' @param tstart column name in x indicating start time
#' @param tend column name in x indicating end time
#' @param method method of interpolation:
#' * 1: linear up - linear down
#' * 2: linear up - logarithmic down
#'
#' @return a dataset with time deviation corrections applied (timevar and depvar adapted). The following variables are added:
#'
#'  **Variable** | **Description**
#'   ---- | -----
#'   create.nr        |         is a missing record created?
#'   create.txt       |         explanation of what is created
#'   trule.nr         |         correction rule number
#'   trule.txt        |         text explaining what was altered
#'   applies.to.time   |        lists all AUCS to which the time deviation rule applies
#'   time.tau, conc.tau  |      time and conc, corrected for AUCtau calculation
#'   time.teval, conc.teval |   time and conc, corrected for AUCteval calculation (AUC0-teval)
#'   time.part, conc.part   |   time and conc, corrected for partial AUC calculation (AUCstart-end, start>0)
#'   time.lastall, conc.lastall | time and conc, corrected for AUClast and AUCall calculation
#'   t0.flag, tau.flag, tstart.flag, tend.flag, teval.flag | flags for what timepoint the correction was needed
#'
#' @export
correct.time <- function(
  x,by="subject",nomtimevar="ntad",timevar="time",
  depvar="dv",tau=NA,tstart=NA,tend=NA,teval=NA,th=NA,reg="SD",method=1
){
  for(arg in c('tau','tstart','tend','teval','reg','method')){
    if(arg %in% names(x)){
      if(!eval(substitute(missing(arg)))){
        warning(arg,' supplied as column overrides like-named argument')
      }
      assign(arg,unique(x[[arg]]))
      x[[arg]] <- NULL
    }
    if(length(get(arg)) > 1) {
      warning(arg, ' has length > 1; only first value will be used')
      assign(arg, get(arg)[[1]])
    }
  }

  data_in=x

  if (!missing(th)) { data_in=left_join(data_in,th%>%select(-no.points,-intercept,-r.squared,-adj.r.squared,-thalf),by=by) }
  data_in$includeCmax.x <- NULL
  data_in$includeCmax.y <- NULL

  data_in = data_in %>%
    mutate(depvar=x[[depvar]],                    # dependent variable                      (internal)
           timevar=x[[timevar]],                  # actual time variable                    (internal)
           ptime=x[[nomtimevar]],                 # nominal time                            (internal)
           create.nr="",                          # is missing record created?
           create.txt="",                         # explanation of what is created
           trule.nr="",                           # correction rule number
           trule.txt="",                          # explanation of time correction
           applies.to.time="",                    # lists all AUCS to which the rule applies
           t0.flag=0,tau.flag=0,tstart.flag=0,tend.flag=0,teval.flag=0, # flags for what timepoint the correction was needed
           missflag=0,                            # flag for missing records
           misstime=NA,                           # time of missing record
           lambda_z=ifelse("lambda_z"%in%names(.),lambda_z,NA)) %>%
    filter(!is.na(x[[nomtimevar]]))               # remove records with no nominal time (must be corrected before)

  data_in=data_in %>% mutate_cond(condition=is.na(timevar),timevar=ptime,
                                  trule.nr="-",
                                  trule.txt=paste("Concentration (",depvar,
                                                  ") at missing actual time at t=",ptime,
                                                  " set to missing",sep=""),
                                  depvar=NA
  )

  # create extra records for missing critical timepoints

  do( for(i in c(0,tau,tstart,tend,teval)) {

    if (!is.na(i)) {
      data_in=data_in %>% mutate(diff=abs(ptime-i)) %>%   # take closest neighbouring records as base for extra record
        mutate_cond(condition=!is.element(i,ptime)&ptime==first(ptime[diff==min(diff)]),
                    misstime=i
        )

      misstime=data_in %>% filter(misstime==i) %>%
        mutate(ptime=i,
               timevar=i,
               depvar=NA,
               missflag=1,
               create.nr="-",
               create.txt=paste("Missing record at t=",misstime," added",sep="")
        )
      data_in=bind_rows(data_in,misstime)
    }
  }

  )

  # Estimate lagging and leading time points and concentrations for each time point

  data_in=lag_lead(data_in,nomtimevar1="ptime",depvar1="depvar",timevar1="timevar",lagc="lagdv",lagt="lagtime",leadc="leaddv",leadt="leadtime")

  # Start time corrections

  data_in=data_in %>%

    # TAU
    mutate(conc.tau=depvar,time.tau=timevar) %>%
    mutate_cond(condition = !is.na(tau)&ptime==tau&timevar<ptime&!is.na(leaddv),

                conc.tau=interpol(c1=depvar, c2=leaddv, t1=tau, t2=timevar, t3=leadtime, method=method), #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.tau=depvar+((leaddv-depvar)*(tau-timevar)/(leadtime-timevar)), #interpolate
                time.tau=tau,
                tau.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tau,
                                " ) corrected to ",round(conc.tau,3),
                                " by interpolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TAU ")
    ) %>%
    mutate_cond(condition = !is.na(tau)&ptime==tau&timevar>ptime&!is.na(lagdv),

                conc.tau=interpol(c1=lagdv, c2=depvar, t1=tau, t2=lagtime, t3=timevar,method=method), #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.tau=lagdv+((depvar-lagdv)*(tau-lagtime)/(timevar-lagtime)),     #interpolate
                time.tau=tau,
                tau.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tau,
                                " ) corrected to ",round(conc.tau,3),
                                " by interpolation (sample taken too late)",sep=""),
                applies.to.time=paste(applies.to.time,"TAU ")
    ) %>%
    mutate_cond(condition = !is.na(tau)&ptime==tau&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
                conc.tau=depvar*exp(-1*lambda_z*(tau-timevar)),                      # extrapolate
                time.tau=tau,
                tau.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tau,
                                " ) corrected to ",round(conc.tau,3),
                                " by extrapolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TAU ")
    ) %>%
    # TEVAL
    mutate(conc.teval=depvar,time.teval=timevar) %>%
    mutate_cond(condition = !is.na(teval)&ptime==teval&timevar<ptime&!is.na(leaddv),

                conc.teval=interpol(c1=depvar, c2=leaddv, t1=teval, t2=timevar, t3=leadtime, method=method), #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.teval=depvar+((leaddv-depvar)*(teval-timevar)/(leadtime-timevar)),  #interpolate
                time.teval=teval,
                teval.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",teval,
                                " ) corrected to ",round(conc.teval,3),
                                " by interpolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TEVAL ")
    ) %>%
    mutate_cond(condition = !is.na(teval)&ptime==teval&timevar>ptime&!is.na(lagdv),

                conc.teval=interpol(c1=lagdv, c2=depvar, t1=teval, t2=lagtime, t3=timevar, method=method), #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.teval=lagdv+((depvar-lagdv)*(teval-lagtime)/(timevar-lagtime)),      #interpolate
                time.teval=teval,
                teval.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",teval,
                                " ) corrected to ",round(conc.teval,3),
                                " by interpolation (sample taken too late)",sep=""),
                applies.to.time=paste(applies.to.time,"TEVAL ")
    ) %>%
    mutate_cond(condition = !is.na(teval)&ptime==teval&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
                conc.teval=depvar*exp(-1*lambda_z*(teval-timevar)),                       # extrapolate
                time.teval=teval,
                teval.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",teval,
                                " ) corrected to ",round(conc.teval,3),
                                " by extrapolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TEVAL ")
    ) %>%
    # PARTIAL
    #   TSTART
    mutate(conc.part=depvar,time.part=timevar) %>%
    mutate_cond(condition = !is.na(tstart)&!is.na(tend)&ptime==tstart&timevar<ptime&!is.na(leaddv),

                conc.part=interpol(c1=depvar, c2=leaddv, t1=tstart, t2=timevar, t3=leadtime, method=method), #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.part=depvar+((leaddv-depvar)*(tstart-timevar)/(leadtime-timevar)),  #interpolate
                time.part=tstart,
                tstart.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tstart,
                                " ) corrected to ",round(conc.part,3),
                                " by interpolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TSTART ")
    ) %>%
    mutate_cond(condition = !is.na(tstart)&!is.na(tend)&ptime==tstart&timevar>ptime&!is.na(lagdv),

                conc.part=interpol(c1=lagdv, c2=depvar, t1=tstart, t2=lagtime, t3=timevar, method=method), #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.part=lagdv+((depvar-lagdv)*(tstart-lagtime)/(timevar-lagtime)),      #interpolate
                time.part=tstart,
                tstart.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tstart,
                                " ) corrected to ",round(conc.part,3),
                                " by interpolation (sample taken too late)",sep=""),
                applies.to.time=paste(applies.to.time,"TSTART ")
    ) %>%
    mutate_cond(condition = !is.na(tstart)&!is.na(tend)&ptime==tstart&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
                conc.part=depvar*exp(-1*lambda_z*(tstart-timevar)),                       # extrapolate
                time.part=tstart,
                tstart.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tstart,
                                " ) corrected to ",round(conc.part,3),
                                " by extrapolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TSTART ")
    ) %>%
    #   TEND
    mutate_cond(condition = !is.na(tstart)&!is.na(tend)&ptime==tend&timevar<ptime&!is.na(leaddv),

                conc.part=interpol(c1=depvar, c2=leaddv, t1=tend, t2=timevar, t3=leadtime, method=method), #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.part=depvar+((leaddv-depvar)*(tend-timevar)/(leadtime-timevar)),    #interpolate
                time.part=tend,
                tend.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tend,
                                " ) corrected to ",round(conc.part,3),
                                " by interpolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TEND ")
    ) %>%
    mutate_cond(condition = !is.na(tstart)&!is.na(tend)&ptime==tend&timevar>ptime&!is.na(lagdv),

                conc.part=interpol(c1=lagdv, c2=depvar, t1=tend, t2=lagtime, t3=timevar, method=method),   #interpolate
                #c1+((c2-c1)*(t1-t2)/(t3-t2))
                #conc.part=lagdv+((depvar-lagdv)*(tend-lagtime)/(timevar-lagtime)),        #interpolate
                time.part=tend,
                tend.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tend,
                                "  ) corrected to ",round(conc.part,3),
                                " by interpolation (sample taken too late)",sep=""),
                applies.to.time=paste(applies.to.time,"TEND ")
    ) %>%
    mutate_cond(condition = !is.na(tstart)&!is.na(tend)&ptime==tend&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
                conc.part=depvar*exp(-1*lambda_z*(tend-timevar)),                         # extrapolate
                time.part=tend,
                tend.flag=1,
                trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
                trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tend,
                                " ) corrected to ",round(conc.part,3),
                                " by extrapolation (sample taken too early)",sep=""),
                applies.to.time=paste(applies.to.time,"TEND ")
    ) %>%
    # T=0            (also needed for AUCpartial since this value may be used for substitution at t=TAU in case TEND=TAU)
    mutate(conc.lastall=depvar,time.lastall=timevar,newdepvar=NA) %>%
    mutate_cond(condition = ptime==0&timevar!=ptime&tolower(reg)=="sd",
                t0.flag=1,
                time.tau=0,time.part=0,time.teval=0,time.lastall=0,                       # set predose time to 0
                conc.tau=depvar,conc.part=depvar,conc.teval=depvar,conc.lastall=depvar,
                trule.nr="SDT-1",
                trule.txt=paste("Deviating time at predose ( t=",timevar," ) set to 0",sep=""),
                applies.to.time="PREDOSE"
    ) %>%


    mutate_cond(condition = ptime==0&timevar<ptime&!is.na(depvar)&depvar>0&tolower(reg)=="md"&!is.na(lambda_z),
                t0.flag=1,
                time.tau=0,time.part=0,time.teval=0,time.lastall=0,
                newdepvar=depvar*exp(-1*lambda_z*(0-timevar)),                            # extrapolate
                conc.tau=newdepvar,conc.part=newdepvar,conc.teval=newdepvar,conc.lastall=newdepvar,
                trule.nr="MDT-3",
                trule.txt=paste("Concentration at predose ( ",depvar," ), taken before dosing ( t=",timevar,
                                " ) corrected to ",round(newdepvar,3),
                                " by extrapolation (sample taken too early)",sep=""),
                applies.to.time="PREDOSE"
    ) %>%

    mutate_cond(condition = ptime==0&timevar<ptime&(is.na(depvar)|depvar==0)&tolower(reg)=="md",
                t0.flag=1,
                time.tau=0,time.part=0,time.teval=0,time.lastall=0,
                newdepvar=depvar,
                conc.tau=newdepvar,conc.part=newdepvar,conc.teval=newdepvar,conc.lastall=newdepvar,
                trule.nr="MDT-3a",
                trule.txt="Time of LOQ or NA value at predose taken before dosing set to 0",
                applies.to.time="PREDOSE"
    ) %>%

    mutate_cond(condition = ptime==0&timevar>ptime&tolower(reg)=="md",
                t0.flag=1,
                time.tau=0,time.part=0,time.teval=0,time.lastall=0,
                conc.tau=NA,conc.part=NA,conc.teval=NA,conc.lastall=NA,                   # set conc to NA
                trule.nr="MDT-1",
                trule.txt=paste("Concentration at predose ( ",depvar," ), taken after dosing ( t=",timevar,
                                " ) set to NA (sample taken too late)",sep=""),
                applies.to.time="PREDOSE"
    ) #%>%
  #       select(-lambda_z,-leaddv,-lagdv,-leadtime,-lagtime, -missflag, -misstime, -diff, -newdepvar)

  result=data_in

  names(result)[names(result)=="ptime"]=nomtimevar
  names(result)[names(result)=="timevar"]=timevar   # to be sure time value of created time points are copied to original time variable
  names(result)[names(result)=="depvar"]=depvar     # to be sure conc value of created time points are copied to original conc variable
  return(result)

}

mutate_cond <- function (.data, condition, ..., envir = parent.frame()){
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
locf <- function(x){
  good <- !is.na(x)
  positions <- seq(length(x))
  good.positions <- good * positions
  last.good.position <- cummax(good.positions)
  last.good.position[last.good.position == 0] <- NA
  x[last.good.position]
}
