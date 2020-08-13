#' Corrects missing concentration at critical time points (e.g, predose, TAU, start and end of user selected AUC interval)
#' @importFrom metrumrg snap locf
#' @description
#' \itemize{Use interpolation if there is a measurable concentration BEFORE and AFTER the missing concentration}\cr
#' \itemize{Use extrapolation if there is NO measurable concentration AFTER the missing concentration}\cr
#' \itemize{Set missing concentration at predose to 0 (SD, non-endogenous) or value at t=TAU (steady state only)}\cr
#' \itemize{Set missing concentration at t=TAU to value at t=0 (steady state only)}\cr
#' \itemize{The following Concentration Deviation Correction Rules will be applied to critical time points (t=0, tau, tstart, tend, teval), if needed:}
#' \tabular{llll}{
#' Rule \tab Regimen \tab Description \tab Applied to \cr
#' SDC-1  \tab   sd   \tab       Set concentration to 0 (only non-endogenous compounds)      \tab      t=0 \cr
#' SDC-2   \tab    sd  \tab      impute missing concentration by interpolation   \tab                  t=tau,tstart,tend,teval\cr
#' SDC-3   \tab    sd    \tab    impute missing concentration by extrapolation   \tab                  t=tau,tend,teval\cr
#' SDC-4   \tab    sd (IV)  \tab impute missing concentration by back-extrapolation   \tab             t=0\cr
#' MDC-1   \tab    md  \tab      impute missing concentration by existing conc at t=0 or t=tau*  \tab  t=0,tau\cr
#' MDC-2   \tab    md    \tab    impute missing concentration by interpolation   \tab                   t=tau,tstart,tend,teval\cr
#' MDC-3   \tab    md    \tab    impute missing concentration by extrapolation   \tab                  t=tau,tend,teval\cr
#' MDC-4   \tab    md (IV)  \tab  impute missing concentration by back-extrapolation  \tab              t=0\cr
#' * only if steady state has been reached\cr
#' }
#' @param x input dataset name input dataset name (contains all data, including LOQ (set conc to zero for these))
#' @param nomtimevar variable name containing the nominal sampling time
#' @param tau dosing interval (for multiple dosing), if single dose, leave empty
#' @param tstart start time of partial AUC (start>0), if not requested, leave empty
#' @param tend end time of partial AUC, if not requested, leave empty
#' @param teval user selected AUC interval, if not requested, leave empty
#' @param th file name of file with lamdba_z information for each curve (can be derived from est.thalf)
#' @param reg regimen, "sd" or "md"
#' @param ss is steady state reached (y/n)
#' @param route route of drug administration ("po","iv")
#' @param method of interpolation: \cr
#'             1: linear up - linear down \cr
#'             2: linear up - logarithmic down \cr
#'
#' @return  a dataset with missing concentrations imputed. The following variables are added:
#' \tabular{cc}{
#'  crule.nr     \tab    correction rule number \cr
#'  crule.txt     \tab   text explaining what was altered \cr
#'  applies.to.conc  \tab  lists all AUCS to which the concentration correction rule applies \cr
#' }
#'
#' @export
#'
correct.conc <- function(
  x,by="subject",nomtimevar="ntad",tau=NA,tstart=NA,
  tend=NA,teval=NA,th=NA,reg="SD",ss="N",route="EV",method=1) {

  #tau=tau,tstart=tstart,tend=tend,teval=teval,reg=reg,ss=ss,route=route,method=method
  for(arg in c('tau','tstart','tend','teval','reg','method','route','ss')){
    if(arg %in% names(x)){
      if(!missing(arg)){
        warning(arg,' supplied as column overrides like-named argument')
      }
      assign(arg,unique(x[[arg]]))
    }
    if(length(get(arg)) > 1) {
      warning(arg, ' has length > 1; only first value will be used')
      assign(arg, get(arg)[[1]])
    }
  }

  data_in=x

  if (!missing(th)) { data_in=left_join(data_in,th%>%select(-no.points,-intercept,-r.squared,-adj.r.squared,-thalf),by=by) }

  data_in=data_in %>%
          mutate(ptime=x[[nomtimevar]],                 # nominal time                            (internal)
                 crule.nr="",                           # correction rule number
                 crule.txt="",                          # explanation of concentration substitution
                 applies.to.conc="",                     # lists all AUCS to which the concentration correction rule applies
                 lambda_z=ifelse("lambda_z"%in%names(.),lambda_z,NA)
                )

# create lead and lag variables for each AUC

   #TAU
if (!is.na(tau)) {
  data_in=lag.lead(data_in,nomtimevar1="ptime",depvar1="conc.tau",timevar1="time.tau",
                   lagc="lag.ctau",lagt="lag.ttau",leadc="lead.ctau",leadt="lead.ttau")
}
   #PARTIAL
if (!is.na(tstart)&!is.na(tend)) {
  data_in=lag.lead(data_in,nomtimevar1="ptime",depvar1="conc.part",timevar1="time.part",
                   lagc="lag.cpart",lagt="lag.tpart",leadc="lead.cpart",leadt="lead.tpart")
   }
   #TEVAL
if (!is.na(teval)) {
  data_in=lag.lead(data_in,nomtimevar1="ptime",depvar1="conc.teval",timevar1="time.teval",
                   lagc="lag.cteval",lagt="lag.tteval",leadc="lead.cteval",leadt="lead.tteval")
   }


  # Start concentration substitutions

  # Create some variables

  data_in=data_in %>% mutate(t0val=ifelse(is.element(0,ptime),conc.tau[ptime==0],NA),
                             tauval=ifelse(is.element(tau,ptime),conc.tau[ptime==tau],NA)
  )

  #TAU
  if (!is.na(tau)) {
    data_in=data_in %>% mutate_cond(condition = ptime==tau&is.na(conc.tau)&!is.na(t0val)
                                    &tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="ivb",
                                    conc.tau=t0val,                                        # take value at t=0
                                    time.tau=tau,
                                    tau.flag=1,
                                    crule.nr="MDC-1",
                                    crule.txt=paste("Missing concentration at t=",tau," corrected to ",round(conc.tau,3),
                                                    " by imputation with existing concentration at t=0",sep=""),
                                    applies.to.conc=paste(applies.to.conc,"TAU ")
    ) %>%
      mutate_cond(condition = ptime==tau&is.na(conc.tau)&!is.na(lag.ctau)&!is.na(lead.ctau),

                  conc.tau=interpol(c1=lag.ctau, c2=lead.ctau, t1=tau, t2=lag.ttau, t3=lead.ttau, method=method), #interpolate
                  #c1+((c2-c1)*(t1-t2)/(t3-t2))
                  #conc.tau=lag.ctau +
                  #         (lead.ctau-lag.ctau)*(tau-lag.ttau)/(lead.ttau-lag.ttau),  # interpolate
                  time.tau=tau,
                  tau.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
                  crule.txt=paste("Missing concentration at t=",tau," corrected to ",round(conc.tau,3),
                                  " by interpolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TAU ")
      ) %>%
      mutate_cond(condition = ptime==tau&is.na(conc.tau)&!is.na(lag.ctau)
                  &is.na(lead.ctau)&!is.na(lambda_z),
                  conc.tau=lag.ctau*exp(-1*lambda_z*(tau-lag.ttau)),  # extrapolate
                  time.tau=tau,
                  tau.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
                  crule.txt=paste("Missing concentration at t=",tau," corrected to ",round(conc.tau,3),
                                  " by extrapolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TAU ")
      ) %>%
      mutate(tauval=conc.tau[time.tau==tau])
  }
  #TSTART and TEND
  #TSTART
  if (!is.na(tstart)&!is.na(tend)) {
    data_in=data_in %>%  mutate_cond(condition = ptime==tstart&!is.na(tau)&tstart==tau&is.na(conc.part)
                                     &!is.na(t0val)&tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="ivb",
                                     conc.part=t0val,                              # take value at t=0 if TSTART=TAU
                                     time.part=tstart,
                                     tstart.flag=1,
                                     crule.nr="MDC-1",
                                     crule.txt=paste("Missing concentration at t=",tstart,
                                                     " corrected to ",round(conc.part,3),
                                                     " by imputation with existing concentration at t=0",sep=""),
                                     applies.to.conc=paste(applies.to.conc,"TSTART ")
    ) %>%
      mutate_cond(condition = ptime==tstart&is.na(conc.part)&!is.na(lag.cpart)&!is.na(lead.cpart),
                  conc.part=interpol(c1=lag.cpart, c2=lead.cpart, t1=tstart, t2=lag.tpart, t3=lead.tpart, method=method), #interpolate
                  #c1+((c2-c1)*(t1-t2)/(t3-t2))
                  #conc.part=lag.cpart +                                         # interpolate
                  #          (lead.cpart-lag.cpart)*(tstart-lag.tpart)/(lead.tpart-lag.tpart),
                  time.part=tstart,
                  tstart.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
                  crule.txt=paste("Missing concentration at t=",tstart,
                                  " corrected to ",round(conc.part,3),
                                  " by interpolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TSTART ")
      ) %>%
      mutate_cond(condition = ptime==tstart&is.na(conc.part)&!is.na(lag.cpart)
                  &is.na(lead.cpart)&!is.na(lambda_z),
                  conc.part=lag.cpart*exp(-1*lambda_z*(tstart-lag.tpart)),       # extrapolate
                  time.part=tstart,
                  tstart.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
                  crule.txt=paste("Missing concentration at t=",tstart,
                                  " corrected to ",round(conc.part,3),
                                  " by extrapolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TSTART ")
      )
  }
  #TEND
  if (!is.na(tstart)&!is.na(tend)) {
    data_in=data_in %>%  mutate_cond(condition = ptime==tend&!is.na(tau)&tend==tau&is.na(conc.part)
                                     &!is.na(t0val)&tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="ivb",
                                     conc.part=t0val,                              # take value at t=0 if TEND=TAU
                                     time.part=tend,
                                     tend.flag=1,
                                     crule.nr="MDC-1",
                                     crule.txt=paste("Missing concentration at t=",tend," corrected to ",
                                                     round(conc.part,3),
                                                     " by imputation with existing concentration at t=0",sep=""),
                                     applies.to.conc=paste(applies.to.conc,"TEND ")
    ) %>%
      mutate_cond(condition = ptime==tend&is.na(conc.part)&!is.na(lag.cpart)&!is.na(lead.cpart),
                  conc.part=interpol(c1=lag.cpart, c2=lead.cpart, t1=tend, t2=lag.tpart, t3=lead.tpart, method=method), #interpolate
                  #c1+((c2-c1)*(t1-t2)/(t3-t2))
                  #conc.part=lag.cpart +                                    # interpolate
                  #          (lead.cpart-lag.cpart)*(tend-lag.tpart)/(lead.tpart-lag.tpart),
                  time.part=tend,
                  tend.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
                  crule.txt=paste("Missing concentration at t=",tend,
                                  " corrected to ",round(conc.part,3),
                                  " by interpolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TEND ")
      ) %>%
      mutate_cond(condition = ptime==tend&is.na(conc.part)&!is.na(lag.cpart)
                  &is.na(lead.cpart)&!is.na(lambda_z),
                  conc.part=lag.cpart*exp(-1*lambda_z*(tend-lag.tpart)),   # extrapolate
                  time.part=tend,
                  tend.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
                  crule.txt=paste("Missing concentration at t=",tend,
                                  " corrected to ",round(conc.part,3),
                                  " by extrapolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TEND ")
      )
  }
  #TEVAL
  if (!is.na(teval)) {
    data_in=data_in %>%  mutate_cond(condition = ptime==teval&!is.na(tau)&teval==tau&is.na(conc.teval)
                                     &!is.na(t0val)&tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="ivb",
                                     conc.teval=t0val,                         # take value at t=0 if TEVAL=TAU
                                     time.teval=teval,
                                     teval.flag=1,
                                     crule.nr="MDC-1",
                                     crule.txt=paste("Missing concentration at t=",teval,
                                                     " corrected to ",round(conc.part,3),
                                                     " by imputation with existing concentration at t=0",sep=""),
                                     applies.to.conc=paste(applies.to.conc,"TEVAL ")
    ) %>%
      mutate_cond(condition = ptime==teval&is.na(conc.teval)&!is.na(lag.cteval)&!is.na(lead.cteval),
                  conc.teval=interpol(c1=lag.cteval, c2=lead.cteval, t1=teval, t2=lag.tteval, t3=lead.tteval, method=method), #interpolate
                  #c1+((c2-c1)*(t1-t2)/(t3-t2))
                  #conc.teval=lag.cteval +                             # interpolate
                  #(lead.cteval-lag.cteval)*(teval-lag.tteval)/(lead.tteval-lag.tteval),
                  time.teval=teval,
                  teval.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
                  crule.txt=paste("Missing concentration at t=",teval,
                                  " corrected to ",round(conc.teval,3),
                                  " by interpolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TEVAL ")
      ) %>%
      mutate_cond(condition = ptime==teval&is.na(conc.teval)&!is.na(lag.cteval)
                  &is.na(lead.cteval)&!is.na(lambda_z),
                  conc.teval=lag.cteval*exp(-1*lambda_z*(teval-lag.tteval)), # extrapolate
                  time.teval=teval,
                  teval.flag=1,
                  crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
                  crule.txt=paste("Missing concentration at t=",teval,
                                  " corrected to ",round(conc.teval,3),
                                  " by extrapolation",sep=""),
                  applies.to.conc=paste(applies.to.conc,"TEVAL ")
      )
  }

  # correct t=0 conc for all aucs where t=0 is needed (for PARTIAL this is NOT needed as it does not start at t=0)

  data_in=data_in %>%

    # back-extrapolate t=0 concentration for all AUCs if route is IVB
    mutate(back_extrap=0,
           lc1=lead(conc.lastall,1),lc2=lead(conc.lastall,2),
           lt1=lead(time.lastall,1),lt2=lead(time.lastall,2),
           firstmeasc=conc.lastall[which(conc.lastall>0)][1],
           firstmeast=ptime[which(conc.lastall>0)][1]) %>%
    # if there are NAs or LOQs between t=0 and first measurable conc, set these equal to first measurable conc
    mutate_cond(condition=tolower(route)=="ivb"&ptime>0&ptime<firstmeast&(is.na(conc.lastall)|conc.lastall==0),
                conc.lastall=firstmeasc,conc.tau=firstmeasc,conc.teval=firstmeasc) %>%
    mutate_cond(condition=tolower(route)=="ivb"&ptime==0,
                back_extrap=ifelse(!is.na(lc1)&!is.na(lc2)&lc1>0&lc2>0&
                                     lc1>lc2,1,0),
                conc.lastall=ifelse(back_extrap==1,exp(log(lc1)+
                                                         (0-lt1)/(lt2-lt1)*
                                                         (log(lc2)-log(lc1))),
                                    firstmeasc),
                t0.flag=1,
                crule.nr=ifelse(tolower(reg)=="sd","SDC-4","MDC-4"),
                crule.txt=paste("Concentration at t=0 back-extrapolated to ",round(conc.lastall,3)," (IV bolus)",sep=""),
                time.tau=0, time.teval=0,
                conc.tau=conc.lastall, conc.teval=conc.lastall
    ) %>%

    # ALL and LAST
    mutate_cond(condition = ptime==0&(is.na(conc.lastall)|conc.lastall>0)&tolower(reg)=="sd"&tolower(route)!="ivb",
                conc.lastall=0,                                   # set conc to 0
                time.lastall=0,
                t0.flag=1,
                crule.nr="SDC-1",
                crule.txt=paste("Missing or measurable concentration at (SD) PREDOSE set to 0",sep="")
    ) %>%
    mutate_cond(condition = ptime==0&is.na(conc.lastall)&!is.na(tauval)
                &tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="ivb",
                conc.lastall=tauval,                              # take value at t=tau
                time.lastall=0,
                t0.flag=1,
                crule.nr="MDC-1",
                crule.txt=paste("Missing concentration at PREDOSE corrected to ",round(conc.tau,3),
                                " by imputation with existing concentration at t=TAU",sep=""),
                applies.to.conc=paste("PREDOSE")
    ) %>%
    #TAU
    mutate_cond(condition = ptime==0&(is.na(conc.tau)|conc.tau>0)&tolower(reg)=="sd"&tolower(route)!="ivb",
                conc.tau=0,                                       # set conc to 0
                time.tau=0,
                t0.flag=1,
                crule.nr="SDC-1",
                crule.txt=paste("Missing or measurable concentration at (SD) PREDOSE set to 0",sep=""),
                applies.to.conc=paste("PREDOSE")
    ) %>%
    mutate_cond(condition = ptime==0&is.na(conc.tau)&!is.na(tauval)&tolower(reg)=="md"
                &tolower(ss)=="y"&tolower(route)!="ivb",
                conc.tau=tauval,                                  # take value at t=tau
                time.tau=0,
                t0.flag=1,
                crule.nr="MDC-1",
                crule.txt=paste("Missing concentration at PREDOSE corrected to ",round(conc.tau,3),
                                " by imputation with existing concentration at t=TAU",sep=""),
                applies.to.conc=paste("PREDOSE")
    ) %>%
    #TEVAL
    mutate_cond(condition = ptime==0&(is.na(conc.teval)|conc.teval>0)&tolower(reg)=="sd"&tolower(route)!="ivb",
                conc.teval=0,                                      # set conc to 0
                time.teval=0,
                t0.flag=1,
                crule.nr="SDC-1",
                crule.txt=paste("Missing or measurable concentration at (SD) PREDOSE set to 0",sep=""),
                applies.to.conc=paste("PREDOSE")
    ) %>%
    mutate_cond(condition = ptime==0&is.na(conc.teval)&!is.na(tauval)&tolower(reg)=="md"
                &tolower(ss)=="y"&tolower(route)!="ivb",
                conc.teval=tauval,                                 # take value at t=tau
                time.teval=0,
                t0.flag=1,
                crule.nr="MDC-1",
                crule.txt=paste("Missing concentration at PREDOSE corrected to ",round(conc.teval,3),
                                " by imputation with existing concentration at t=TAU",sep=""),
                applies.to.conc=paste("PREDOSE")
    )
  #                                 %>%
  #                       select(-lag.ctau,-lag.ttau,-lead.ctau,-lead.ttau,
  #                              -lag.cteval,-lag.tteval,-lead.cteval,-lead.tteval,
  #                              -lag.cpart,-lag.tpart,-lead.cpart,-lead.tpart,
  #                              -lambda_z,-ptime, -t0val, -tauval)

  result=data_in
  return(result)
}
