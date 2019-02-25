#Original script for this version of qPNCA ported from svn. Commenting out as separated into one script for each function.

# #######################################################################################
# #
# # VERSION 1.0.11
# #
# # Changed in this version:
# #
# # Added route option (PO/IV)
# #
# #######################################################################################
#
# #### Functions ########################################################################
# #
# # CALC.CTMAX: calculate Cmax, Tmax from raw data for each PK curve (define using group_by)
# #
# # Input dataset:
# # - contains all uncorrected data, including LOQ
# #
# #
# # USAGE:
# #
# # calc.ctmax(x, timevar="time", depvar="dv")
# #
# # ARGUMENTS:
# #
# # x:       input dataset name (if called within dplyr: .)
# # timevar: variable name containing the sampling time
# # depvar:  variable name containing the dependent variable (e.g., concentration)
# #
# # METHOD:
# #
# # Estimate first occurence of maximum concentration for each PK curve
# # If all concentrations are NA, set Cmax and Tmax also to NA
# #
# # OUTPUT:
# #
# # dataset with estimates for the following parameters, one observation per subject:
# #
# #   - cmax     : maximum concentration
# #   - tmax     : time of first occurence of cmax
#
# calc.ctmax <- function(x,timevar="time",depvar="dv") {
#   result=x %>% mutate(depvar=x[[depvar]],        # calculated dependent variable           (internal)
#                       timevar=x[[timevar]]       # calculated time variable                (internal)
#   ) %>%
#     filter(!is.na(depvar))            # na.rm=T in the max function does not work?
#
#   if (dim(result)[1]==0) {                          # if all concentrations are NA, set Cmax and Tmax to NA
#     result=result%>%summarise(cmax=NA,
#                               tmax=NA)
#   }
#   else {
#     result=result%>% summarise(cmax=max(depvar),
#                                tmax=first(timevar[depvar==cmax]) # there might be more than 1 timepoint with Cmax
#     )
#   }
#
#   return(result)
# }
#
# # CALC.PAR: calculate AUCs, tlast, clast.obs for each PK curve (define using group_by)
# #
# # Input dataset:
# # - contains all data after time/concentration deviation corrections
# #
# #
# # USAGE:
# #
# # calc.par(x, timevar="time", depvar="dv", tau=, tend=, tend=, teval=, method=)
# #
# # ARGUMENTS:
# #
# # x:       input dataset name (if called within dplyr: .)
# # tau:     dosing interval (for multiple dosing), if single dose, leave empty
# # tstart:  starting time of user defined interval, if not requested, leave empty
# # tend:    end time of user defined interval, if not requested, leave empty
# # teval:   user selected AUC interval, if not requested, leave empty
# # route:   route of drug administration ("po","iv")
# # method:  method for trapezoidal rule:
# #          1: linear up - linear down
# #          2: linear up - logarithmic down
# #          3: linear before first Tmax, logarithmic after first Tmax
# #
# # OUTPUT:
# #
# # dataset with estimates for the following parameters, one observation per subject:
# #
# #   - t0.ok     : flags if t=0 concentration could be corrected/imputes. If not, no AUCs starting at t=0 are calculated
# #   - tlast.ok  : flags if there is at least one measurable concentration. If not, no AUClast can be calculated
# #   - tlast     : time of last sample with measurable concentration
# #   - clast.obs : observed concentration at tlast
# #   - aucall    : auc calculated over all observations, including values below LOQ (which are set to 0)
# #   - auclast   : auc calculated using all observations up to and including the last measurable concentration (clast.obs at tlast)
# #   - aumcall   : aumc calculated over all observations, including values below LOQ (which are set to 0)
# #   - aumclast  : aumc calculated using all observations up to and including the last measurable concentration (clast.obs at tlast)
# #   - tau       : the dosing interval (if specified)
# #   - calc.tau  : flags if AUCtau could be calculated
# #   - auctau    : auc calculated over the dosing interval, only calculated if tau is specified
# #   - aumctau   : aumc calculated over the dosing interval, only calculated if tau is specified
# #   - teval     : user selected AUC interval starting at t=0 (if specified)
# #   - calc.teval: flags if AUCteval could be calculated
# #   - aucxx     : auc calculated from t=0 up to/including teval, only calculated if teval is specified (xx is substituted by teval)
# #   - calc.part : flags if AUCpart could be calculated
# #   - tstart    : start time of partial AUC (if specified)
# #   - tend      : end time of partial AUC (if specified)
# #   - aucx_y    : partial auc from time=x up to/including time=y, where x>0, only calculated if tstart and tend are specified
# #   - c0        : back-extrapolated concentration at t=0 for IV bolus administration
#
# calc.par <- function(x,tau=NA,tstart=NA,tend=NA,teval=NA,route="po",method=1){
#
#   # Calculate parameters
#
#   # check for each AUC if start and end time/conc are available
#
#   tlast.ok=0         # internal variable
#   t0.ok=0            # internal variable
#   tau.ok=0           # internal variable
#   teval.ok=0         # internal variable
#   part.ok=0          # internal variable
#
#   par = x %>% filter(!is.na(time.lastall))     # if critical time point, it's not NA anymore, so all time=NA are deleted
#
#   if (max(par$conc.lastall,na.rm=T)>0) { tlast.ok=1 }
#
#   if (is.element(0,par$ptime)) { t0.ok=ifelse(par$time.lastall[par$ptime==0]==0
#                                               &!is.na(par$conc.lastall[par$ptime==0]),1,t0.ok) }
#
#   if (!is.na(tau)&t0.ok==1&is.element(tau,par$time.tau))       { tau.ok= ifelse(is.na(par$conc.tau[par$time.tau==0]) |
#                                                                                   is.na(par$conc.tau[par$time.tau==tau]),0,1)  }
#   if (!is.na(teval)&t0.ok==1&is.element(teval,par$time.teval)) { teval.ok=
#     ifelse(is.na(par$conc.teval[par$time.teval==0]) |
#              is.na(par$conc.teval[par$time.teval==teval]),0,1)}
#   if (!is.na(tstart)&!is.na(tend)&is.element(tstart,par$time.part)&is.element(tend,par$time.part))
#   {  part.ok=  ifelse(is.na(par$conc.part[par$time.part==tstart]) |
#                         is.na(par$conc.part[par$time.part==tend]),0,1) }
#
#   par=par %>%    summarise( route=route,
#                             method=method,
#                             tlast=ifelse(tlast.ok==1,last(time.lastall[!is.na(conc.lastall)&conc.lastall>0]),NA),
#                             clast.obs=ifelse(tlast.ok==1,conc.lastall[time.lastall==tlast],NA),
#                             tlast.ok=tlast.ok,
#                             t0.ok=t0.ok,
#
#                             aucall=ifelse(t0.ok==1,
#                                           trap(x=time.lastall[!is.na(conc.lastall)],
#                                                y=conc.lastall[!is.na(conc.lastall)], method=method),NA),
#
#                             auclast=ifelse(t0.ok==1&tlast.ok==1,
#                                            trap(x=time.lastall[!is.na(conc.lastall)&time.lastall<=tlast],
#                                                 y=conc.lastall[!is.na(conc.lastall)&time.lastall<=tlast], method=method),NA),
#
#                             aumcall=ifelse(t0.ok==1&tlast.ok==1,
#                                            trap(x=time.lastall[!is.na(conc.lastall)],
#                                                 y=time.lastall[!is.na(conc.lastall)] *
#                                                   conc.lastall[!is.na(conc.lastall)], method=method),NA),
#
#                             aumclast=ifelse(t0.ok==1&tlast.ok==1,
#                                             trap(x=time.lastall[!is.na(conc.lastall)&time.lastall<=tlast],
#                                                  y=conc.lastall[!is.na(conc.lastall)&time.lastall<=tlast]*
#                                                    time.lastall[!is.na(conc.lastall)&time.lastall<=tlast], method=method),NA),
#
#                             calc.tau=tau.ok,
#                             auctau=ifelse(t0.ok==1&tau.ok==1,
#                                           trap(x=time.tau[!is.na(conc.tau)&time.tau<=tau],
#                                                y=conc.tau[!is.na(conc.tau)&time.tau<=tau], method=method),NA),
#                             aumctau=ifelse(t0.ok==1&tau.ok==1,
#                                            trap(x=time.tau[!is.na(conc.tau)&time.tau<=tau],
#                                                 y=conc.tau[!is.na(conc.tau)&time.tau<=tau]*
#                                                   time.tau[!is.na(conc.tau)&time.tau<=tau], method=method),NA),
#                             tau=ifelse(!is.na(tau),tau,NA),
#
#                             calc.teval=teval.ok,
#                             aucteval=ifelse(t0.ok==1&teval.ok==1,
#                                             trap(x=time.teval[!is.na(conc.teval)&time.teval<=teval],
#                                                  y=conc.teval[!is.na(conc.teval)&time.teval<=teval], method=method),NA),
#                             teval=ifelse(!is.na(teval),teval,NA),
#
#                             calc.part=part.ok,
#                             aucpart=ifelse(part.ok==1,
#                                            trap(x=time.part[!is.na(conc.part)&time.part>=tstart&time.part<=tend],
#                                                 y=conc.part[!is.na(conc.part)&time.part>=tstart&time.part<=tend], method=method),NA),
#                             tstart=ifelse(!is.na(tstart),tstart,NA),
#                             tend=ifelse(!is.na(tend),tend,NA),
#                             c0=ifelse(tolower(route)=="iv",conc.lastall[ptime==0],NA),
#
#
#                             area.back.extr=ifelse(tolower(route)=="iv",
#                                                   trap(x=time.lastall[time.lastall<=firstmeast&!is.na(conc.lastall)],
#                                                        y=conc.lastall[time.lastall<=firstmeast&!is.na(conc.lastall)],
#                                                        method=method),NA)
#   )
#
#   result=par
#
#   names(result)[names(result)=="aucteval"]=ifelse(!is.na(teval),paste("auc",teval,sep=""),"aucteval") # rename aucteval to aucXX
#   names(result)[names(result)=="aucpart"]=ifelse(!is.na(tend),paste("auc",tstart,"_",tend,sep=""),"aucpart") # rename aucpart to aucXX_XX
#   if(is.na(tau))   result = subset(result, select = -c(tau,calc.tau,auctau,aumctau) )    # drop au(m)ctau and tau if not requested
#   if(is.na(teval)) result = subset(result, select = -c(teval,calc.teval,aucteval) )      # drop aucteval and teval if not requested
#   if(is.na(tend))  result = subset(result, select = -c(tstart,tend,calc.part,aucpart) )  # drop aucpart, tstart and tend if not requested
#
#   return(result)
# }
#
# # CALC.PAR.TH: calculate PK parameters that need lambda_z
# #
# # Input dataset:
# # - result parameter dataset from calc.par
# # - result dataset from est.thalf
# # - co-variates dataset (containing at least dose for CL calculation)
# #
# # USAGE:
# #
# # calc.par.th(x=par,th=th,cov=cov,dose=dose,factor=1,reg="sd")
# #
# # ARGUMENTS:
# #
# # x     : result parameter dataset from calc.par (if called within dplyr: .)
# # th    : result dataset from est.thalf
# # cov   : co-variates dataset
# # dose  : variable containing the dose amount
# # factor: conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000)
# # reg   : regimen, "sd" or "md"
# # ss    : is steady state reached (y/n)
# #
# # METHOD:
# #
# # -
# #
# # OUTPUT:
# #
# # dataset with estimates for the following parameters, one observation per subject:
# #
# #   - all parameters calculated in th
# #   - all parameters calculated in par
# #   - clast.pred: predicted concentration at tlast
# #   - aucinf.obs: aucinf based on observed concentration at tlast
# #   - aucinf.pred: aucinf based on predicted concentration at tlast
# #   - aumcinf.obs: area under the first moment curve extrapolated to infinity, based on observed concentration at tlast
# #   - aumcinf.pred: area under the first moment curve extrapolated to infinity, based on predicted concentration at tlast
# #   - cl.f.obs: clearance based on aucinf.obs, at steady state based on auctau
# #   - cl.f.pred: clearance based on aucinf.pred
# #   - mrt.obs: Mean residence time based on aumcinf.obs and aucinf.obs
# #   - mrt.pred: Mean residence time based on aumcinf.pred and aucinf.pred
# #   - vz.f.obs: distribution volume based on cl.f.obs, at steady state based on auctau
# #   - vz.f.pred: distribution based on cl.f.pred
# #   - vss.obs: Steady-state volume based on cl.obs and mrt.obs
# #   - vss.pred: Steady-state volume based on cl.pred and mrt.pred
# #   - regimen (reg)
# #   - steady state reached Y/N? (ss)
# #
# # WARNING: ctmax must be merged separately as these were calculated from uncorrected data
# #
# calc.par.th <- function(x=par,th=th,cov=cov,dosevar="dose",factor=1, reg="sd", ss="n") {
#
#   result=left_join(x,th) %>%
#     left_join(cov) %>%
#     mutate(dosevar=.[[dosevar]],
#            reg=tolower(reg),
#            ss=tolower(ss),
#            clast.pred=exp(intercept-lambda_z*tlast),
#            aucinf.obs=auclast+clast.obs/lambda_z,
#            aucinf.pred=auclast+clast.pred/lambda_z,
#            aumcinf.obs=aumclast+tlast*clast.obs/lambda_z + clast.obs/lambda_z^2,
#            aumcinf.pred=aumclast+tlast*clast.pred/lambda_z + clast.pred/lambda_z^2,
#            cl.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/aucinf.obs, dosevar*factor/auctau),
#            cl.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/aucinf.pred, dosevar*factor/auctau),
#            mrt.obs= ifelse(tolower(ss)=="n", aumcinf.obs/aucinf.obs,
#                            (aumctau + tau*(aucinf.obs-auctau))/auctau),
#            mrt.pred= ifelse(tolower(ss)=="n", aumcinf.pred/aucinf.pred,
#                             (aumctau + tau*(aucinf.pred-auctau))/auctau),
#            vz.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.obs),
#                              dosevar*factor/(lambda_z*auctau)),
#            vz.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.pred),NA),
#            vss.obs= mrt.obs*cl.f.obs,
#            vss.pred= mrt.pred*cl.f.pred,
#            pctextr.obs=(clast.obs/lambda_z)/aucinf.obs*100,
#            pctextr.pred=(clast.pred/lambda_z)/aucinf.pred*100,
#            pctback.obs=area.back.extr/aucinf.obs*100,
#            pctback.pred=area.back.extr/aucinf.pred*100
#     ) %>%
#     select(-dosevar)
#
#   return(result)
#
# }
#
# # EST.THALF: calculate lambda_z and thalf for each PK curve (define using group_by)
# #
# # Input dataset:
# # - contains all uncorrected data, including LOQ
# #
# #
# # USAGE:
# #
# # est.thalf(x, timevar="time", depvar="dv", includeCmax="Y")
# #
# # ARGUMENTS:
# #
# # x:       input dataset name (if called within dplyr: .)
# # timevar: variable name containing the sampling time
# # depvar:  variable name containing the dependent variable (e.g., concentration)
# # includeCmax: include results of regression including Cmax in selection? (y/n)
# #
# # METHOD:
# #
# # The function starts with the last three sample points and performs log-linear regression on it. It then adds one sampling point
# # at a time (including and ending at tmax) and performs the regression again.
# # Internal variable EST checks whether there are at least 3 timepoints for estimation after
# # removal of LOQs and NAs
# #
# # OUTPUT:
# #
# # dataset with estimates for each regression analysis in one observation. The following parameters are available:
# #
# #   - no.points:     number of data points used in the regression analysis
# #   - intercept:     estimated intercept
# #   - lambda_z :     -1*estimated slope
# #   - r.squared:     square of the correlation coefficient
# #   - adj.r.squared: adjusted square of the correlation coefficient
# #   - thalf:         elimination half-life
# #   - start_th:      time of first sample included in the thalf estimation
# #   - end_th:        time of last sample included in the thalf estimation
# #   - includeCmax:   include results of regression including Cmax in selection? (y/n)
#
# est.thalf <- function(x,timevar="time",depvar="dv",includeCmax="Y"){
#   data_in = x %>% mutate(timevar=x[[timevar]],
#                          depvar=x[[depvar]]) %>%
#     filter(!is.na(depvar)&depvar>0) %>%
#     filter(timevar>=first(timevar[depvar==max(depvar,na.rm=T)]))      # include Cmax
#
#   if (tolower(includeCmax)=="n") {
#     data_in=data_in %>% filter(timevar>first(timevar[depvar==max(depvar,na.rm=T)]))   # exclude Cmax
#   }
#
#   est=0
#   i=length(data_in$timevar)-2
#   result = data.frame(matrix(ncol=7,nrow = 1))
#   if (i>=1) {
#     est=1
#     result = data.frame(matrix(ncol=7,nrow = i))
#   }
#
#   while (i>=1) {
#     ## make subset of data frame for each number of data points
#     xx = data_in[i:nrow(data_in), ]
#     ## execute loglinear model fit
#     lmcall=lm(log(xx$depvar)~xx$timevar)
#     lmcall.estimates = as.numeric(summary(lmcall)$coef[,"Estimate"])
#     ## save results of loglin lm fit
#     result[i,1]=length(data_in$timevar)-i+1
#     result[i,2]=lmcall.estimates[1]
#     result[i,3]=lmcall.estimates[2]*-1
#     result[i,4]=summary(lmcall)$r.squared
#     result[i,5]=summary(lmcall)$adj.r.squared
#     result[i,6]=data_in$timevar[i]
#     result[i,7]=last(data_in$timevar)
#     i=i-1
#   }
#
#   names(result) = Cs(no.points,intercept,lambda_z,r.squared,adj.r.squared,start_th,end_th)
#   if (est==1) {
#     result=result %>% mutate(sel=no.points[adj.r.squared==max(adj.r.squared)],
#                              thalf=log(2)/lambda_z,
#                              includeCmax=includeCmax
#     ) %>%
#       filter(sel==no.points) %>%
#       select(-sel)
#   }
#   else {
#     result=result %>% mutate(no.points=NA,intercept=NA,lambda_z=NA,r.squared=NA,adj.r.squared=NA,thalf=NA,
#                              start_th=NA,end_th=NA)
#   }
#   return(result)
# }
#
# # CORRECT.TIME: correct concentration at critical, but deviating, time points (e.g, predose, TAU, start and end of user selected AUC
# #               interval), AND add missing records at these critical time points
# #
# # Input dataset:
# #
# # - contains all data, including LOQ (set conc to zero for these)
# # - dataset containing lambda.z, created using the est.thalf function
# #
# # USAGE:
# #
# # correct.time(x,nomtimevar="ntad",timevar="time",depvar="dv",tau=NA,tstart=NA,tend=NA,teval=NA,th=th,reg="sd",method=1)
# #
# # ARGUMENTS:
# #
# # x         : input dataset name (if called within dplyr: .)
# # nomtimevar: variable name containing the nominal sampling time
# # timevar   : variable name containing the actual sampling time
# # depvar    : variable name containing the dependent variable (e.g., concentration)
# # tau       : dosing interval (for multiple dosing), if single dose, leave empty
# # teval     : user selected AUC interval, if not requested, leave empty
# # th        : file name of file with lamdba_z information for each curve
# # reg       : regimen, "sd" or "md"
# # method    : method of interpolation:
# #             1: linear up - linear down
# #             2: linear up - logarithmic down
# #
# # METHOD:
# #
# # Records with missing NOMINAL time will be removed, this must be corrected before the function is called
# # If record at the critical time point is missing, add it and set time to nominal time and dv conc to NA
# # If there is a measurable concentration AFTER the nominal time point (i.e. sample is taken too late), use interpolation
# # If there is NO measurable concentration AFTER the nominal time point (i.e. sample is taken too early), use extrapolation
# # Set deviating time at predose to 0
# # Original time and conc will be kept in original variables
# #
# # OUTPUT:
# #
# # dataset with time deviation corrections applied (timevar and depvar adapted).
# # The following variables are added:
# #
# #   - create.nr                 : is a missing record created?
# #   - create.txt                : explanation of what is created
# #   - trule.nr                  : correction rule number
# #   - trule.txt                 : text explaining what was altered
# #   - applies.to.time           : lists all AUCS to which the time deviation rule applies
# #   - time.tau, conc.tau        : time and conc, corrected for AUCtau calculation
# #   - time.teval, conc.teval    : time and conc, corrected for AUCteval calculation (AUC0-teval)
# #   - time.part, conc.part      : time and conc, corrected for partial AUC calculation (AUCstart-end, start>0)
# #   - time.lastall, conc.lastall: time and conc, corrected for AUClast and AUCall calculation
# #
# #   - t0.flag, tau.flag, tstart.flag, tend.flag, teval.flag: flags for what timepoint the correction was needed
# #
# # The following Time Deviation Correction Rules will be applied to critical time points (t=0, tau, tstart, tend, teval), if needed:
# #
# # Rule nr.   Regimen    Description                                                                 Applied to
# #
# # SDT-1      sd         Set actual time to 0                                                        t=0
# # SDT-2      sd         Correct concentration at deviating time by interpolation                    t=tau,tstart,tend,teval
# # SDT-3      sd         Correct concentration at deviating time by extrapolation                    t=tau,tend,teval
# #
# # MDT-1      md         if predose sample taken after dosing, set actual time to 0 and conc to NA   t=0
# # MDT-2      md         Correct concentration at deviating time by interpolation (too late)         t=tau,tstart,tend,teval
# # MDT-3      md         Correct concentration at deviating time by extrapolation (too early)        t=0,tau,tend,teval
# # MDT-3a     md         Set actual time to zero if concentration is BLOQ (too early)                t=0
#
# correct.time <- function(x,nomtimevar="ntad",timevar="time",depvar="dv",tau=NA,tstart=NA,tend=NA,teval=NA,th=NA,reg="sd",method=1) {
#
#   data_in=x
#
#   if (!missing(th)) { data_in=left_join(data_in,th%>%select(-no.points,-intercept,-r.squared,-adj.r.squared,-thalf)) }
#
#   data_in = data_in %>%
#     mutate(depvar=x[[depvar]],                    # dependent variable                      (internal)
#            timevar=x[[timevar]],                  # actual time variable                    (internal)
#            ptime=x[[nomtimevar]],                 # nominal time                            (internal)
#            create.nr="",                          # is missing record created?
#            create.txt="",                         # explanation of what is created
#            trule.nr="",                           # correction rule number
#            trule.txt="",                          # explanation of time correction
#            applies.to.time="",                    # lists all AUCS to which the rule applies
#            t0.flag=0,tau.flag=0,tstart.flag=0,tend.flag=0,teval.flag=0, # flags for what timepoint the correction was needed
#            missflag=0,                            # flag for missing records
#            misstime=NA,                           # time of missing record
#            lambda_z=ifelse("lambda_z"%in%names(.),lambda_z,NA)) %>%
#     filter(!is.na(x[[nomtimevar]]))               # remove records with no nominal time (must be corrected before)
#
#   data_in=data_in %>% mutate_cond(condition=is.na(timevar),timevar=ptime,
#                                   trule.nr="-",
#                                   trule.txt=paste("Concentration (",depvar,
#                                                   ") at missing actual time at t=",ptime,
#                                                   " set to missing",sep=""),
#                                   depvar=NA
#   )
#
#   # create extra records for missing critical timepoints
#
#   do( for(i in c(0,tau,tstart,tend,teval)) {
#
#     if (!is.na(i)) {
#       data_in=data_in %>% mutate(diff=abs(ptime-i)) %>%   # take closest neighbouring records as base for extra record
#         mutate_cond(condition=!is.element(i,ptime)&ptime==first(ptime[diff==min(diff)]),
#                     misstime=i
#         )
#
#       misstime=data_in %>% filter(misstime==i) %>%
#         mutate(ptime=i,
#                timevar=i,
#                depvar=NA,
#                missflag=1,
#                create.nr="-",
#                create.txt=paste("Missing record at t=",misstime," added",sep="")
#         )
#       data_in=bind_rows(data_in,misstime)
#     }
#   }
#
#   )
#
#   # Estimate lagging and leading time points and concentrations for each time point
#
#   data_in=lag.lead(data_in,nomtimevar1="ptime",depvar1="depvar",timevar1="timevar",lagc="lagdv",lagt="lagtime",leadc="leaddv",leadt="leadtime")
#
#   # Start time corrections
#
#   data_in=data_in %>%
#
#     # TAU
#     mutate(conc.tau=depvar,time.tau=timevar) %>%
#     mutate_cond(condition = !is.na(tau)&ptime==tau&timevar<ptime&!is.na(leaddv),
#
#                 conc.tau=interpol(c1=depvar, c2=leaddv, t1=tau, t2=timevar, t3=leadtime, method=method), #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.tau=depvar+((leaddv-depvar)*(tau-timevar)/(leadtime-timevar)), #interpolate
#                 time.tau=tau,
#                 tau.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tau,
#                                 " ) corrected to ",round(conc.tau,3),
#                                 " by interpolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TAU ")
#     ) %>%
#     mutate_cond(condition = !is.na(tau)&ptime==tau&timevar>ptime&!is.na(lagdv),
#
#                 conc.tau=interpol(c1=lagdv, c2=depvar, t1=tau, t2=lagtime, t3=timevar,method=method), #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.tau=lagdv+((depvar-lagdv)*(tau-lagtime)/(timevar-lagtime)),     #interpolate
#                 time.tau=tau,
#                 tau.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tau,
#                                 " ) corrected to ",round(conc.tau,3),
#                                 " by interpolation (sample taken too late)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TAU ")
#     ) %>%
#     mutate_cond(condition = !is.na(tau)&ptime==tau&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
#                 conc.tau=depvar*exp(-1*lambda_z*(tau-timevar)),                      # extrapolate
#                 time.tau=tau,
#                 tau.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tau,
#                                 " ) corrected to ",round(conc.tau,3),
#                                 " by extrapolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TAU ")
#     ) %>%
#     # TEVAL
#     mutate(conc.teval=depvar,time.teval=timevar) %>%
#     mutate_cond(condition = !is.na(teval)&ptime==teval&timevar<ptime&!is.na(leaddv),
#
#                 conc.teval=interpol(c1=depvar, c2=leaddv, t1=teval, t2=timevar, t3=leadtime, method=method), #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.teval=depvar+((leaddv-depvar)*(teval-timevar)/(leadtime-timevar)),  #interpolate
#                 time.teval=teval,
#                 teval.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",teval,
#                                 " ) corrected to ",round(conc.teval,3),
#                                 " by interpolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TEVAL ")
#     ) %>%
#     mutate_cond(condition = !is.na(teval)&ptime==teval&timevar>ptime&!is.na(lagdv),
#
#                 conc.teval=interpol(c1=lagdv, c2=depvar, t1=teval, t2=lagtime, t3=timevar, method=method), #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.teval=lagdv+((depvar-lagdv)*(teval-lagtime)/(timevar-lagtime)),      #interpolate
#                 time.teval=teval,
#                 teval.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",teval,
#                                 " ) corrected to ",round(conc.teval,3),
#                                 " by interpolation (sample taken too late)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TEVAL ")
#     ) %>%
#     mutate_cond(condition = !is.na(teval)&ptime==teval&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
#                 conc.teval=depvar*exp(-1*lambda_z*(teval-timevar)),                       # extrapolate
#                 time.teval=teval,
#                 teval.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",teval,
#                                 " ) corrected to ",round(conc.teval,3),
#                                 " by extrapolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TEVAL ")
#     ) %>%
#     # PARTIAL
#     #   TSTART
#     mutate(conc.part=depvar,time.part=timevar) %>%
#     mutate_cond(condition = !is.na(tstart)&ptime==tstart&timevar<ptime&!is.na(leaddv),
#
#                 conc.part=interpol(c1=depvar, c2=leaddv, t1=tstart, t2=timevar, t3=leadtime, method=method), #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.part=depvar+((leaddv-depvar)*(tstart-timevar)/(leadtime-timevar)),  #interpolate
#                 time.part=tstart,
#                 tstart.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tstart,
#                                 " ) corrected to ",round(conc.part,3),
#                                 " by interpolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TSTART ")
#     ) %>%
#     mutate_cond(condition = !is.na(tstart)&ptime==tstart&timevar>ptime&!is.na(lagdv),
#
#                 conc.part=interpol(c1=lagdv, c2=depvar, t1=tstart, t2=lagtime, t3=timevar, method=method), #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.part=lagdv+((depvar-lagdv)*(tstart-lagtime)/(timevar-lagtime)),      #interpolate
#                 time.part=tstart,
#                 tstart.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tstart,
#                                 " ) corrected to ",round(conc.part,3),
#                                 " by interpolation (sample taken too late)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TSTART ")
#     ) %>%
#     mutate_cond(condition = !is.na(tstart)&ptime==tstart&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
#                 conc.part=depvar*exp(-1*lambda_z*(tstart-timevar)),                       # extrapolate
#                 time.part=tstart,
#                 tstart.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tstart,
#                                 " ) corrected to ",round(conc.part,3),
#                                 " by extrapolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TSTART ")
#     ) %>%
#     #   TEND
#     mutate_cond(condition = !is.na(tend)&ptime==tend&timevar<ptime&!is.na(leaddv),
#
#                 conc.part=interpol(c1=depvar, c2=leaddv, t1=tend, t2=timevar, t3=leadtime, method=method), #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.part=depvar+((leaddv-depvar)*(tend-timevar)/(leadtime-timevar)),    #interpolate
#                 time.part=tend,
#                 tend.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tend,
#                                 " ) corrected to ",round(conc.part,3),
#                                 " by interpolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TEND ")
#     ) %>%
#     mutate_cond(condition = !is.na(tend)&ptime==tend&timevar>ptime&!is.na(lagdv),
#
#                 conc.part=interpol(c1=lagdv, c2=depvar, t1=tend, t2=lagtime, t3=timevar, method=method),   #interpolate
#                 #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                 #conc.part=lagdv+((depvar-lagdv)*(tend-lagtime)/(timevar-lagtime)),        #interpolate
#                 time.part=tend,
#                 tend.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-2","MDT-2"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tend,
#                                 "  ) corrected to ",round(conc.part,3),
#                                 " by interpolation (sample taken too late)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TEND ")
#     ) %>%
#     mutate_cond(condition = !is.na(tend)&ptime==tend&timevar<ptime&is.na(leaddv)&!is.na(lambda_z),
#                 conc.part=depvar*exp(-1*lambda_z*(tend-timevar)),                         # extrapolate
#                 time.part=tend,
#                 tend.flag=1,
#                 trule.nr=ifelse(tolower(reg)=="sd","SDT-3","MDT-3"),
#                 trule.txt=paste("Concentration ( ",depvar," ) at deviating time ( t=",timevar," instead of ",tend,
#                                 " ) corrected to ",round(conc.part,3),
#                                 " by extrapolation (sample taken too early)",sep=""),
#                 applies.to.time=paste(applies.to.time,"TEND ")
#     ) %>%
#     # T=0            (also needed for AUCpartial since this value may be used for substitution at t=TAU in case TEND=TAU)
#     mutate(conc.lastall=depvar,time.lastall=timevar,newdepvar=NA) %>%
#     mutate_cond(condition = ptime==0&timevar!=ptime&tolower(reg)=="sd",
#                 t0.flag=1,
#                 time.tau=0,time.part=0,time.teval=0,time.lastall=0,                       # set predose time to 0
#                 conc.tau=depvar,conc.part=depvar,conc.teval=depvar,conc.lastall=depvar,
#                 trule.nr="SDT-1",
#                 trule.txt=paste("Deviating time at predose ( t=",timevar," ) set to 0",sep=""),
#                 applies.to.time="PREDOSE"
#     ) %>%
#
#
#     mutate_cond(condition = ptime==0&timevar<ptime&!is.na(depvar)&depvar>0&tolower(reg)=="md"&!is.na(lambda_z),
#                 t0.flag=1,
#                 time.tau=0,time.part=0,time.teval=0,time.lastall=0,
#                 newdepvar=depvar*exp(-1*lambda_z*(0-timevar)),                            # extrapolate
#                 conc.tau=newdepvar,conc.part=newdepvar,conc.teval=newdepvar,conc.lastall=newdepvar,
#                 trule.nr="MDT-3",
#                 trule.txt=paste("Concentration at predose ( ",depvar," ), taken before dosing ( t=",timevar,
#                                 " ) corrected to ",round(newdepvar,3),
#                                 " by extrapolation (sample taken too early)",sep=""),
#                 applies.to.time="PREDOSE"
#     ) %>%
#
#     mutate_cond(condition = ptime==0&timevar<ptime&(is.na(depvar)|depvar==0)&tolower(reg)=="md",
#                 t0.flag=1,
#                 time.tau=0,time.part=0,time.teval=0,time.lastall=0,
#                 newdepvar=depvar,
#                 conc.tau=newdepvar,conc.part=newdepvar,conc.teval=newdepvar,conc.lastall=newdepvar,
#                 trule.nr="MDT-3a",
#                 trule.txt="Time of LOQ or NA value at predose taken before dosing set to 0",
#                 applies.to.time="PREDOSE"
#     ) %>%
#
#     mutate_cond(condition = ptime==0&timevar>ptime&tolower(reg)=="md",
#                 t0.flag=1,
#                 time.tau=0,time.part=0,time.teval=0,time.lastall=0,
#                 conc.tau=NA,conc.part=NA,conc.teval=NA,conc.lastall=NA,                   # set conc to NA
#                 trule.nr="MDT-1",
#                 trule.txt=paste("Concentration at predose ( ",depvar," ), taken after dosing ( t=",timevar,
#                                 " ) set to NA (sample taken too late)",sep=""),
#                 applies.to.time="PREDOSE"
#     ) #%>%
#   #       select(-lambda_z,-leaddv,-lagdv,-leadtime,-lagtime, -missflag, -misstime, -diff, -newdepvar)
#
#   result=data_in
#
#   names(result)[names(result)=="ptime"]=nomtimevar
#   names(result)[names(result)=="timevar"]=timevar   # to be sure time value of created time points are copied to original time variable
#   names(result)[names(result)=="depvar"]=depvar     # to be sure conc value of created time points are copied to original conc variable
#   return(result)
#
# }
#
# # CORRECT.CONC: correct missing concentration at critical time points (e.g, predose, TAU, start and end of user selected AUC interval)
# #
# # Input dataset:
# #
# # - dataset corrected for time deviations, created using the correct.time function
# # - dataset containing lambda.z, created using the est.thalf function
# #
# # USAGE:
# #
# # correct.conc(x,nomtimevar="ntad",timevar="time",depvar="dv",tau=NA,teval=NA,th=NA,reg="sd",method=1,route="po")
# #
# # ARGUMENTS:
# #
# # x         : input dataset name (if called within dplyr: .)
# # nomtimevar: variable name containing the nominal sampling time
# # tau       : dosing interval (for multiple dosing), if single dose, leave empty
# # teval     : user selected AUC interval, if not requested, leave empty
# # tstart    : start time of partial AUC (start>0), if not requested, leave empty
# # tend      : emd time of partial AUC, if not requested, leave empty
# # th        : file name of file with lamdba_z information for each curve
# # reg       : regimen, "sd" or "md"
# # ss        : is steady state reached (y/n)
# # method    : method of interpolation:
# #             1: linear up - linear down
# #             2: linear up - logarithmic down
# # route     : route of drug administration ("po","iv")
# #
# # METHOD:
# #
# # If there is a measurable concentration BEFORE and AFTER the missing concentration, use interpolation
# # If there is NO measurable concentration AFTER the missing concentration, use extrapolation
# # Set missing concentration at predose to 0 (SD, non-endogenous) or value at t=TAU (steady state only)
# # Set missing concentration at t=TAU to value at t=0 (steady state only)
# #
# # OUTPUT:
# #
# # dataset with missing concentrations imputed.
# # The following variables are added:
# #
# #   - crule.nr         : correction rule number
# #   - crule.txt        : text explaining what was altered
# #   - applies.to.conc  : lists all AUCS to which the concentration correction rule applies
# #
# # The following Concentration Deviation Correction Rules will be applied to critical time points (t=0, tau, tstart, tend, teval),
# # if needed:
# #
# # Rule nr.   Regimen    Description                                                      Applied to
# #
# # SDC-1      sd         Set concentration to 0 (only non-endogenous compounds)           t=0
# # SDC-2      sd         impute missing concentration by interpolation                    t=tau,tstart,tend,teval
# # SDC-3      sd         impute missing concentration by extrapolation                    t=tau,tend,teval
# # SDC-4      sd (IV)    impute missing concentration by back-extrapolation               t=0
# #
# # MDC-1      md         impute missing concentration by existing conc at t=0 or t=tau*   t=0,tau
# # MDC-2      md         impute missing concentration by interpolation                    t=tau,tstart,tend,teval
# # MDC-3      md         impute missing concentration by extrapolation                    t=tau,tend,teval
# # MDC-4      md (IV)    impute missing concentration by back-extrapolation               t=0
# #
# # * only if steady state has been reached
#
# correct.conc <- function(x,nomtimevar="ntad",tau=NA,tstart=NA,tend=NA,teval=NA,th=NA,reg="sd",ss="n",route="po",method=1) {
#
#   data_in=x
#
#   if (!missing(th)) { data_in=left_join(data_in,th%>%select(-no.points,-intercept,-r.squared,-adj.r.squared,-thalf)) }
#
#   data_in=data_in %>%
#     mutate(ptime=x[[nomtimevar]],                 # nominal time                            (internal)
#            crule.nr="",                           # correction rule number
#            crule.txt="",                          # explanation of concentration substitution
#            applies.to.conc="",                     # lists all AUCS to which the concentration correction rule applies
#            lambda_z=ifelse("lambda_z"%in%names(.),lambda_z,NA)
#     )
#
#   # create lead and lag variables for each AUC
#
#   #TAU
#   if (!is.na(tau)) {
#     data_in=lag.lead(data_in,nomtimevar1="ptime",depvar1="conc.tau",timevar1="time.tau",
#                      lagc="lag.ctau",lagt="lag.ttau",leadc="lead.ctau",leadt="lead.ttau")
#   }
#   #PARTIAL
#   if (!is.na(tstart)&!is.na(tend)) {
#     data_in=lag.lead(data_in,nomtimevar1="ptime",depvar1="conc.part",timevar1="time.part",
#                      lagc="lag.cpart",lagt="lag.tpart",leadc="lead.cpart",leadt="lead.tpart")
#   }
#   #TEVAL
#   if (!is.na(teval)) {
#     data_in=lag.lead(data_in,nomtimevar1="ptime",depvar1="conc.teval",timevar1="time.teval",
#                      lagc="lag.cteval",lagt="lag.tteval",leadc="lead.cteval",leadt="lead.tteval")
#   }
#
#
#   # Start concentration substitutions
#
#   # Create some variables
#
#   data_in=data_in %>% mutate(t0val=ifelse(is.element(0,ptime),conc.tau[ptime==0],NA),
#                              tauval=ifelse(is.element(tau,ptime),conc.tau[ptime==tau],NA)
#   )
#
#   #TAU
#   if (!is.na(tau)) {
#     data_in=data_in %>% mutate_cond(condition = ptime==tau&is.na(conc.tau)&!is.na(t0val)
#                                     &tolower(reg)=="md"&tolower(ss)=="y",
#                                     conc.tau=t0val,                                        # take value at t=0
#                                     time.tau=tau,
#                                     tau.flag=1,
#                                     crule.nr="MDC-1",
#                                     crule.txt=paste("Missing concentration at t=",tau," corrected to ",round(conc.tau,3),
#                                                     " by imputation with existing concentration at t=0",sep=""),
#                                     applies.to.conc=paste(applies.to.conc,"TAU ")
#     ) %>%
#       mutate_cond(condition = ptime==tau&is.na(conc.tau)&!is.na(lag.ctau)&!is.na(lead.ctau),
#
#                   conc.tau=interpol(c1=lag.ctau, c2=lead.ctau, t1=tau, t2=lag.ttau, t3=lead.ttau, method=method), #interpolate
#                   #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                   #conc.tau=lag.ctau +
#                   #         (lead.ctau-lag.ctau)*(tau-lag.ttau)/(lead.ttau-lag.ttau),  # interpolate
#                   time.tau=tau,
#                   tau.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
#                   crule.txt=paste("Missing concentration at t=",tau," corrected to ",round(conc.tau,3),
#                                   " by interpolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TAU ")
#       ) %>%
#       mutate_cond(condition = ptime==tau&is.na(conc.tau)&!is.na(lag.ctau)
#                   &is.na(lead.ctau)&!is.na(lambda_z),                                                                         conc.tau=lag.ctau*exp(-1*lambda_z*(tau-lag.ttau)),                  # extrapolate
#                   conc.tau=lag.ctau*exp(-1*lambda_z*(tau-lag.ttau)),  # extrapolate
#                   time.tau=tau,
#                   tau.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
#                   crule.txt=paste("Missing concentration at t=",tau," corrected to ",round(conc.tau,3),
#                                   " by extrapolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TAU ")
#       ) %>%
#       mutate(tauval=conc.tau[time.tau==tau])
#   }
#   #TSTART and TEND
#   #TSTART
#   if (!is.na(tstart)) {
#     data_in=data_in %>%  mutate_cond(condition = ptime==tstart&!is.na(tau)&tstart==tau&is.na(conc.part)
#                                      &!is.na(t0val)&tolower(reg)=="md"&tolower(ss)=="y",
#                                      conc.part=t0val,                              # take value at t=0 if TSTART=TAU
#                                      time.part=tstart,
#                                      tstart.flag=1,
#                                      crule.nr="MDC-1",
#                                      crule.txt=paste("Missing concentration at t=",tstart,
#                                                      " corrected to ",round(conc.part,3),
#                                                      " by imputation with existing concentration at t=0",sep=""),
#                                      applies.to.conc=paste(applies.to.conc,"TSTART ")
#     ) %>%
#       mutate_cond(condition = ptime==tstart&is.na(conc.part)&!is.na(lag.cpart)&!is.na(lead.cpart),
#                   conc.part=interpol(c1=lag.cpart, c2=lead.cpart, t1=tstart, t2=lag.tpart, t3=lead.tpart, method=method), #interpolate
#                   #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                   #conc.part=lag.cpart +                                         # interpolate
#                   #          (lead.cpart-lag.cpart)*(tstart-lag.tpart)/(lead.tpart-lag.tpart),
#                   time.part=tstart,
#                   tstart.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
#                   crule.txt=paste("Missing concentration at t=",tstart,
#                                   " corrected to ",round(conc.part,3),
#                                   " by interpolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TSTART ")
#       ) %>%
#       mutate_cond(condition = ptime==tstart&is.na(conc.part)&!is.na(lag.cpart)
#                   &is.na(lead.cpart)&!is.na(lambda_z),
#                   conc.part=lag.cpart*exp(-1*lambda_z*(tstart-lag.tpart)),       # extrapolate
#                   time.part=tstart,
#                   tstart.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
#                   crule.txt=paste("Missing concentration at t=",tstart,
#                                   " corrected to ",round(conc.part,3),
#                                   " by extrapolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TSTART ")
#       )
#   }
#   #TEND
#   if (!is.na(tend)) {
#     data_in=data_in %>%  mutate_cond(condition = ptime==tend&!is.na(tau)&tend==tau&is.na(conc.part)&!is.na(t0val)&tolower(reg)=="md"&tolower(ss)=="y",
#                                      conc.part=t0val,                              # take value at t=0 if TEND=TAU
#                                      time.part=tend,
#                                      tend.flag=1,
#                                      crule.nr="MDC-1",
#                                      crule.txt=paste("Missing concentration at t=",tend," corrected to ",
#                                                      round(conc.part,3),
#                                                      " by imputation with existing concentration at t=0",sep=""),
#                                      applies.to.conc=paste(applies.to.conc,"TEND ")
#     ) %>%
#       mutate_cond(condition = ptime==tend&is.na(conc.part)&!is.na(lag.cpart)&!is.na(lead.cpart),
#                   conc.part=interpol(c1=lag.cpart, c2=lead.cpart, t1=tend, t2=lag.tpart, t3=lead.tpart, method=method), #interpolate
#                   #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                   #conc.part=lag.cpart +                                    # interpolate
#                   #          (lead.cpart-lag.cpart)*(tend-lag.tpart)/(lead.tpart-lag.tpart),
#                   time.part=tend,
#                   tend.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
#                   crule.txt=paste("Missing concentration at t=",tend,
#                                   " corrected to ",round(conc.part,3),
#                                   " by interpolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TEND ")
#       ) %>%
#       mutate_cond(condition = ptime==tend&is.na(conc.part)&!is.na(lag.cpart)
#                   &is.na(lead.cpart)&!is.na(lambda_z),
#                   conc.part=lag.cpart*exp(-1*lambda_z*(tend-lag.tpart)),   # extrapolate
#                   time.part=tend,
#                   tend.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
#                   crule.txt=paste("Missing concentration at t=",tend,
#                                   " corrected to ",round(conc.part,3),
#                                   " by extrapolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TEND ")
#       )
#   }
#   #TEVAL
#   if (!is.na(teval)) {
#     data_in=data_in %>%  mutate_cond(condition = ptime==teval&!is.na(tau)&teval==tau&is.na(conc.teval)&!is.na(t0val)&tolower(reg)=="md"&tolower(ss)=="y",
#                                      conc.teval=t0val,                         # take value at t=0 if TEVAL=TAU
#                                      time.teval=teval,
#                                      teval.flag=1,
#                                      crule.nr="MDC-1",
#                                      crule.txt=paste("Missing concentration at t=",teval,
#                                                      " corrected to ",round(conc.part,3),
#                                                      " by imputation with existing concentration at t=0",sep=""),
#                                      applies.to.conc=paste(applies.to.conc,"TEVAL ")
#     ) %>%
#       mutate_cond(condition = ptime==teval&is.na(conc.teval)&!is.na(lag.cteval)&!is.na(lead.cteval),
#                   conc.teval=interpol(c1=lag.cteval, c2=lead.cteval, t1=teval, t2=lag.tteval, t3=lead.tteval, method=method), #interpolate
#                   #c1+((c2-c1)*(t1-t2)/(t3-t2))
#                   #conc.teval=lag.cteval +                             # interpolate
#                   #(lead.cteval-lag.cteval)*(teval-lag.tteval)/(lead.tteval-lag.tteval),
#                   time.teval=teval,
#                   teval.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-2","MDC-2"),
#                   crule.txt=paste("Missing concentration at t=",teval,
#                                   " corrected to ",round(conc.teval,3),
#                                   " by interpolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TEVAL ")
#       ) %>%
#       mutate_cond(condition = ptime==teval&is.na(conc.teval)&!is.na(lag.cteval)
#                   &is.na(lead.cteval)&!is.na(lambda_z),
#                   conc.teval=lag.cteval*exp(-1*lambda_z*(teval-lag.tteval)), # extrapolate
#                   time.teval=teval,
#                   teval.flag=1,
#                   crule.nr=ifelse(tolower(reg)=="sd","SDC-3","MDC-3"),
#                   crule.txt=paste("Missing concentration at t=",teval,
#                                   " corrected to ",round(conc.teval,3),
#                                   " by extrapolation",sep=""),
#                   applies.to.conc=paste(applies.to.conc,"TEVAL ")
#       )
#   }
#
#   # correct t=0 conc for all aucs where t=0 is needed (for PARTIAL this is NOT needed as it does not start at t=0)
#
#   data_in=data_in %>%
#
#     # back-extrapolate t=0 concentration for all AUCs if route is IV
#     mutate(back_extrap=0,
#            lc1=lead(conc.lastall,1),lc2=lead(conc.lastall,2),
#            lt1=lead(time.lastall,1),lt2=lead(time.lastall,2),
#            firstmeasc=conc.lastall[which(conc.lastall>0)][1],
#            firstmeast=ptime[which(conc.lastall>0)][1]) %>%
#     # if there are NAs or LOQs between t=0 and first measurable conc, set these equal to first measurable conc
#     mutate_cond(condition=ptime>0&ptime<firstmeast&(is.na(conc.lastall)|conc.lastall==0),
#                 conc.lastall=firstmeasc,conc.tau=firstmeasc,conc.teval=firstmeasc) %>%
#     mutate_cond(condition=tolower(route)=="iv"&ptime==0,
#                 back_extrap=ifelse(!is.na(lc1)&!is.na(lc2)&lc1>0&lc2>0&
#                                      lc1>lc2,1,0),
#                 conc.lastall=ifelse(back_extrap==1,exp(log(lc1)+
#                                                          (0-lt1)/(lt2-lt1)*
#                                                          (log(lc2)-log(lc1))),
#                                     firstmeasc),
#                 t0.flag=1,
#                 crule.nr=ifelse(tolower(reg)=="sd","SDC-4","MDC-4"),
#                 crule.txt=paste("Concentration at t=0 back-extrapolated to ",round(conc.lastall,3)," (IV bolus)",sep=""),
#                 time.tau=0, time.teval=0,
#                 conc.tau=conc.lastall, conc.teval=conc.lastall
#     ) %>%
#
#     # ALL and LAST
#     mutate_cond(condition = ptime==0&(is.na(conc.lastall)|conc.lastall>0)&tolower(reg)=="sd"&tolower(route)!="iv",
#                 conc.lastall=0,                                   # set conc to 0
#                 time.lastall=0,
#                 t0.flag=1,
#                 crule.nr="SDC-1",
#                 crule.txt=paste("Missing or measurable concentration at (SD) PREDOSE set to 0",sep="")
#     ) %>%
#     mutate_cond(condition = ptime==0&is.na(conc.lastall)&!is.na(tauval)&tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="iv",
#                 conc.lastall=tauval,                              # take value at t=tau
#                 time.lastall=0,
#                 t0.flag=1,
#                 crule.nr="MDC-1",
#                 crule.txt=paste("Missing concentration at PREDOSE corrected to ",round(conc.tau,3),
#                                 " by imputation with existing concentration at t=TAU",sep=""),
#                 applies.to.conc=paste("PREDOSE")
#     ) %>%
#     #TAU
#     mutate_cond(condition = ptime==0&(is.na(conc.tau)|conc.tau>0)&tolower(reg)=="sd"&tolower(route)!="iv",
#                 conc.tau=0,                                       # set conc to 0
#                 time.tau=0,
#                 t0.flag=1,
#                 crule.nr="SDC-1",
#                 crule.txt=paste("Missing or measurable concentration at (SD) PREDOSE set to 0",sep=""),
#                 applies.to.conc=paste("PREDOSE")
#     ) %>%
#     mutate_cond(condition = ptime==0&is.na(conc.tau)&!is.na(tauval)&tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="iv",
#                 conc.tau=tauval,                                  # take value at t=tau
#                 time.tau=0,
#                 t0.flag=1,
#                 crule.nr="MDC-1",
#                 crule.txt=paste("Missing concentration at PREDOSE corrected to ",round(conc.tau,3),
#                                 " by imputation with existing concentration at t=TAU",sep=""),
#                 applies.to.conc=paste("PREDOSE")
#     ) %>%
#     #TEVAL
#     mutate_cond(condition = ptime==0&(is.na(conc.teval)|conc.teval>0)&tolower(reg)=="sd"&tolower(route)!="iv",
#                 conc.teval=0,                                      # set conc to 0
#                 time.teval=0,
#                 t0.flag=1,
#                 crule.nr="SDC-1",
#                 crule.txt=paste("Missing or measurable concentration at (SD) PREDOSE set to 0",sep=""),
#                 applies.to.conc=paste("PREDOSE")
#     ) %>%
#     mutate_cond(condition = ptime==0&is.na(conc.teval)&!is.na(tauval)&tolower(reg)=="md"&tolower(ss)=="y"&tolower(route)!="iv",
#                 conc.teval=tauval,                                 # take value at t=tau
#                 time.teval=0,
#                 t0.flag=1,
#                 crule.nr="MDC-1",
#                 crule.txt=paste("Missing concentration at PREDOSE corrected to ",round(conc.teval,3),
#                                 " by imputation with existing concentration at t=TAU",sep=""),
#                 applies.to.conc=paste("PREDOSE")
#     )
#   #                                 %>%
#   #                       select(-lag.ctau,-lag.ttau,-lead.ctau,-lead.ttau,
#   #                              -lag.cteval,-lag.tteval,-lead.cteval,-lead.tteval,
#   #                              -lag.cpart,-lag.tpart,-lead.cpart,-lead.tpart,
#   #                              -lambda_z,-ptime, -t0val, -tauval)
#
#   result=data_in
#   return(result)
# }
#
# # TAB.CORR: tabulate for each subject what records were added, time deviations and concentration imputations were applied
# #
# # Input dataset:
# # - result concentration dataset created by the correct.time and correct.conc functions, containing time and conc corrected data
# #
# # USAGE:
# #
# # tab.corr(x=pc,nomtimevar="time",by="subject")
# #
# # ARGUMENTS:
# #
# # x         : concentration dataset with time and conc corrected data (if called within dplyr: .)
# # nomtimevar: variable containing the nominal time
# # by        : by-variable(s), e.g. c("subject","day")
# #
# # METHOD:
# #
# # -
# #
# # OUTPUT:
# #
# # dataset with applied corrections (rule number and rule text) listed by by-variable(s) and nominal time
# #
#
# tab.corr <- function(x,nomtimevar="time",by="subject") {
#   create=x%>%dplyr::rename(rule.nr=create.nr,rule.txt=create.txt)
#   trules=x%>%dplyr::rename(rule.nr=trule.nr,rule.txt=trule.txt,applies.to=applies.to.time)
#   crules=x%>%dplyr::rename(rule.nr=crule.nr,rule.txt=crule.txt,applies.to=applies.to.conc)
#
#   result=bind_rows(create,trules,crules) %>%
#     filter(rule.nr!="") %>%
#     select(one_of(by,nomtimevar),applies.to,rule.nr,rule.txt) %>%
#     arrange_(by,nomtimevar)
#
#   return(result)
# }
#
# # LAG.LEAD: estimate lagging and leading time point and concentration for each time point, used by
# #           interpolation and extrapolation rules
# #
# # Input dataset:
# # - called from the correct.time and correct.conc functions
# #
# # USAGE:
# #
# # lag.lead(x,nomtimevar="ntad",depvar="dv",timevar="tad",lagc="lagdv",lagt="lagtime",leadc="leaddv",leadt="leadtime")
# #
# # ARGUMENTS:
# #
# # x          : concentration dataset with time and conc corrected data (if called within dplyr: .)
# # nomtimevar1: variable containing the nominal time
# # depvar1    : variable name containing the dependent variable (e.g., concentration)
# # timevar 1  : variable name containing the actual sampling time
# # lagc       : variable name that will contain the lagging concentration value in the output dataset
# # lagt       : variable name that will contain the lagging time value in the output dataset
# # leadc      : variable name that will contain the leading concentration value in the output dataset
# # leadt      : variable name that will contain the leading concentration value in the output dataset
# #
# # METHOD:
# #
# # -
# #
# # OUTPUT:
# #
# # dataset with lagging and leading time points and concentration added in separate variables
# #
#
# lag.lead <- function(x,nomtimevar1=NA,depvar1=NA,timevar1=NA,lagc=NA,lagt=NA,leadc=NA,leadt=NA) {
#
#   original=x %>% mutate(depvar=x[[depvar1]],                    # dependent variable                      (internal)
#                         timevar=x[[timevar1]],                  # actual time variable                    (internal)
#                         ptime=x[[nomtimevar1]]                  # nominal time                            (internal)
#   ) %>%
#     mutate(flag=ifelse(!is.na(depvar),0,1)) %>%   # flags type of missing value (in between or at the end)
#     mutate(flag=ifelse(is.na(depvar)&timevar>last(timevar[!is.na(depvar)]),2,flag))
#
#   #1 delete NA's
#
#   no.na=original %>%filter(!is.na(depvar))
#
#   #2 calc lead and lag
#
#   no.na=no.na %>% arrange(ptime) %>%
#     mutate(leadc=lead(depvar),              # concentration at next sampling time     (internal)
#            lagc=lag(depvar),                # concentration at previous sampling time (internal)
#            leadt=lead(timevar),             # next sampling time                      (internal)
#            lagt=lag(timevar)                # previous sampling time                  (internal)
#     ) %>%
#     select(ptime,leadc,lagc,leadt,lagt)
#
#   #3 merge with original
#
#   newdata=left_join(original,no.na,by="ptime")
#
#   newdata=newdata %>% arrange(ptime) %>%
#     mutate(leadc  =ifelse(flag==1,locf(leadc),leadc),
#            leadt  =ifelse(flag==1,locf(leadt),leadt),
#            lagc   =ifelse(flag==2,last(depvar[!is.na(depvar)]),lagc),
#            lagt   =ifelse(flag==2,last(timevar[!is.na(depvar)]),lagt)
#     ) %>%
#     arrange(-ptime) %>%
#     mutate(lagc   =ifelse(flag==1,locf(lagc),lagc),
#            lagt   =ifelse(flag==1,locf(lagt),lagt)
#     ) %>%
#     arrange(ptime) %>%
#     mutate(leadc  =ifelse(ptime==last(ptime),NA,leadc),
#            leadt  =ifelse(ptime==last(ptime),NA,leadt)
#     )
#
#   names(newdata)[names(newdata)=="lagc"]=lagc
#   names(newdata)[names(newdata)=="lagt"]=lagt
#   names(newdata)[names(newdata)=="leadc"]=leadc
#   names(newdata)[names(newdata)=="leadt"]=leadt
#
#   return(newdata)
#
# }
#
# # INTERPOL: interpolation of time deviations and missing concentrations
# #
# # used for interpolation of concentrations at critical time points if the sample
# # was taken too early or too late, or if concentration values at critical time points are missing
# #
# # if sample taken too early: corrconc = conc +((leadc-conc)*(desttime-time)/(leadt-time))
# # if sample taken too late : corrconc = lagc +((conc-lagc) *(desttime-lagt)/(time-lagt))
# # if sample concentration is missing: corrconc = lagc +((leadc-lagc) *(desttime-lagt)/(leadt-lagt))
# #
# # where: conc=actual conc, time=actual time measured/recorded
# #        desttime=critical time point (TAU, TEVAL, TSTART, TEND)
# #        corrconc=corrected concentration at the critical timepoint
# # USAGE:
# #
# # interpol(c1=NA, c2=NA, t1=NA, t2=NA, t3=NA, method=)
# #
# # ARGUMENTS:
# #
# # c1 - t3    : concentration and time values used in interpolation: c1+((c2-c1)*(t1-t2)/(t3-t2))
# # method     : 1: all linear interpolation
# #              2: lin up / log down
# #
# # remark: if c1 or c2 is 0, always linear interpolation will be performed
#
# interpol = function(c1=NA, c2=NA, t1=NA, t2=NA, t3=NA, method=1) {
#   ifelse(method==2&c1>c2&c1>0&c2>0,
#          exp(log(c1)+((log(c2)-log(c1))*(t1-t2)/(t3-t2))), # log down
#          c1+((c2-c1)*(t1-t2)/(t3-t2)))                     # lin or c1|c2==0
#
# }
#
# # TRAP: trapezoidal rules used in AUC calculations
# #
# # USAGE:
# #
# # trap(x=, y=, method=)
# #
# # ARGUMENTS:
# #
# # x, y    : variable names of x and y coordinates
# # method  : 1: all linear trapedoizal rule
# #           2: linear trap. rule up / logarithmic trap. rule down
# #           3: linear trap. rule before first Tmax, logarithmic trap. rule after first Tmax
# #
# # remark: if c1 or c2 is 0, always linear interpolation will be performed
#
# trap=function(x=NA, y=NA, method=1) {
#   cm=max(df1$dv,na.rm=T)
#   tmax=first(df1$tad[df1$dv==cm&!is.na(df1$dv)])
#   if (method==1)   {
#     z=  sum((x-lag(x))*(y+lag(y))/2,na.rm=T)
#   }
#   if (method==2) {
#     z=  sum(ifelse(lag(y)>y&lag(y)>0&y>0,
#                    (x-lag(x))*(y-lag(y))/log(y/lag(y)),
#                    (x-lag(x))*(y+lag(y))/2
#     ),na.rm=T)
#   }
#   if (method==3) {
#     z=  sum(ifelse(x>tmax&lag(y)>0&y>0&lag(y)!=y,
#                    (x-lag(x))*(y-lag(y))/log(y/lag(y)),
#                    (x-lag(x))*(y+lag(y))/2
#     ),na.rm=T)
#   }
#
#   return(z)
# }
#
# ##############################################
#
# ###### end of code
#
#
#
#
