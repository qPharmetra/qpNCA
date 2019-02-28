# Calculate parameters

calc.par <- function(x,tau=NA,tstart=NA,tend=NA,teval=NA,route="po",method=1){

# check for each AUC if start and end time/conc are available
  tlast.ok=0         # internal variable
  t0.ok=0            # internal variable
  tau.ok=0           # internal variable
  teval.ok=0         # internal variable
  part.ok=0          # internal variable

  par = x %>% filter(!is.na(time.lastall))     # if critical time point, it's not NA anymore, so all time=NA are deleted

  if (max(par$conc.lastall,na.rm=T)>0) { tlast.ok=1 }

  if (is.element(0,par$ptime)) { t0.ok=ifelse(par$time.lastall[par$ptime==0]==0
                                              &!is.na(par$conc.lastall[par$ptime==0]),1,t0.ok) }

  if (!is.na(tau)&t0.ok==1&is.element(tau,par$time.tau))       { tau.ok= ifelse(is.na(par$conc.tau[par$time.tau==0]) |
                                                                                  is.na(par$conc.tau[par$time.tau==tau]),0,1)  }
  if (!is.na(teval)&t0.ok==1&is.element(teval,par$time.teval)) { teval.ok=
    ifelse(is.na(par$conc.teval[par$time.teval==0]) |
             is.na(par$conc.teval[par$time.teval==teval]),0,1)}
  if (!is.na(tstart)&!is.na(tend)&is.element(tstart,par$time.part)&is.element(tend,par$time.part))
  {  part.ok=  ifelse(is.na(par$conc.part[par$time.part==tstart]) |
                        is.na(par$conc.part[par$time.part==tend]),0,1) }

  par=par %>%    summarise( route=route,
                            method=method,
                            tlast=ifelse(tlast.ok==1,last(time.lastall[!is.na(conc.lastall)&conc.lastall>0]),NA),
                            clast.obs=ifelse(tlast.ok==1,conc.lastall[time.lastall==tlast],NA),
                            tlast.ok=tlast.ok,
                            t0.ok=t0.ok,

                            aucall=ifelse(t0.ok==1,
                                          trap(x=time.lastall[!is.na(conc.lastall)],
                                               y=conc.lastall[!is.na(conc.lastall)], method=method),NA),

                            auclast=ifelse(t0.ok==1&tlast.ok==1,
                                           trap(x=time.lastall[!is.na(conc.lastall)&time.lastall<=tlast],
                                                y=conc.lastall[!is.na(conc.lastall)&time.lastall<=tlast], method=method),NA),

                            aumcall=ifelse(t0.ok==1&tlast.ok==1,
                                           trap(x=time.lastall[!is.na(conc.lastall)],
                                                y=time.lastall[!is.na(conc.lastall)] *
                                                  conc.lastall[!is.na(conc.lastall)], method=method),NA),

                            aumclast=ifelse(t0.ok==1&tlast.ok==1,
                                            trap(x=time.lastall[!is.na(conc.lastall)&time.lastall<=tlast],
                                                 y=conc.lastall[!is.na(conc.lastall)&time.lastall<=tlast]*
                                                   time.lastall[!is.na(conc.lastall)&time.lastall<=tlast], method=method),NA),

                            calc.tau=tau.ok,
                            auctau=ifelse(t0.ok==1&tau.ok==1,
                                          trap(x=time.tau[!is.na(conc.tau)&time.tau<=tau],
                                               y=conc.tau[!is.na(conc.tau)&time.tau<=tau], method=method),NA),
                            aumctau=ifelse(t0.ok==1&tau.ok==1,
                                           trap(x=time.tau[!is.na(conc.tau)&time.tau<=tau],
                                                y=conc.tau[!is.na(conc.tau)&time.tau<=tau]*
                                                  time.tau[!is.na(conc.tau)&time.tau<=tau], method=method),NA),
                            tau=ifelse(!is.na(tau),tau,NA),

                            calc.teval=teval.ok,
                            aucteval=ifelse(t0.ok==1&teval.ok==1,
                                            trap(x=time.teval[!is.na(conc.teval)&time.teval<=teval],
                                                 y=conc.teval[!is.na(conc.teval)&time.teval<=teval], method=method),NA),
                            teval=ifelse(!is.na(teval),teval,NA),

                            calc.part=part.ok,
                            aucpart=ifelse(part.ok==1,
                                           trap(x=time.part[!is.na(conc.part)&time.part>=tstart&time.part<=tend],
                                                y=conc.part[!is.na(conc.part)&time.part>=tstart&time.part<=tend], method=method),NA),
                            tstart=ifelse(!is.na(tstart),tstart,NA),
                            tend=ifelse(!is.na(tend),tend,NA),
                            c0=ifelse(tolower(route)=="iv",conc.lastall[ptime==0],NA),


                            area.back.extr=ifelse(tolower(route)=="iv",
                                                  trap(x=time.lastall[time.lastall<=firstmeast&!is.na(conc.lastall)],
                                                       y=conc.lastall[time.lastall<=firstmeast&!is.na(conc.lastall)],
                                                       method=method),NA)
  )

  result=par

  names(result)[names(result)=="aucteval"]=ifelse(!is.na(teval),paste("auc",teval,sep=""),"aucteval") # rename aucteval to aucXX
  names(result)[names(result)=="aucpart"]=ifelse(!is.na(tend),paste("auc",tstart,"_",tend,sep=""),"aucpart") # rename aucpart to aucXX_XX
  if(is.na(tau))   result = subset(result, select = -c(tau,calc.tau,auctau,aumctau) )    # drop au(m)ctau and tau if not requested
  if(is.na(teval)) result = subset(result, select = -c(teval,calc.teval,aucteval) )      # drop aucteval and teval if not requested
  if(is.na(tend))  result = subset(result, select = -c(tstart,tend,calc.part,aucpart) )  # drop aucpart, tstart and tend if not requested

  return(result)
}
