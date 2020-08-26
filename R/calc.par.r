#' Calculate NCA Parameters
#'
#' Calculates AUCs, tlast, clast.obs for each
#' PK curve (define using group_by).
#' @importFrom dplyr arrange mutate summarize filter group_by do select
#' @importFrom tidyr drop_na
#' @import magrittr
#' @param x contains all data after time/concentration deviation corrections obtained from correct.time and correct.conc
#' @param tau  dosing interval (for multiple dosing), if single dose, leave empty
#' @param tstart starting time of user defined interval, if not requested, leave empty
#' @param tend end time of user defined interval, if not requested, leave empty
#' @param teval user selected AUC interval, if not requested, leave empty
#' @param route route of drug administration ("EV","IVB","IVI")
#' @param method method for trapezoidal rule:\cr
#'          1: linear up - linear down \cr
#'          2: linear up - logarithmic down \cr
#'          3: linear before first Tmax, logarithmic after first Tmax \cr
#' @return A dataset with estimates for the following parameters, one observation per subject:
#'\tabular{ll}{
#' t0.ok     \tab flags if t=0 concentration could be corrected/imputes. If not, no AUCs starting at t=0 are calculated \cr
#' tlast.ok  \tab flags if there is at least one measurable concentration. If not, no AUClast can be calculated\cr
#' tlast     \tab time of last sample with measurable concentration\cr
#' clast.obs \tab observed concentration at tlast\cr
#' aucall    \tab auc calculated over all observations, including values below LOQ (which are set to 0)\cr
#' auclast   \tab auc calculated using all observations up to and including the last measurable concentration (clast.obs at tlast)\cr
#' aumcall   \tab aumc calculated over all observations, including values below LOQ (which are set to 0)\cr
#' aumclast  \tab aumc calculated using all observations up to and including the last measurable concentration (clast.obs at tlast)\cr
#' tau       \tab the dosing interval (if specified)\cr
#' calc.tau  \tab flags if AUCtau could be calculated\cr
#' auctau    \tab auc calculated over the dosing interval, only calculated if tau is specified\cr
#' aumctau   \tab aumc calculated over the dosing interval, only calculated if tau is specified\cr
#' teval     \tab user selected AUC interval starting at t=0 (if specified)\cr
#' calc.teval\tab flags if AUCteval could be calculated\cr
#' aucxx     \tab auc calculated from t=0 up to/including teval, only calculated if teval is specified (xx is substituted by teval)\cr
#' calc.part \tab flags if AUCpart could be calculated\cr
#' tstart    \tab start time of partial AUC (if specified)\cr
#' tend      \tab end time of partial AUC (if specified)\cr
#' aucx_y    \tab partial auc from time=x up to/including time=y, where x>0, only calculated if tstart and tend are specified\cr
#' c0        \tab back-extrapolated concentration at t=0 for IV bolus administration\cr
#' area.back.extr \tab area back-extrapolated to 0\cr
#' }
#' @export
#' @importFrom utils read.csv
calc.par <- function(x,tau=NA,tstart=NA,tend=NA,teval=NA,route="EV",method=1){
  for(arg in c('tau','tstart','tend','teval','route','method')){
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

  # Calculate parameters

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

  par=par %>% summarise(
    route=route,
    method = method,
    tlast = ifelse(
      tlast.ok == 1,
      last(time.lastall[!is.na(conc.lastall) & conc.lastall > 0 & bloqvar1 == 0]),
      NA
    ),
    clast.obs = ifelse(tlast.ok == 1, conc.lastall[time.lastall ==
                                                     tlast], NA),
    tlast.ok = tlast.ok,
    t0.ok = t0.ok,

    aucall = ifelse(
      t0.ok == 1,
      trap(
        x = time.lastall[!is.na(conc.lastall)],
        y = conc.lastall[!is.na(conc.lastall)],
        method = method
      ),
      NA
    ),

    auclast = ifelse(t0.ok == 1 & tlast.ok == 1,
                     trap(x = time.lastall[!is.na(conc.lastall) &
                                             time.lastall <= tlast],
                          y = conc.lastall[!is.na(conc.lastall) &
                                             time.lastall <= tlast], method = method),
                     NA),

    aumcall = ifelse(t0.ok == 1 & tlast.ok == 1,
                     trapm(x = time.lastall[!is.na(conc.lastall)],
                           y = conc.lastall[!is.na(conc.lastall)], method =
                             method),
                     NA),

    aumclast = ifelse(t0.ok == 1 &
                        tlast.ok == 1,
                      trapm(x = time.lastall[!is.na(conc.lastall) &
                                               time.lastall <= tlast],
                            y = conc.lastall[!is.na(conc.lastall) &
                                               time.lastall <= tlast], method = method),
                      NA),

    mrtall = ifelse(t0.ok == 1 &
                      tlast.ok == 1, aumcall / aucall, NA),
    mrtlast = ifelse(t0.ok == 1 &
                       tlast.ok == 1, aumclast / auclast, NA),

    calc.tau = tau.ok,

    auctau = ifelse(t0.ok == 1 & tau.ok == 1,
                    trap(x = time.tau[!is.na(conc.tau) &
                                        time.tau <= tau],
                         y = conc.tau[!is.na(conc.tau) &
                                        time.tau <= tau], method = method),
                    NA),
    aumctau = ifelse(t0.ok == 1 & tau.ok == 1,
                     trapm(x = time.tau[!is.na(conc.tau) &
                                          time.tau <= tau],
                           y = conc.tau[!is.na(conc.tau) &
                                          time.tau <= tau], method = method),
                     NA),
    tau = ifelse(!is.na(tau), tau, NA),

    calc.teval = teval.ok,
    aucteval = ifelse(t0.ok == 1 &
                        teval.ok == 1,
                      trap(x = time.teval[!is.na(conc.teval) &
                                            time.teval <= teval],
                           y = conc.teval[!is.na(conc.teval) &
                                            time.teval <= teval], method = method),
                      NA),
    teval = ifelse(!is.na(teval), teval, NA),

    calc.part = part.ok,
    aucpart = ifelse(part.ok == 1,
                     trap(x = time.part[!is.na(conc.part) &
                                          time.part >= tstart & time.part <= tend],
                          y = conc.part[!is.na(conc.part) &
                                          time.part >= tstart & time.part <= tend], method = method),
                     NA),
    tstart = ifelse(!is.na(tstart), tstart, NA),
    tend = ifelse(!is.na(tend), tend, NA),
    c0 = ifelse(tolower(route) == "ivb", conc.lastall[ptime ==
                                                        0], NA),


    area.back.extr = ifelse(tolower(route) == "ivb",
                            trap(x = time.lastall[time.lastall <=
                                                    firstmeast & !is.na(conc.lastall)],
                                 y = conc.lastall[time.lastall <=
                                                    firstmeast & !is.na(conc.lastall)],
                                 method = method),
                            NA)
  )

  result=par

  names(result)[names(result)=="aucteval"]=ifelse(!is.na(teval),paste("auc",teval,sep=""),"aucteval") # rename aucteval to aucXX
  names(result)[names(result)=="aucpart"]=ifelse(!is.na(tend),paste("auc",tstart,"_",tend,sep=""),"aucpart") # rename aucpart to aucXX_XX
  #  if(is.na(tau))   result = subset(result, select = -c(tau,calc.tau,auctau,aumctau) )    # drop au(m)ctau and tau if not requested
  #  if(is.na(teval)) result = subset(result, select = -c(teval,calc.teval,aucteval) )      # drop aucteval and teval if not requested
  #  if(is.na(tend))  result = subset(result, select = -c(tstart,tend,calc.part,aucpart) )  # drop aucpart, tstart and tend if not requested

  return(result)
}
