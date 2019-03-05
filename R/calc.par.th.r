#' Calculates PK parameters that need lambda_z
#'
#' @param x result parameter dataset from calc.par
#' @param th result dataset from est.thalf
#' @param cov covariates dataset (containing at least dose for CL calculation)
#' @param dosevar variable containing the dose amount
#' @param factor conversion factor for CL and V calculation (e.g. dose in mg, conc in ng/mL, factor=1000)
#' @param reg regimen, "sd" or "md"
#' @param ss is steady state reached (y/n)
#'
#' @return dataset with estimates for the following parameters, one observation per subject: \cr
#'   all parameters calculated in th \cr
#'   all parameters calculated in par \cr
#'   clast.pred: predicted concentration at tlast \cr
#'   aucinf.obs: aucinf based on observed concentration at tlast \cr
#'   aucinf.pred: aucinf based on predicted concentration at tlast \cr
#'   aumcinf.obs: area under the first moment curve extrapolated to infinity, based on observed concentration at tlast \cr
#'   aumcinf.pred: area under the first moment curve extrapolated to infinity, based on predicted concentration at tlast \cr
#'   cl.f.obs: clearance based on aucinf.obs, at steady state based on auctau \cr
#'   cl.f.pred: clearance based on aucinf.pred \cr
#'   mrt.obs: Mean residence time based on aumcinf.obs and aucinf.obs \cr
#'   mrt.pred: Mean residence time based on aumcinf.pred and aucinf.pred \cr
#'   vz.f.obs: distribution volume based on cl.f.obs, at steady state based on auctau \cr
#'   vz.f.pred: distribution based on cl.f.pred \cr
#'   vss.obs: Steady-state volume based on cl.obs and mrt.obs \cr
#'   vss.pred: Steady-state volume based on cl.pred and mrt.pred \cr
#'   regimen (reg) \cr
#'   steady state reached Y/N? (ss) \cr
#' NOTE: ctmax must be merged separately as those were calculated from uncorrected data \cr
#'
#' @examples
#' library(dplyr)
#' #'# We need half-lives for corrections, so first let's get that.
#' th = Theoph %>%
#'  group_by(Subject=as.numeric(Subject)) %>%
#'  do(est.thalf(.,timevar="Time",depvar="conc",includeCmax="Y")) %>%
#'  ungroup()
#'
#'# We need nominal time variable as well, so let's generate that.
#' ID <- as.numeric(Theoph$Subject)
#' NTAD <- c(0,0.3,0.5,1,2,4,5,7,9,12,24)
#' Theoph1 <- Theoph %>% mutate(NTAD=metrumrg::snap(Time, NTAD))
#'
#' #let's say we want AUC0-8. We only have 7 and 9 hr concentrations, so we need to interpolate conc for 8 hr.
#' tc = Theoph1 %>%
#' group_by(Subject=as.numeric(Subject)) %>%
#' do(correct.time(.,nomtimevar="NTAD",timevar="Time",depvar="conc",
#'                  tau=,tstart=,tend=,teval=8,th=th,reg="sd")) %>%
#' #above step added timepoints that we will add interpolated concentrations to with this next step
#' do(correct.conc(.,nomtimevar="NTAD",tau=,tstart=,tend=,teval=8,
#' th=th,reg="sd",ss="n")) %>%
#' ungroup()
#' #Now get parameters that do not require LAMBDAZ
#' par <- tc %>%
#'      group_by(Subject=as.numeric(Subject)) %>%
#'      ##tstart by default is 0, so here we will get AUC0-8hr and AUC0-24 as sampling ends at 24 hr
#'      do(calc.par(.,tau=NA, tstart=NA, tend=NA, teval=8, route="po", method=1)) %>%
#'      ungroup()
#' cov <- data.frame(Subject=as.numeric(Theoph1$Subject), DOSE=Theoph$Dose) %>%
#'        distinct(.,.keep_all = T)
#' par <- calc.par.th(x=par,th=th ,cov=cov,
#'                  dose="DOSE",factor=1,
#'                  reg="sd",ss="n")
#' @export
calc.par.th <- function(x=par,th=th,cov=cov,dosevar="dose",factor=1, reg="sd", ss="n") {

  result=left_join(x,th) %>%
    left_join(cov) %>%
    mutate(dosevar=.[[dosevar]],
           reg=tolower(reg),
           ss=tolower(ss),
           clast.pred=exp(intercept-lambda_z*tlast),
           aucinf.obs=auclast+clast.obs/lambda_z,
           aucinf.pred=auclast+clast.pred/lambda_z,
           aumcinf.obs=aumclast+tlast*clast.obs/lambda_z + clast.obs/lambda_z^2,
           aumcinf.pred=aumclast+tlast*clast.pred/lambda_z + clast.pred/lambda_z^2,
           cl.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/aucinf.obs, dosevar*factor/auctau),
           cl.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/aucinf.pred, dosevar*factor/auctau),
           mrt.obs= ifelse(tolower(ss)=="n", aumcinf.obs/aucinf.obs,
                           (aumctau + tau*(aucinf.obs-auctau))/auctau),
           mrt.pred= ifelse(tolower(ss)=="n", aumcinf.pred/aucinf.pred,
                            (aumctau + tau*(aucinf.pred-auctau))/auctau),
           vz.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.obs),
                             dosevar*factor/(lambda_z*auctau)),
           vz.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.pred),NA),
           vss.obs= mrt.obs*cl.f.obs,
           vss.pred= mrt.pred*cl.f.pred,
           pctextr.obs=(clast.obs/lambda_z)/aucinf.obs*100,
           pctextr.pred=(clast.pred/lambda_z)/aucinf.pred*100,
           pctback.obs=area.back.extr/aucinf.obs*100,
           pctback.pred=area.back.extr/aucinf.pred*100
    ) %>%
    select(-dosevar)

  return(result)

}
