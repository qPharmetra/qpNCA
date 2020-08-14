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
#' @return A dataset with estimates for the following parameters, one observation per subject: \cr
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
#' NOTE: ctmax must be merged separately as those were calculated from uncorrected data.
#'
#' @export
calc.par.th <- function(x=par,by="subject",th=th,covfile=covfile,dose="dose",factor=1, reg="SD", ss="N", route="EV") {
  for(arg in c('reg','ss','factor','route')){
    if(arg %in% names(x)){
      if(!eval(substitute(missing(arg)))){
        warning(arg,' supplied as column overrides like-named argument')
      }
      assign(arg,x[[arg]])
    }
    # if(length(get(arg)) > 1) {
    #   warning(arg, ' has length > 1; only first value will be used')
    #   assign(arg, get(arg)[[1]])
    # }
  }

    if(!missing(covfile)){
      if(is.character(covfile)){
        if(file.exists(covfile)){
          covfile <- read.csv(covfile)
        }else{
          covfile=get(covfile) # to convert the covfile string to the real data frame
        }
      }
    }
                   result=left_join(x,th,by=by) %>%
                      left_join(covfile,by=by) %>%
                      mutate(dosevar=dose,
                             factor=ifelse(is.na(factor),1,factor),
                             reg=tolower(reg),
                             ss=tolower(ss),
                             route=tolower(route),
                             clast.pred=exp(log(intercept)-lambda_z*tlast), # intercept was already exponentiated in calc_thalf
                             aucinf.obs=auclast+clast.obs/lambda_z,
                             aucinf.pred=auclast+clast.pred/lambda_z,
                             aumcinf.obs=aumclast+tlast*clast.obs/lambda_z + clast.obs/lambda_z^2,
                             aumcinf.pred=aumclast+tlast*clast.pred/lambda_z + clast.pred/lambda_z^2,
                             cl.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/aucinf.obs, dosevar*factor/auctau),
                             cl.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/aucinf.pred, dosevar*factor/auctau),
                             mrtinf.obs= ifelse(tolower(ss)=="n", aumcinf.obs/aucinf.obs,
                                             (aumctau + tau*(aucinf.obs-auctau))/auctau),
                             mrtinf.pred= ifelse(tolower(ss)=="n", aumcinf.pred/aucinf.pred,
                                              (aumctau + tau*(aucinf.pred-auctau))/auctau),
                             vz.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.obs),
                                               dosevar*factor/(lambda_z*auctau)),
                             vz.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.pred),
                                               dosevar*factor/(lambda_z*auctau)),
#                             vss.obs= mrtinf.obs*cl.f.obs,
#                             vss.pred= mrtinf.pred*cl.f.pred,
                             vss.obs= ifelse(route=="ivb"|route=="ivi",mrtinf.obs*cl.f.obs,NA),
                             vss.pred= ifelse(route=="ivb"|route=="ivi",mrtinf.pred*cl.f.pred,NA),
                             pctextr.obs=(clast.obs/lambda_z)/aucinf.obs*100,
                             pctextr.pred=(clast.pred/lambda_z)/aucinf.pred*100,
                             pctback.obs=area.back.extr/aucinf.obs*100,
                             pctback.pred=area.back.extr/aucinf.pred*100
                      ) %>%
                      select(-dosevar)

               if (tolower(route)=="ivb"|tolower(route)=="ivi") {    # rename CL and V after IVB/IVI

                 result = result %>% mutate(cl.obs=cl.f.obs, cl.pred=cl.f.pred,
                                            vz.obs=vz.f.obs, vz.pred=vz.f.pred) %>%
                                     mutate(cl.f.obs=NA, cl.f.pred=NA,
                                            vz.f.obs=NA, vz.f.pred=NA)

               }

               else {

                 result = result %>% mutate(cl.obs=NA, cl.pred=NA,
                                            vz.obs=NA, vz.pred=NA)

               }

               if (tolower(ss)=="y") {          #rename CL at steady state

                 if (tolower(route)=="ivb"|tolower(route)=="ivi") {
                 result = result %>% mutate(cl.ss=cl.obs , cl.obs=NA, cl.pred=NA,cl.f.ss=NA)

                 }

                 else {

                 result = result %>% mutate(cl.f.ss=cl.f.obs, cl.f.obs=NA,cl.f.pred=NA,cl.ss=NA)

                 }
               }

               else {

              result = result %>% mutate(cl.f.ss=NA, cl.ss=NA)

               }
               return(result)

}
