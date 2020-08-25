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
      x[[arg]] <- NULL
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
  mutate(
    dosevar=dose,
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
    mrtinf.obs= ifelse(
      tolower(ss)=="n",
      aumcinf.obs/aucinf.obs,
      (aumctau + tau*(aucinf.obs-auctau))/auctau
    ),
    mrtinf.pred= ifelse(
      tolower(ss)=="n",
      aumcinf.pred/aucinf.pred,
      (aumctau + tau*(aucinf.pred-auctau))/auctau
    ),
    vz.f.obs=  ifelse(
      tolower(ss)=="n",
      dosevar*factor/(lambda_z*aucinf.obs),
      dosevar*factor/(lambda_z*auctau)
    ),
    vz.f.pred= ifelse(
      tolower(ss)=="n",
      dosevar*factor/(lambda_z*aucinf.pred),
      dosevar*factor/(lambda_z*auctau)
    ),
    #   vss.obs= mrtinf.obs*cl.f.obs,
    #   vss.pred= mrtinf.pred*cl.f.pred,
    vss.obs= ifelse(route=="ivb"|route=="ivi",mrtinf.obs*cl.f.obs,NA),
    vss.pred= ifelse(route=="ivb"|route=="ivi",mrtinf.pred*cl.f.pred,NA),
    pctextr.obs=(clast.obs/lambda_z)/aucinf.obs*100,
    pctextr.pred=(clast.pred/lambda_z)/aucinf.pred*100,
    pctback.obs=area.back.extr/aucinf.obs*100,
    pctback.pred=area.back.extr/aucinf.pred*100
  ) %>%
  select(-dosevar)

  # rename CL and V after IVB/IVI

  qpiv <- tolower(route)=="ivb"|tolower(route)=="ivi"
  result %<>% mutate( cl.obs =    ifelse( qpiv, cl.f.obs,  NA))
  result %<>% mutate( cl.pred =   ifelse( qpiv, cl.f.pred, NA))
  result %<>% mutate( vz.obs =    ifelse( qpiv, vz.f.obs,  NA))
  result %<>% mutate( vz.pred =   ifelse( qpiv, vz.f.pred, NA))
  result %<>% mutate( cl.f.obs =  ifelse(!qpiv, cl.f.obs,  NA))
  result %<>% mutate( cl.f.pred = ifelse(!qpiv, cl.f.pred, NA))
  result %<>% mutate( vz.f.obs =  ifelse(!qpiv, vz.f.obs,  NA))
  result %<>% mutate( vz.f.pred = ifelse(!qpiv, vz.f.pred, NA))

# rename CL at steady state
  ssiv <- tolower(ss) == 'y' &  qpiv
  ssev <- tolower(ss) == 'y' & !qpiv

  result %<>% mutate(cl.f.ss = NA, cl.ss = NA)

  result %<>% mutate(cl.ss     = ifelse( ssiv, cl.obs,   cl.ss   ))
  result %<>% mutate(cl.obs    = ifelse(!ssiv, cl.obs,   NA      ))
  result %<>% mutate(cl.pred   = ifelse(!ssiv, cl.pred,  NA      ))

  result %<>% mutate(cl.f.ss   = ifelse( ssev, cl.f.obs, cl.f.ss ))
  result %<>% mutate(cl.f.obs  = ifelse(!ssev, cl.f.obs, NA      ))
  result %<>% mutate(cl.f.pred = ifelse(!ssev, cl.f.pred,  NA      ))

  # if (tolower(ss)=="y") {
  #   if (tolower(route)=="ivb"|tolower(route)=="ivi") {
  #     result = result %>% mutate( cl.ss   = cl.obs,   cl.obs   = NA, cl.pred   = NA, cl.f.ss = NA )
  #   }else{
  #     result = result %>% mutate( cl.f.ss = cl.f.obs, cl.f.obs = NA, cl.f.pred = NA, cl.ss = NA )
  #   }
  # }else{
  #   result = result %>% mutate(cl.f.ss=NA, cl.ss=NA)
  # }

  return(result)
}

# sessionInfo()
# R version 3.6.2 (2019-12-12)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
#
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so
#
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] knitr_1.28     qpNCA_1.0.24   tidyr_1.0.2    testthat_2.3.1
# [5] magrittr_1.5   dplyr_0.8.4
#
# loaded via a namespace (and not attached):
#   [1] metrumrg_5.57       splines_3.6.2       foreach_1.4.8
# [4] gam_1.16.1          ellipse_0.4.1       gtools_3.8.1
# [7] Formula_1.2-3       assertthat_0.2.1    latticeExtra_0.6-29
# [10] pillar_1.4.3        backports_1.1.5     lattice_0.20-38
# [13] glue_1.3.1          digest_0.6.25       RColorBrewer_1.1-2
# [16] checkmate_2.0.0     qpToolkit_0.2.7     colorspace_1.4-1
# [19] htmltools_0.4.0     Matrix_1.2-18       plyr_1.8.5
# [22] XML_3.99-0.3        pkgconfig_2.0.3     purrr_0.3.3
# [25] scales_1.1.0        gdata_2.18.0        jpeg_0.1-8.1
# [28] htmlTable_1.13.3    tibble_2.1.3        ggplot2_3.2.1
# [31] withr_2.1.2         nnet_7.3-12         lazyeval_0.2.2
# [34] cli_2.0.1           survival_3.1-8      crayon_1.3.4
# [37] fansi_0.4.1         nlme_3.1-142        MASS_7.3-51.4
# [40] foreign_0.8-72      tools_3.6.2         data.table_1.12.8
# [43] hms_0.5.3           lifecycle_0.1.0     stringr_1.4.0
# [46] xpose4_4.7.0        munsell_0.5.0       cluster_2.1.0
# [49] packrat_0.5.0       compiler_3.6.2      snowfall_1.84-6.1
# [52] rlang_0.4.4         grid_3.6.2          iterators_1.0.12
# [55] rstudioapi_0.11     htmlwidgets_1.5.1   base64enc_0.1-3
# [58] gtable_0.3.0        codetools_0.2-16    deSolve_1.27.1
# [61] reshape_0.8.8       reshape2_1.4.3      R6_2.4.1
# [64] gridExtra_2.3       utf8_1.1.4          Hmisc_4.3-1
# [67] readr_1.3.1         stringi_1.4.6       Rcpp_1.0.3
# [70] vctrs_0.2.3         rpart_4.1-15        acepack_1.4.1
# [73] png_0.1-7           tidyselect_1.0.0    xfun_0.12
