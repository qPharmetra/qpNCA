#' Calculates lambda_z and thalf for each PK curve (defined using group_by)
#' @importFrom Hmisc Cs
#' @description The function starts with the last three sample points and performs log-linear regression on it. It then adds one sampling point at a time (including and ending at tmax) and performs the regression again. Internal variable EST checks whether there are at least 3 timepoints for estimation after removal of LOQs and NAs
#' @param x a dataset (not needed to be corrected time and conc)
#' @param timevar variable name containing the sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param includeCmax include results of regression including Cmax in selection? (y/n)
#' @return a dataset with estimates for each regression analysis in one observation. The following parameters are available:
#' no.points = number of data points used in the regression analysis; \cr
#' intercept = estimated intercept;\cr
#' lambda_z = -1*estimated slope;\cr
#' r.squared = square of the correlation coefficient;\cr
#' adj.r.squared = adjusted square of the correlation coefficient;\cr
#' thalf = elimination half-life;\cr
#' start_th = time of first sample included in the thalf estimation;\cr
#' end_th = time of last sample included in the thalf estimation;\cr
#' includeCmax = include results of regression including Cmax in selection? (y/n)\cr
#' @examples
#' th = Theoph %>%
#'  group_by(Subject=as.numeric(Subject)) %>%
#'  do(est.thalf(.,timevar="Time",depvar="conc",includeCmax="Y")) %>%
#'  ungroup()
#' @export
est.thalf <- function(x,timevar="time",depvar="dv",includeCmax="Y"){
  data_in = x %>% mutate(timevar=x[[timevar]],
                         depvar=x[[depvar]]) %>%
    filter(!is.na(depvar)&depvar>0) %>%
    filter(timevar>=first(timevar[depvar==max(depvar,na.rm=T)]))      # include Cmax

  if (tolower(includeCmax)=="n") {
    data_in=data_in %>% filter(timevar>first(timevar[depvar==max(depvar,na.rm=T)]))   # exclude Cmax
  }

  est=0
  i=length(data_in$timevar)-2
  result = data.frame(matrix(ncol=7,nrow = 1))
  if (i>=1) {
    est=1
    result = data.frame(matrix(ncol=7,nrow = i))
  }

  while (i>=1) {
    ## make subset of data frame for each number of data points
    xx = data_in[i:nrow(data_in), ]
    ## execute loglinear model fit
    lmcall=lm(log(xx$depvar)~xx$timevar)
    lmcall.estimates = as.numeric(summary(lmcall)$coef[,"Estimate"])
    ## save results of loglin lm fit
    result[i,1]=length(data_in$timevar)-i+1
    result[i,2]=lmcall.estimates[1]
    result[i,3]=lmcall.estimates[2]*-1
    result[i,4]=summary(lmcall)$r.squared
    result[i,5]=summary(lmcall)$adj.r.squared
    result[i,6]=data_in$timevar[i]
    result[i,7]=last(data_in$timevar)
    i=i-1
  }

  names(result) = Cs(no.points,intercept,lambda_z,r.squared,adj.r.squared,start_th,end_th)
  if (est==1) {
    result=result %>% mutate(sel=no.points[adj.r.squared==max(adj.r.squared)],
                             thalf=log(2)/lambda_z,
                             includeCmax=includeCmax
    ) %>%
      filter(sel==no.points) %>%
      select(-sel)
  }
  else {
    result=result %>% mutate(no.points=NA,intercept=NA,lambda_z=NA,r.squared=NA,adj.r.squared=NA,thalf=NA,
                             start_th=NA,end_th=NA)
  }
  return(result)
}
