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
