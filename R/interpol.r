interpol <- function(c1=NA, c2=NA, t1=NA, t2=NA, t3=NA, method=1) {
  ifelse(method==2&c1>c2&c1>0&c2>0,
         exp(log(c1)+((log(c2)-log(c1))*(t1-t2)/(t3-t2))), # log down
         c1+((c2-c1)*(t1-t2)/(t3-t2)))                     # lin or c1|c2==0

}
