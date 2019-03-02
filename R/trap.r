#' Trapezoid function
#'
#' @param x
#' @param y
#' @param method
#'
#' @return
#' @export
#'
#' @examples
trap <- function(df1, x=NA, y=NA, method=1) { #added df1 here by Krina on 1 march 2019 due to error message without it
  cm=max(df1$dv,na.rm=T)
  tmax=first(df1$tad[df1$dv==cm&!is.na(df1$dv)])
  if (method==1)   {
    z=  sum((x-lag(x))*(y+lag(y))/2,na.rm=T)
  }
  if (method==2) {
    z=  sum(ifelse(lag(y)>y&lag(y)>0&y>0,
                   (x-lag(x))*(y-lag(y))/log(y/lag(y)),
                   (x-lag(x))*(y+lag(y))/2
    ),na.rm=T)
  }
  if (method==3) {
    z=  sum(ifelse(x>tmax&lag(y)>0&y>0&lag(y)!=y,
                   (x-lag(x))*(y-lag(y))/log(y/lag(y)),
                   (x-lag(x))*(y+lag(y))/2
    ),na.rm=T)
  }

  return(z)
}
