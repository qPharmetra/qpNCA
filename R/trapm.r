#' Trapezoidal rules used in AUMC calculations
#' trapm(x=, y=, method=)
#' remark: if c1 or c2 is 0, always linear interpolation will be performed
#' @param x variable names of x coordinates
#' @param y variable names of y coordinates
#' @param method
#'           1: all linear trapedoizal rule\cr
#'           2: linear trap. rule up / logarithmic trap. rule down\cr
#'           3: linear trap. rule before first Tmax, logarithmic trap. rule after first Tmax \cr
#' @export
trapm=function(x=NA, y=NA, method=1) {
  cm=max(y,na.rm=T)
  tmax=first(x[y==cm&!is.na(y)])
  if (method==1)   {
    z=  sum((x-lag(x))*(y*x+lag(y)*lag(x))/2,na.rm=T)
  }
  if (method==2) {
    z=  sum(ifelse(lag(y)>y&lag(y)>0&y>0,
                   (x-lag(x))*(y*x-lag(y)*lag(x))/log(y/lag(y))-(x-lag(x))^2*(y-lag(y))/(log(y/lag(y)))^2,
                   (x-lag(x))*(y*x+lag(y)*lag(x))/2
    ),na.rm=T)
  }
  if (method==3) {
    z=  sum(ifelse(x>tmax&lag(y)>0&y>0&lag(y)!=y,
                   (x-lag(x))*(y*x-lag(y)*lag(x))/log(y/lag(y))-(x-lag(x))^2*(y-lag(y))/(log(y/lag(y)))^2,
                   (x-lag(x))*(y*x+lag(y)*lag(x))/2
    ),na.rm=T)
  }

  return(z)
}
