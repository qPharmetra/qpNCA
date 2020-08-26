#' Calculate Area Under the Moment Curve Using Trapezoids
#'
#' Calculates AUMC using the trapezoidal method.
#' Assumes data represent a single profile.
#' Despite choice of method, only linear interpolation is used
#' for areas of intervals beginning or ending with y: 0.

#' @param x variable names of x coordinates
#' @param y variable names of y coordinates
#' @param method method:
#' * 1: linear up - linear down
#' * 2: linear up - logarithmic down
#' * 3: linear before Tmax, logarithmic after Tmax
#' @return area (length-one numeric)
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
