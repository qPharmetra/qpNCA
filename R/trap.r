#' Calculate Area Under the Curve Using Trapezoids
#'
#' Calculates AUC using the trapezoidal method. Assumes data represent a single profile.
#' Despite choice of method, only linear interpolation is used
#' for areas of intervals beginning or ending with y: 0.
#'
#' @param x x variable, i.e. time
#' @param y y variable, i.e. concentration
#' @param method method:
#' * 1: linear up - linear down
#' * 2: linear up - logarithmic down
#' * 3: linear before Tmax, logarithmic after Tmax
#' @return area (length-one numeric)

trap <- function(x = NA, y = NA, method = 1){
  stopifnot(length(x) == length(y))
  cm <- max(y, na.rm = T )
  tmax <- first(x[y == cm & !is.na(y) ])
  if (method == 1){
    z <- sum( (x-lag(x)) * (y+lag(y)) / 2, na.rm = T )
  }
  if (method==2) {
    z <- sum(
      ifelse(
        lag(y) > y & lag(y) > 0 & y > 0,
        (x - lag(x)) * (y - lag(y)) / log(y/lag(y)), # See Rowland & Tozer 471
        (x - lag(x)) * (y + lag(y)) / 2
      ),
      na.rm=T
    )
  }
  if (method==3) {
    z <- sum(
      ifelse(

        x > tmax    &
        lag(y) > 0  &
        y > 0       &
        lag(y) != y,

        (x-lag(x))*(y-lag(y))/log(y/lag(y)),
        (x-lag(x))*(y+lag(y))/2
      ),
      na.rm=T
    )
  }

  return(z)
}
