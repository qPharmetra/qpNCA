#' Trapezoid function
#'
#' @param x x variable, i.e. time
#' @param y y variable, i.e. concentration
#' @param method (1) linear up - linear down; (2) linear up - logarithmic down; or (3) linear before Tmax, logarithmic after Tmax

#' @export
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
