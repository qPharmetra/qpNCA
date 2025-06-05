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
#' @importFrom dplyr first lag

trap <- function(x = NA, y = NA, method = 1){
  stopifnot(length(x) == length(y))
  cm <- max(y, na.rm = T )
  tmax <- first(x[y == cm & !is.na(y) ])
  if (method == 1){
    z <- sum( (x-lag(x)) * (y+lag(y)) / 2, na.rm = T )
  }
  if (method == 2) {
    dx <- diff(x)
    y2 <- y[-1]
    y1 <- y[-length(y)]
    down <- which(y1 > y2 & y1 > 0 & y2 > 0)
    up   <- setdiff(seq_along(dx), down)
    z <- 0
    if (length(down)) {
      z <- z + sum(dx[down] * (y2[down] - y1[down]) / log(y2[down]/y1[down]), na.rm = TRUE)
    }
    if (length(up)) {
      z <- z + sum(dx[up] * (y2[up] + y1[up]) / 2, na.rm = TRUE)
    }
  }
  if (method == 3) {
    dx <- diff(x)
    y2 <- y[-1]
    y1 <- y[-length(y)]
    logidx <- which(x[-1] > tmax & y1 > 0 & y2 > 0 & y1 != y2)
    linidx <- setdiff(seq_along(dx), logidx)
    z <- 0
    if (length(logidx)) {
      z <- z + sum(dx[logidx] * (y2[logidx] - y1[logidx]) / log(y2[logidx]/y1[logidx]), na.rm = TRUE)
    }
    if (length(linidx)) {
      z <- z + sum(dx[linidx] * (y2[linidx] + y1[linidx]) / 2, na.rm = TRUE)
    }
  }

  return(z)
}
