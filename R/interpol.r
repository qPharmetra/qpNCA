#' Interpolate Concentrations
#'
#' Interpolates conentrations.
#' Used by correct.xx functions to interpolate concentrations.
#' Uses linear interpolation unless method is 2 (log down), c1 > c2,
#' and both concentrations are non-zero.
#' @param c1 concentration 1
#' @param c2 concentration 2
#' @param t1 time 1
#' @param t2 time 2
#' @param t3 time 3
#' @param method calculation method (1, 2, or 3)
interpol <- function(c1=NA, c2=NA, t1=NA, t2=NA, t3=NA, method=1) {
  ifelse(
    method == 2 & c1>c2 & c1>0 & c2>0,
    exp(log(c1)+((log(c2)-log(c1))*(t1-t2)/(t3-t2))), # log down
    c1+((c2-c1)*(t1-t2)/(t3-t2))                      # lin or c1|c2==0
  )

}
