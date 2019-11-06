#' Create title for regression plots from by-variables (in plot.reg function)
#'
#' @param df dataset containing concentration-time information of the current curve
#' @param by by-variable(s), e.g. c("subject","day")
#'
#' @return
#' @export
#'
#' @examples
titlefun <- function(df,by) {
  plottitle=""
  for (i in 1:length(by)) {
    plottitle <- paste0(plottitle,by[i],": ",unique(df[[by[i]]])," ")
  }
  plottitle=substr(plottitle,1,(nchar(plottitle)-1))
  return(plottitle)
}
