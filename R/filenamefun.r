#' Create file name for regression plots (*.png) from by-variables (in plot.reg function)
#' This function is used by qpNCA wrapper function.
#' @param df A dataframe
#' @param by A grouping variable
#' @export
filenamefun <- function(df,by) {
  filename=""
  for (i in 1:length(by)) {
    filename <- paste0(filename,by[i],"_",df[[by[i]]],"_")
  }
  filename=substr(filename,1,(nchar(filename)-1))
  return(filename)
}
