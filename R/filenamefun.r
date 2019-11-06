#' Create file name for regression plots (*.png) from by-variables (in plot.reg function)
#'
#' @param df
#' @param by
#'
#' @return
#' @export
#'
#' @examples
filenamefun <- function(df,by) {
  filename=""
  for (i in 1:length(by)) {
    filename <- paste0(filename,by[i],"_",df[[by[i]]],"_")
  }
  filename=substr(filename,1,(nchar(filename)-1))
  return(filename)
}
