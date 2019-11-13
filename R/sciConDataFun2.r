#' Similar to ConDataFun2 from qpToolkit, but prints very small (<0.0001) or very large (>9999) numbers in scientific notation.
#'
#' @param y
#' @param nSignif
#' @param na.rm
#'
#' @export
sciConDataFun2 = function (y, nSignif, na.rm) {
  if(median(y, na.rm=TRUE) < 0.0001 | median(y, na.rm = TRUE) > 9999){
    paste(formatC(signif(median(y, na.rm=na.rm), nSignif), format="e", digits=2),
          " (", formatC(signif(min(y, na.rm=na.rm), nSignif),  format="e", digits=2),
          " - ", formatC(signif(max(y, na.rm=na.rm), nSignif),  format="e", digits=2),
          ")", sep = "")
  } else paste(signif(median(y, na.rm=na.rm), nSignif), " (", signif(min(y, na.rm=na.rm), nSignif),
               " - ", signif(max(y, na.rm=na.rm), nSignif), ")", sep = "")}
