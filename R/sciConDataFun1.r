#' Similar to ConDataFun1 from qpToolkit, but prints very small (<0.0001) or very large (>9999) numbers in scientific notation.
#'
#' @param y numeric
#' @param nSignif passed to signif()
#' @param na.rm passed to mean
#'
#' @export
sciConDataFun1 = function (y, nSignif, na.rm) {
  if(mean(y, na.rm=TRUE) < 0.0001 | mean(y, na.rm=TRUE) > 9999){
    paste(formatC(signif(mean(y, na.rm=na.rm), nSignif), format="e", digits=2),
          " (", formatC(signif(sqrt(var(y, na.rm=na.rm)), nSignif),  format="e", digits=2),
          ")", sep = "")
  } else  paste(signif(mean(y, na.rm=na.rm), nSignif), " (", signif(sqrt(var(y, na.rm=na.rm)), nSignif), ")", sep = "")
}
