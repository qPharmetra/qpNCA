#' Similar to tabSummarize from qpToolkit, but prints very small (<0.0001) or very large (>9999) numbers in scientific notation.
#'
#' @param formula
#' @param data
#' @param nSignif
#' @param extra.blank.line
#' @param ndigits.categorical
#' @param na.rm
#'
#' @return
#' @export
tabSummarize_sci = function (formula, data, nSignif = 3, extra.blank.line = TRUE,
                             ndigits.categorical = 1)
{
  allX = all.vars(nlme::getResponseFormula(formula)[[2]])
  allY = all.vars(nlme::getCovariateFormula(formula)[[2]])
  BY = lapply(1:length(allX), function(x, allX, data) eval(as.name(allX[[x]]),
                                                           data), data = data, allX = allX)
  names(BY) = allX
  YYY = lapply(allY, function(yyy, data) eval(as.name(yyy),
                                              data), data = data)
  names(YYY) = allY
  theData = do.call("rbind", lapply(1:length(YYY), function(z,
                                                            YYY, BY, extra.blank.line, nSignif, ndigits.categorical) {
    stats = tabStats_sci(x = YYY[[z]], BY = BY, parName = names(YYY)[z],
                         nSignif = nSignif, ndigits.categorical = ndigits.categorical)
    if (extra.blank.line == TRUE) {
      EBL = stats[1, ]
      EBL[1, ] = rep("", ncol(stats))
      stats = rbind(EBL, stats)
    }
    return(stats)
  },
  YYY = YYY, BY = BY, extra.blank.line = extra.blank.line,
  nSignif = nSignif, ndigits.categorical = ndigits.categorical))
  row.names(theData) = 1:nrow(theData)
  theData$parameter = full.names(theData$parameter)
  names.order = as.character(unique(eval(as.name(allX[1]),
                                         data)))
  theData = theData[, c("parameter", names.order, "All")]
  ndf = theData[1, ]
  NNN = tapply(YYY[[1]], BY, length)
  NNN = NNN[names.order]
  ndf[1, ] = c("", paste("(N=", c(as.numeric(NNN), length(YYY[[1]])),
                         ")", sep = ""))
  theData = rbind(ndf, theData)
  return(theData)
}
