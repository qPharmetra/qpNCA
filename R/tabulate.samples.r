tabulate.samples = function(data = Mprot.969, analyte = "M.protein", oddCode = Cs(BLOQ,NS,M))
{
  if(any(is.na(data[, analyte]))) {message("column'",analyte,"' contains NA records. Solve this first.");return()}
  resTab = table(data[, analyte])
  resTab = resTab[names(resTab) %in% oddCode]
  missingItems = oddCode %nin% names(resTab)
  if(any(missingItems)){mstuff = rep(0, length(oddCode[missingItems])); names(mstuff) = oddCode[missingItems]; resTab = c(resTab, mstuff)}
  resTab = resTab[match(oddCode, names(resTab))]
  resTab = c(resTab, Numeric.Values = length(data[data[, analyte] %nin% oddCode, analyte]))
  myAll = sum(resTab)
  resTab = c(resTab, All = myAll)
  return(resTab)
}
