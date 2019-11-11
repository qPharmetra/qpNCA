#' Format the table NCA parameters table
#' Base R way of NCA tabulation functions (see vignettes for more help)
#' @param ds a resulting data frame from nca.summary.table
#' @param units a data frame containing units
#' @export
Format.Table = function(ds, units=parUnits){
  ncaTab = ds
  names(ncaTab) = paste(units$dispName[match(tolower(names(ncaTab)), units$param)],
                        units$unit[match(tolower(names(ncaTab)), units$param)],
                        sep="")
  names(ncaTab) = ifelse(names(ncaTab)=="NANA", "~", names(ncaTab))
  return(ncaTab)
}