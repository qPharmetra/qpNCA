#' Creates a table of individual parameterss table with summary statistics
#'
#' @param ds
#' @param subjVar
#' @param by
#' @param nsig
#' @param vars_ignore
#' @param keep
#'
#' @examples
#' @export
Listing.Sumstat.Param = function(ds, subjVar, by, nsig, vars_ignore="", keep="",...){

  #arg = rlang::enexprs(...)
  subjVar=rlang::sym(subjVar)
  by1= rlang::syms(by)
  by=c(by1)

  out1 = Listing.Param(ds, subjVar=subjVar, nsig=nsig)

  #Do not allow summary stats by grouping in combined table for now.
  out2 = Sumstat.Param(ds, subjVar, by=by, nsig=nsig, vars_ignore=vars_ignore, keep=keep)

  out = bind_rows(out1, out2) %>%
    dplyr::select(!!subjVar, Statistics, everything()) %>%
    mutate(Subject = paste(as.character(!!subjVar), Statistics, sep=""),
           Subject = gsub("NA", "", Subject)) %>%
    dplyr::select(Subject, everything()) %>%
    select(-Statistics)

  return(out)
}

