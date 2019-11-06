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
    mutate_all(funs(make.true.NA)) %>%
    dplyr::select(Subject, everything(), -Statistics,-!!subjVar)

  return(out)
}

#Listing.Sumstat.Param(subset(data_parent,formulation=="A"), subjVar="ID", by=NULL, nsig=3, vars_ignore=ignore_var,keep=c("adj.r.squared"),adj.r.squared>=0.7)
#out1=Listing.Param(ds, subjVar='ID', nsig=3)
#out2 = Sumstat.Param(ds, subjVar='ID', by=NULL, nsig=3, vars_ignore=ignore_var, keep=c("adj.r.squared"),adj.r.squared>=0.7)

