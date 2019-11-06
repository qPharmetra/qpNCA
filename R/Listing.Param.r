#' Parameters Tables for NCA Report
#'
#' @param ds
#' @param subjVar
#' @param nsig
#' @examples
# library(dplyr)
#
# Listing_Param_theo = Listing.Param(NCAResults$pkpar, subjVar="ID", nsig=3)
#' @export
Listing.Param = function(ds, subjVar, nsig=3){

    subjVar = rlang::sym(subjVar)

  paramlist = ds %>%
    mutate_if(is.numeric, qpToolkit::formatted.signif, digits=nsig) %>% #Set significant digits for the numeric columns
    dplyr::select(order(colnames(.))) %>%
    dplyr::select(!!subjVar,everything())

  return(paramlist)
}



