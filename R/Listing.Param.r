#' Generate NCA Parameter Table
#'
#' Generates NCA parameter table.
#'
#' @param ds data.frame
#' @param subjVar column name in ds indicating subject
#' @param nsig number of significant digits
#' @importFrom qpToolkit formatted.signif
#' @import dplyr
#' @export
#' @examples
#' library(dplyr)
#' # Listing_Param_theo = Listing.Param(NCAResults$pkpar, subjVar="ID", nsig=3)
Listing.Param = function(ds, subjVar, nsig=3){
  subjVar = rlang::sym(subjVar)
  paramlist = ds %>%
    mutate_if(is.numeric, qpToolkit::formatted.signif, digits=nsig) %>% #Set significant digits for the numeric columns
    dplyr::select(order(colnames(.))) %>%
    dplyr::select(!!subjVar,everything())
  return(paramlist)
}



