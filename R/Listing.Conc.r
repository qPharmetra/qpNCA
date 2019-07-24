#' Wide Concentrations vs. Nominal Time Tables for NCA Report
#'
#' @param ds
#' @param subjVar
#' @param timeVar
#' @param concVar
#' @param blqVar
#' @param carryAlong
#' @param nsig
#' @examples
#'library(dplyr)
#'NTAD <- c(0,0.3,0.5,1,2,4,5,7,9,12,24)
#'Theoph1 <- Theoph %>%
#'  mutate(NTAD=metrumrg::snap(Time, NTAD)) %>%
#'  mutate(Subject=as.numeric(as.character(Subject)),     #converting from factor to numeric
#'         BQL = ifelse(conc<=0.25, 1, 0),                #just adding few BLQs to demonstrate functionality
#'         conc= ifelse(conc<=0.25, NA, round(conc, 2)))  #just adding digits to demonstrate nsig functionality
#'
#'test <- Listing.Conc(Theoph1, subjVar="Subject", timeVar="NTAD", dvVar="conc",
#'                     blqVar="BQL", carryAlong="Wt", nsig=3)
#' @export
Listing.Conc = function(ds, subjVar, timeVar, dvVar,
                        blqVar, carryAlong=NULL, nsig=3){

  subjVar = rlang::sym(subjVar)
  timeVar=rlang::sym(timeVar)
  dvVar=rlang::sym(dvVar)
  blqVar=rlang::sym(blqVar)
  carryAlong = rlang::syms(carryAlong)

  conclist = ds %>%
    dplyr::select(!!subjVar, !!timeVar, !!dvVar, !!!carryAlong, !!blqVar) %>%
    dplyr::mutate(Concentrations = qpToolkit::formatted.signif(!!dvVar, nsig)) %>%
    dplyr::mutate(Concentrations = ifelse(!!blqVar == 1, "BLQ", Concentrations)) %>%
    select(-!!dvVar, -!!blqVar) %>%
    do(tibble::rowid_to_column(.)) %>%
    group_by(!!subjVar) %>%
    tidyr::spread(., key=!!timeVar, value="Concentrations", fill=NA)%>%
    ungroup() %>%
    select(-rowid) %>%
    group_by(!!subjVar) %>%
    mutate_all(., .funs=metrumrg::forbak) %>%
    slice(n()) %>%
    ungroup() %>%
    mutate_if(is.character, funs(ifelse(is.na(.), "M", .)))

  return(conclist)
}

