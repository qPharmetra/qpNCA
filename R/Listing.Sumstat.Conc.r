#' Creates a table of individual concentrations table with summary statistics by time variable
#'
#' @param ds data.frame
#' @param subjVar column names in ds indicating subject
#' @param timeVar column name in ds indicating time
#' @param dvVar column name in ds indicating concentration
#' @param blqVar column name in ds indicating below limit of quantitation
#' @param by column names in ds the interaction of which indicate data subsets
#' @param carryAlong names of columns in ds to preserve in the result
#' @param nsig number of significant digits
#' @param na.rm whether to remove NA values when computing statistics
#'
#' @examples
#' library(dplyr)
#' NTAD <- c(0,0.3,0.5,1,2,4,5,7,9,12,24)
#' Theoph1 <- Theoph %>%
#'  mutate(NTAD=metrumrg::snap(Time, NTAD)) %>%
#'  mutate(
#'  Subject=as.numeric(as.character(Subject)),
#'  #converting from factor to numeric
#'  BQL = ifelse(conc<=0.25, 1, 0),
#'  #just adding few BLQs to demonstrate functionality
#'  conc= ifelse(conc<=0.25, NA, round(conc, 2))
#'  #just adding digits to demonstrate nsig functionality
#' )
#'
#' test1 <- Listing.Sumstat.Conc(Theoph1, subjVar="Subject", timeVar="NTAD",
#'                               dvVar="conc", blqVar="BQL", carryAlong="Wt", nsig=3,
#'                               na.rm=TRUE)
#'
#' @export
Listing.Sumstat.Conc = function(ds, subjVar, timeVar,
                                dvVar, blqVar, by,
                                carryAlong, nsig=3, na.rm){

  out1 = Listing.Conc(ds, subjVar=subjVar, timeVar=timeVar, dvVar=dvVar, blqVar=blqVar, carryAlong=carryAlong, nsig=nsig)

  #Do not allow summary stats by grouping in combined table for now.
  out2 = Sumstat.Conc(ds, timeVar=timeVar, dvVar=dvVar, by=NULL, nsig=nsig, na.rm=na.rm)  %>%
    select(Statistics, everything())

  subjVar = rlang::sym(subjVar)

  out = bind_rows(out1, out2) %>%
    select(!!subjVar, Statistics, !!carryAlong, everything()) %>%
    mutate(Subject = paste(as.character(!!subjVar), Statistics, sep=""),
           Subject = gsub("NA", "", Subject)) %>%
    select(Subject, everything()) %>%
    select(-Statistics)

  return(out)
}
