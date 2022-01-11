#' Calculate Cmax and Tmax
#'
#' Calculates Cmax and Tmax from raw data for each PK curve defined using \code{by}. \cr
#'
#' Input dataset can contain all uncorrected data, including LOQ;
#' estimate first occurence of maximum concentration for each PK curve;
#' if all concentrations are NA, sets Cmax and Tmax also to NA.
#' @param x data.frame
#' @param by character: column names in x indicating grouping variables; default is as.character(dplyr::groups(x))
#' @param timevar variable name containing the actual sampling time after dose
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @return A dataset with estimates for the Cmax (maximum concentration)
#' and Tmax (time of first occurence of cmax) parameters: one observation per subject
#' @import magrittr
#' @importFrom dplyr arrange mutate summarize filter group_by do summarise first rename
#' @export
#' @examples
#' \donttest{
#' library(magrittr)
#' library(dplyr)
#' data(ncx)
#' x <- ncx
#' x %<>% group_by(subject)
#' x %<>% correct.loq
#' x %>% calc.ctmax %>% head
#' }
calc.ctmax <- function(
  x,
  by = NULL,
  timevar="time",
  depvar="dv"
){
  if(is.null(by)) by <- as.character(groups(x))
  x <- group_by_at(x, vars(by))
  x <- do(
    .data = x,
    .calc.ctmax(
      .,
      timevar = timevar,
      depvar = depvar
    )
  )
 # x <- ungroup(x)
  x
}

.calc.ctmax <- function(x,timevar,depvar){
    x %<>% rename(
    depvar = !!depvar,        # calculated dependent variable           (internal)
    timevar = !!timevar       # calculated time variable                (internal)
  )
  x %<>% filter(!is.na(depvar))            # na.rm=T in the max function does not work?
  if (!nrow(x)) {            # if all concentrations are NA, set Cmax and Tmax to NA
    x %<>% summarise(cmax = NA, tmax = NA)
  } else {
    x %<>% summarise(
      cmax=max(depvar),
      tmax=first(timevar[depvar==cmax]) # there might be more than 1 timepoint with Cmax
    )
  }
  x
}
