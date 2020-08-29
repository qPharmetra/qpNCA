#' Calculate Cmax and Tmax
#'
#' Calculates Cmax and Tmax from raw data for each PK curve defined using \code{by}.
#' Input dataset can contain all uncorrected data, including LOQ;
#' estimate first occurence of maximum concentration for each PK curve;
#' if all concentrations are NA, sets Cmax and Tmax also to NA.
#' @param x data.frame
#' @param by column names in x indicating grouping variables
#' @param timevar column name in x indicating time
#' @param depvar column name in x indicating concentration
#' @return A dataset with estimates for the Cmax (maximum concentration)
#' and Tmax (time of first occurence of cmax) parameters: one observation per subject
#' @import magrittr
#' @importFrom dplyr arrange mutate summarize filter group_by do summarise first
#' @export
#' @examples
#' example(est.thalf)
#' ctmax <- x %>% calc.ctmax(by = 'subject')
#' ctmax %>% head
calc.ctmax <- function(x,by = character(0), timevar="time",depvar="dv"){
  x <- group_by_at(x, vars(by))
  x <- do(
    .data = x,
    .calc.ctmax(
      .,
      timevar = timevar,
      depvar = depvar
    )
  )
  x <- ungroup(x)
  x
}

.calc.ctmax <- function(x,timevar,depvar){
    result <- x %>% mutate(
    depvar = x[[depvar]],        # calculated dependent variable           (internal)
    timevar = x[[timevar]]       # calculated time variable                (internal)
  ) %>%
  filter(!is.na(depvar))            # na.rm=T in the max function does not work?
  if (dim(result)[1]==0) {            # if all concentrations are NA, set Cmax and Tmax to NA
    result %<>% summarise(cmax = NA, tmax = NA)
  } else {
    result %<>% summarise(
      cmax=max(depvar),
      tmax=first(timevar[depvar==cmax]) # there might be more than 1 timepoint with Cmax
    )
  }
  return(result)
}
