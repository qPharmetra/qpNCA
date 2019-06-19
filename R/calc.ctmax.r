#' Calculates Cmax, Tmax from raw data for each PK curve (define using group_by)
#' @importFrom dplyr arrange "%>%" mutate summarize filter group_by do
#' @description Input dataset can contain all uncorrected data, including LOQ; estimate first occurence of maximum concentration for each PK curve; if all concentrations are NA, sets Cmax and Tmax also to NA
#' @param x a dataset
#' @param timevar a column name
#' @param depvar a column name
#' @return a dataset with estimates for the Cmax (maximum concentration) and Tmax (time of first occurence of cmax) parameters, one observation per subject
#' @import dplyr
#' @export
calc.ctmax <- function(x,timevar="time",depvar="dv") {
              result=x %>% mutate(depvar=x[[depvar]],        # calculated dependent variable           (internal)
                                  timevar=x[[timevar]]       # calculated time variable                (internal)
              ) %>%
                filter(!is.na(depvar))            # na.rm=T in the max function does not work?

              if (dim(result)[1]==0) {            # if all concentrations are NA, set Cmax and Tmax to NA
                result=result%>%summarise(cmax=NA,
                                          tmax=NA)
              }
              else {
                result=result%>% summarise(cmax=max(depvar),
                                           tmax=first(timevar[depvar==cmax]) # there might be more than 1 timepoint with Cmax
              )
              }

              return(result)
}
