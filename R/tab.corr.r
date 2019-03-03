#'Tabulates what records were added, time deviations and concentration imputations were applied, for each subject
#' @param x concentration dataset created by the correct.time and correct.conc functions, containing time and conc corrected data
#' @param nomtimevar ariable containing the nominal time
#' @param by by-variable(s), e.g. c("subject","day")
#' @return dataset with applied corrections (rule number and rule text) listed by by-variable(s) and nominal time
#' @examples
#' @export
tab.corr <- function(x,nomtimevar="time",by="subject") {
  create=x%>%dplyr::rename(rule.nr=create.nr,rule.txt=create.txt)
  trules=x%>%dplyr::rename(rule.nr=trule.nr,rule.txt=trule.txt,applies.to=applies.to.time)
  crules=x%>%dplyr::rename(rule.nr=crule.nr,rule.txt=crule.txt,applies.to=applies.to.conc)

  result=bind_rows(create,trules,crules) %>%
    filter(rule.nr!="") %>%
    select(one_of(by,nomtimevar),applies.to,rule.nr,rule.txt) %>%
    arrange_(by,nomtimevar)

  return(result)
}
