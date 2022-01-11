#' Tabulate Corrections
#'
#' Tabulates what records were added, time deviations and
#' concentration imputations were applied, for each subject.
#' @param x concentration dataset created by the correct.time and correct.conc functions, containing time and conc corrected data
#' @param by character: column names in x indicating grouping variables; default is as.character(dplyr::groups(x))
#' @param nomtimevar column in x containing the nominal time after dose
#' @return dataset with applied corrections (rule number and rule text) listed by by-variable(s) and nominal time
#' @export
#' @importFrom dplyr bind_rows one_of arrange_at
#' @examples
#' \donttest{
#' library(magrittr)
#' library(dplyr)
#' data(ncx)
#' x <- ncx
#' x %<>% group_by(subject)
#' x %<>% correct.loq
#' x %<>% correct.time
#' x %<>% correct.conc
#' x %>% tab.corr %>% head
#' }
tab.corr <- function(x, by = NULL, nomtimevar = "time") {
  if(is.null(by)) by <- as.character(groups(x))
  loqrules=x%>%dplyr::rename(rule.nr=loqrule.nr,rule.txt=loqrule.txt)
  create=x%>%dplyr::rename(rule.nr=create.nr,rule.txt=create.txt)
  trules=x%>%dplyr::rename(rule.nr=trule.nr,rule.txt=trule.txt,applies.to=applies.to.time)
  crules=x%>%dplyr::rename(rule.nr=crule.nr,rule.txt=crule.txt,applies.to=applies.to.conc)

  result=bind_rows(loqrules,create,trules,crules) %>%
         filter(rule.nr!="") %>%
         select(one_of(by,nomtimevar),applies.to,rule.nr,rule.txt) %>%
         arrange_at(c(by,nomtimevar))

  return(result)
}
