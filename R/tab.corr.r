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
