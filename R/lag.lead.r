#' Estimate lagging and leading time points and concentrations for each time point
#'
#' @param x
#' @param nomtimevar1
#' @param depvar1
#' @param timevar1
#' @param lagc
#' @param lagt
#' @param leadc
#' @param leadt
#'
#' @return
#' @export
#'
#' @examples
lag.lead <- function(x,nomtimevar1=NA,depvar1=NA,timevar1=NA,lagc=NA,lagt=NA,leadc=NA,leadt=NA) {

  original=x %>% mutate(depvar=x[[depvar1]],                    # dependent variable                      (internal)
                        timevar=x[[timevar1]],                  # actual time variable                    (internal)
                        ptime=x[[nomtimevar1]]                  # nominal time                            (internal)
  ) %>%
    mutate(flag=ifelse(!is.na(depvar),0,1)) %>%   # flags type of missing value (in between or at the end)
    mutate(flag=ifelse(is.na(depvar)&timevar>last(timevar[!is.na(depvar)]),2,flag))

  #1 delete NA's

  no.na=original %>%filter(!is.na(depvar))

  #2 calc lead and lag

  no.na=no.na %>% arrange(ptime) %>%
    mutate(leadc=lead(depvar),              # concentration at next sampling time     (internal)
           lagc=lag(depvar),                # concentration at previous sampling time (internal)
           leadt=lead(timevar),             # next sampling time                      (internal)
           lagt=lag(timevar)                # previous sampling time                  (internal)
    ) %>%
    select(ptime,leadc,lagc,leadt,lagt)

  #3 merge with original

  newdata=left_join(original,no.na,by="ptime")

  newdata=newdata %>% arrange(ptime) %>%
    mutate(leadc  =ifelse(flag==1,locf(leadc),leadc),
           leadt  =ifelse(flag==1,locf(leadt),leadt),
           lagc   =ifelse(flag==2,last(depvar[!is.na(depvar)]),lagc),
           lagt   =ifelse(flag==2,last(timevar[!is.na(depvar)]),lagt)
    ) %>%
    arrange(-ptime) %>%
    mutate(lagc   =ifelse(flag==1,locf(lagc),lagc),
           lagt   =ifelse(flag==1,locf(lagt),lagt)
    ) %>%
    arrange(ptime) %>%
    mutate(leadc  =ifelse(ptime==last(ptime),NA,leadc),
           leadt  =ifelse(ptime==last(ptime),NA,leadt)
    )

  names(newdata)[names(newdata)=="lagc"]=lagc
  names(newdata)[names(newdata)=="lagt"]=lagt
  names(newdata)[names(newdata)=="leadc"]=leadc
  names(newdata)[names(newdata)=="leadt"]=leadt

  return(newdata)

}
