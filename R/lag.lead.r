#' Estimate Lagging and Leading Times and Concentrations
#'
#' Estimates lagging and leading times and concentrations.
#' Used by correct.xx functions to estimate lagging and leading timepoints
#' and concentrations for each timepoint.
#' @param x data.frame
#' @param nomtimevar1 column name in x indicating nominal time after dose
#' @param depvar1 column name in x indicating concentration
#' @param timevar1 column name in x indicating actual time after dose
#' @param lagc concentration at previous sampling time
#' @param lagt previous sampling time
#' @param leadc concentration at next sampling time
#' @param leadt next sampling time
#' @param ... ignored
#' @return data.frame
#' @importFrom dplyr last lead lag left_join

lag_lead <- function(
  x,nomtimevar1=NA,depvar1=NA,timevar1=NA,
  lagc=NA,lagt=NA,leadc=NA,leadt=NA,...
){
  # original <- x |>
  #   mutate(depvar = !!depvar1,   # dependent variable   (internal)
  #          timevar = !!timevar1,        # actual time variable (internal)
  #          ptime = !!nomtimevar1       # nominal time         (internal)
  #   ) |>
  #   mutate(flag=ifelse(!is.na(depvar),0,1)) |>   # flags type of missing value (in between or at the end)
  #   mutate(flag=ifelse(is.na(depvar)&timevar>last(timevar[!is.na(depvar)]),2,flag))

  original <- x
  original$depvar <- original[[depvar1]]
  original$timevar <- original[[timevar1]]
  original$ptime <- original[[nomtimevar1]]
  original <- original |>
    mutate(flag=ifelse(!is.na(depvar),0,1)) |>
    mutate(flag=ifelse(is.na(depvar)&timevar>last(timevar[!is.na(depvar)]),2,flag))

  #1 delete NA's

  no.na=original |>filter(!is.na(depvar))

  #2 calc lead and lag

  no.na=no.na |> arrange(ptime) |>
    mutate(leadc=lead(depvar),              # concentration at next sampling time     (internal)
           lagc=lag(depvar),                # concentration at previous sampling time (internal)
           leadt=lead(timevar),             # next sampling time                      (internal)
           lagt=lag(timevar)                # previous sampling time                  (internal)
    ) |>
    select(ptime,leadc,lagc,leadt,lagt)

  #3 merge with original

  newdata=left_join(original,no.na,by="ptime")
  # newdata=left_join(original,no.na)

  newdata = newdata |> arrange(ptime) |>
    mutate(leadc  =ifelse(flag==1,locf(leadc),leadc),
           leadt  =ifelse(flag==1,locf(leadt),leadt),
           lagc   =ifelse(flag==2,last(depvar[!is.na(depvar)]),lagc),
           lagt   =ifelse(flag==2,last(timevar[!is.na(depvar)]),lagt)
    )
  newdata = newdata |> arrange(-ptime)
  newdata = newdata |>
    mutate(lagc   =ifelse(flag==1,locf(lagc),lagc),
           lagt   =ifelse(flag==1,locf(lagt),lagt)
    )
  newdata = newdata |>
    arrange(ptime) |>
    mutate(leadc  =ifelse(ptime==last(ptime),NA,leadc),
           leadt  =ifelse(ptime==last(ptime),NA,leadt)
    )

  names(newdata)[names(newdata)=="lagc"]=lagc
  names(newdata)[names(newdata)=="lagt"]=lagt
  names(newdata)[names(newdata)=="leadc"]=leadc
  names(newdata)[names(newdata)=="leadt"]=leadt

  return(newdata)

}
