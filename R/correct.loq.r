#' Impute Concentrations Below the Limit of Quantitation
#'
#' Imputes LOQ values according to the chosen LOQ substitution rule.
#' 
#' Imputations will be applied to the original depvar (no new concentration
#' variable will be created).
#' @param x input dataset name contains all uncorrected data, including LOQ
#' @param by column names in x indicating grouping variables
#' @param nomtimevar variable name containing the nominal sampling time after dose
#' @param timevar variable name containing the actual sampling time after dose
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param loqvar variable name containing the LOQ value
#' @param loqrule rule number to be applied to the LOQ values in the curve. x$loqrule overrides if provided
#' * 1: 0 before first measurable concentration (FMC); NA after FMC
#' * 2: 0 before FMC; 0 after FMC
#' * 3: 0 before FMC; 0.5xLOQ for first consecutive LOQ after FMC, NA for other LOQ
#' * 4: 0 before FMC; 0.5xLOQ for first consecutive LOQ after FMC, 0 for other LOQ
#' @return A dataset with imputed BLOQ concentrations using the chosen imputation rule
#' @export
#' @importFrom dplyr lag
#' @examples
#' library(magrittr)
#' library(dplyr)
#' library(qpNCA)
#' x <- Theoph
#' ntad <- c(0,0.25,0.5,1,2,4,5,7,9,12,24)
#' for(i in 1:nrow(x)){
#'   time  <- x$Time[[i]]
#'   delta <- abs(ntad - time)
#'   best  <- min(delta)
#'   index <- match(best, delta)
#'   nom   <- ntad[[index]]
#'   x$ntad[[i]] <- nom
#' }
#' rm(list = c('time','delta','best','index','nom', 'i','ntad'))
#' x %<>% rename(time = Time, dv = conc)
#' x %<>% mutate(bloq = ifelse(dv==0,1,0), loq = 0.01, tad = time, loqrule = 1, 
#'               subject=as.numeric(Subject), ntad=as.numeric(ntad))
#' x %>% head
#' x %<>% correct.loq('subject')
#' x %>%  head
#'
correct.loq <- function(
  x,
  by=character(0),
  nomtimevar="ntad",
  timevar="time",
  depvar="dv",
  bloqvar="bloq",
  loqvar="loq",
  loqrule=1
){
  supplied <- character(0)
  if(!missing(loqrule)) supplied <- 'loqrule'
  x <- group_by_at(x, vars(by))
  x <- do(
    .data = x,
    .correct.loq(
      .,
      nomtimevar = nomtimevar,
      timevar = timevar,
      depvar = depvar,
      bloqvar = bloqvar,
      loqvar = loqvar,
      loqrule = loqrule,
      supplied = supplied
    )
  )
  x <- ungroup(x)
  x
}

.correct.loq <- function(
    x,
    nomtimevar,
    timevar,
    depvar,
    bloqvar,
    loqvar,
    loqrule,
    supplied
  ){

if('loqrule' %in% names(x)){
  if('loqrule' %in% supplied){
    warning('loqrule supplied as column overrides like-named argument')
  }
  loqrule <- unique(x$loqrule)
  # x$loqrule <- NULL
}
if(length(loqrule) > 1) {
  warning('loqrule has length > 1; only first value will be used')
  loqrule <- loqrule[[1]]
}

data_in = x

  data_in = data_in %>%
    mutate(depvar1=x[[depvar]],                    # dependent variable                      (internal)
           timevar1=x[[timevar]],                  # actual time variable                    (internal)
           ptime1=x[[nomtimevar]],                 # nominal time                            (internal)
           bloqvar1=x[[bloqvar]],                  # indicated LOQ value (1: yes, 0: no)     (internal)
           loqvar1=x[[loqvar]],                    # LOQ value (concentration)               (internal)
           loqrule.nr="",                          # correction rule number
           loqrule.txt=""                          # explanation of time correction
    )

  data_in = data_in %>%
    mutate(firstmeast=timevar1[which(depvar1>0)][1],
           consecutive=ifelse(bloqvar1==1&lag(bloqvar1)==1,1,0),
           anymeas=ifelse(any(bloqvar1==0,na.rm=T),1,0)
    )

if (any(data_in$anymeas==1)) {  
  if (loqrule==1) {

    data_in = data_in %>%
      mutate_cond(
        condition=(bloqvar1==1&timevar1<firstmeast),
        depvar1=0,
        loqrule.nr="LOQ1",
        loqrule.txt="BLOQ values before first measurable concentration set to 0"
      ) %>%
      mutate_cond(
        condition=(bloqvar1==1&timevar1>firstmeast),
        depvar1=NA,
        loqrule.nr="LOQ1",
        loqrule.txt="BLOQ values after first measurable concentration set to missing"
      )
  }

  if (loqrule==2) {

    data_in = data_in %>%
      mutate_cond(condition=(bloqvar1==1&timevar1<firstmeast),depvar1=0,loqrule.nr="LOQ2",
                  loqrule.txt="BLOQ values before first measurable concentration set to 0") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast),depvar1=0,loqrule.nr="LOQ2",
                  loqrule.txt="BLOQ values after first measurable concentration set to 0")

  }

  if (loqrule==3) {

    data_in = data_in %>%
      mutate_cond(condition=(bloqvar1==1&timevar1<firstmeast),depvar1=0,loqrule.nr="LOQ3",
                  loqrule.txt="BLOQ values before first measurable concentration set to 0") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==0),depvar1=loqvar1/2,bloqvar1=0,
                  loqrule.nr="LOQ3",
                  loqrule.txt="First BLOQ value after first measurable concentration set to 1/2*LOQ") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==1),depvar1=NA,loqrule.nr="LOQ3",
                  loqrule.txt="Consecutive BLOQ values after first measurable concentration set to missing")

  }

  if (loqrule==4) {

    data_in = data_in %>%
      mutate_cond(condition=(bloqvar1==1&timevar1<firstmeast),depvar1=0,loqrule.nr="LOQ4",
                  loqrule.txt="BLOQ values before first measurable concentration set to 0") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==0),depvar1=loqvar1/2,bloqvar1=0,
                  loqrule.nr="LOQ4",
                  loqrule.txt="First BLOQ value after first measurable concentration set to 1/2*LOQ") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==1),depvar1=0,loqrule.nr="LOQ4",
                  loqrule.txt="Consecutive BLOQ values after first measurable concentration set to 0")

  }

} 

  result = data_in %>%
    select(-depvar,-nomtimevar,-timevar)

  names(result)[names(result)=="ptime1"]=nomtimevar
  names(result)[names(result)=="timevar1"]=timevar
  names(result)[names(result)=="depvar1"]=depvar

  return(result)

}
