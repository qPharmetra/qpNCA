#' Imputes LOQ values according to the chosen LOQ substitution rule
#' Imputations will be applied to the original depvar(no new concentration variable will be created)
#' @param x input dataset name (if called within dplyr: .) contains all uncorrected data, including LOQ
#' @param nomtimevar variable name containing the nominal sampling time
#' @param timevar variable name containing the sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param loqvar variable name containing the LOQ value
#' @param loqrule rule number to be applied to the LOQ values in the curve:
#'
#' @return A dataset with imputed BLOQ concentrations using the chosen imputation rule
#' @export
correct.loq <- function(x,nomtimevar="ntad",timevar="time",depvar="dv",bloqvar="bloq",loqvar="loq",loqrule=1) {
  
data_in=x
  
  data_in = data_in %>%
    mutate(depvar1=x[[depvar]],                    # dependent variable                      (internal)
           timevar1=x[[timevar]],                  # actual time variable                    (internal)
           ptime1=x[[nomtimevar]],                 # nominal time                            (internal)
           bloqvar1=x[[bloqvar]],                  # indicated LOQ value (1: yes, 0: no)     (internal)
           loqvar1=x[[loqvar]],                    # LOQ value (concentration)               (internal)
           loqrule.nr="",                         # correction rule number
           loqrule.txt=""                         # explanation of time correction
    )
  
  data_in = data_in %>%
    mutate(firstmeast=timevar1[which(depvar1>0)][1],
           consecutive=ifelse(bloqvar1==1&lag(bloqvar1)==1,1,0)
    )
  
  if (loqrule==1) {
    
    data_in = data_in %>%
      mutate_cond(condition=(bloqvar1==1&timevar1<firstmeast), depvar1=0, loqrule.nr="LOQ1",
                  loqrule.txt="BLOQ values before first measurable concentration set to 0") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast), depvar1=NA, loqrule.nr="LOQ1",
                  loqrule.txt="BLOQ values after first measurable concentration set to missing")
    
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
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==0),depvar1=loqvar1/2,loqrule.nr="LOQ3",
                  loqrule.txt="First BLOQ value after first measurable concentration set to 1/2*LOQ") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==1),depvar1=NA,loqrule.nr="LOQ3",
                  loqrule.txt="Consecutive BLOQ values after first measurable concentration set to missing")
    
  }
  
  if (loqrule==4) {
    
    data_in = data_in %>%
      mutate_cond(condition=(bloqvar1==1&timevar1<firstmeast),depvar1=0,loqrule.nr="LOQ4",
                  loqrule.txt="BLOQ values before first measurable concentration set to 0") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==0),depvar1=loqvar1/2,loqrule.nr="LOQ4",
                  loqrule.txt="First BLOQ value after first measurable concentration set to 1/2*LOQ") %>%
      mutate_cond(condition=(bloqvar1==1&timevar1>firstmeast&consecutive==1),depvar1=0,loqrule.nr="LOQ4",
                  loqrule.txt="Consecutive BLOQ values after first measurable concentration set to 0")
    
  }
  
  result = data_in %>%
    select(-depvar,-nomtimevar,-timevar) 
  
  names(result)[names(result)=="ptime1"]=nomtimevar
  names(result)[names(result)=="timevar1"]=timevar    
  names(result)[names(result)=="depvar1"]=depvar    
  
  return(result)
  
}
