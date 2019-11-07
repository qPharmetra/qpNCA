#' Creates summary statistics of concenctrations by time variable and grouping variable
#'
#' @param ds
#' @param subjVar
#' @param by
#' @param nsig
#' @param vars_ignore
#' @param keep
#' @export
Sumstat.Param = function(ds, subjVar, by=NULL, nsig=3, vars_ignore="", keep="" ,...){

  arg = rlang::enexprs(...)
  subjVar=rlang::sym(subjVar)
  by1= rlang::syms(by)
  by=c(by1)

  # Function to replace all empty cells with NA

  empty_as_na <- function(x){
    ifelse(as.character(x)!="", x, NA)
  }

  # Function to calculate coefficient of variation (CV) percentage

  CV = function(x) (sd(x,na.rm=T)/mean(x,na.rm=T))*100

  #Create dataset without the variables to ignore (vars_ignore)
  dsa = ds %>%
    dplyr::select(-vars_ignore) %>%
    mutate_all(funs(empty_as_na)) %>%
    mutate_all(funs(as.numeric))

  #Create dataset with the variables to ignore (vars_ignore) and with the filtering
  dsb = ds %>%
    dplyr::select(vars_ignore,keep) %>%
    dplyr::filter(!!!arg) %>%
    mutate_all(funs(empty_as_na)) %>%
    mutate_all(funs(as.numeric))

  # Create summary statistics
  # N
  out1a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) sum(!is.na(.)) else NA)) %>%
    mutate(STAT="N") %>%
    dplyr::select(-one_of(keep))

  out1b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) sum(!is.na(.)) else NA))

  out1 = cbind(out1a,out1b)

  # Mean
  out2a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) mean(., na.rm = T) else NA)) %>%
    mutate(STAT="Mean") %>%
    dplyr::select(-one_of(keep))

  out2b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) mean(., na.rm = T) else NA))

  out2 = cbind(out2a,out2b)

  # SD
  out3a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) sd(., na.rm = T) else NA)) %>%
    mutate(STAT="SD") %>%
    dplyr::select(-one_of(keep))

  out3b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) sd(., na.rm = T) else NA))

  out3 = cbind(out3a,out3b)

  # Median
  out4a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) median(., na.rm = T) else NA)) %>%
    mutate(STAT="Median") %>%
    dplyr::select(-one_of(keep))

  out4b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) median(., na.rm = T) else NA))

  out4 = cbind(out4a,out4b)

  # Min
  out5a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) min(., na.rm = T) else NA)) %>%
    mutate(STAT="Min") %>%
    dplyr::select(-one_of(keep))

  out5b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) min(., na.rm = T) else NA))

  out5 = cbind(out5a,out5b)

  # Max
  out6a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) max(., na.rm = T) else NA)) %>%
    mutate(STAT="Max") %>%
    dplyr::select(-one_of(keep))

  out6b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) max(., na.rm = T) else NA))

  out6 = cbind(out6a,out6b)

  LCI95 = function(x){(qt(0.025, length(x)) * (sqrt(var(x, na.rm=T))/sqrt(length(x)))) + mean(x, na.rm=T)}
  UCI95 = function(x){(qt(0.95, length(x)) * (sqrt(var(x, na.rm=T))/sqrt(length(x)))) + mean(x, na.rm=T)}

  # L95CI
  out7a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(LCI95) %>%
    mutate(STAT="L95CI") %>%
    dplyr::select(-one_of(keep))

  out7b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(LCI95)

  out7 = cbind(out7a,out7b)

  # U95CI
  out8a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(UCI95) %>%
    mutate(STAT="U95CI") %>%
    dplyr::select(-one_of(keep))

  out8b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(UCI95)

  out8 = cbind(out8a,out8b)

  # Geometric mean
  out9a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) qpToolkit::geomean(., na.rm = T) else NA)) %>%
    mutate(STAT="GeoMean") %>%
    dplyr::select(-one_of(keep))

  out9b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(funs(if(!is.na(.)) qpToolkit::geomean(., na.rm = T) else NA))

  out9 = cbind(out9a,out9b)

  # Coefficient of variation
  out10a = dsa %>% group_by(!!!by) %>%
    dplyr::summarise_all(CV) %>%
    mutate(STAT="CV") %>%
    dplyr::select(-one_of(keep))

  out10b = dsb %>% group_by(!!!by) %>%
    dplyr::summarise_all(CV)

  out10 = cbind(out10a,out10b)

  # Combine data

  out = rbind(out1,out2,out3,out4,out5,out6,out7,out8,out9,out10)

  out = out %>%
    dplyr::select(STAT, everything(),-!!subjVar) #Move 'STAT' column to the front

  out = out %>%
    mutate_if(is.numeric, qpToolkit::formatted.signif, digits=nsig) #Set significant digits for the numeric columns

  # REPLACE THE SD,L95CI,U95CI,... VALUEs BY NC (NOT CALCULATED) WHEN N<=2 or N<3?

  out = out %>%
    tidyr::gather(statistics, value, -STAT) %>%
    tidyr::spread(STAT, value) %>% # transpose the dataset
    dplyr::select(statistics,N,everything()) %>%
    mutate_cond(!is.na(N)<=2,SD="NC",L95CI="NC",U95CI="NC") %>% #if less than 2 observation then report NC
    tidyr::gather(Statistics, value, -statistics) %>%
    tidyr::spread(statistics, value)

  out$Statistics = factor(out$Statistics,
                          c("N", "Mean", "SD", "Median", "Min", "Max", "CV", "L95CI", "U95CI", "GeoMean"),
                          levels=c("N", "Mean", "SD", "Median", "Min", "Max", "CV", "L95CI", "U95CI", "GeoMean"))


  out = out %>% arrange(!!!by1, Statistics) %>%
    dplyr::select(order(colnames(.))) %>%
    dplyr::select(!!!by1, Statistics, everything())


  # Make true NAs recognized by is.na()

  out = out


  return(out)

}
