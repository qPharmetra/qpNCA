#' Calculate Lambda-Z and T-half
#'
#' Calculates lambda_z and thalf for each PK curve identified using \code{by}.
#' The function starts with the last three sample points and performs
#' log-linear regression on it. It then adds one sampling point at a time
#' (including and ending at tmax) and performs the regression again.
#' Internal variable EST checks whether there are at least 3 timepoints
#' for estimation after removal of LOQs and NAs.
#'
#' @param x a dataset (not needed to be corrected time and conc)
#' @param by column names in x indicating grouping variables
#' @param timevar variable name containing the sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param includeCmax include results of regression including Cmax in selection? (y/n); x$includeCmax overrides if provided
#' @param exclvar a variable name containing information about points to be excluded (these should have exclvar = 1)
#' @importFrom stats lm
#' @importFrom dplyr first last
#' @return
#' a dataset with estimates for each regression analysis in one observation.
#' The following parameters are available.
#' * **no.points** number of data points used in the regression analysis
#' * **intercept** estimated intercept
#' * **lambda_z** -1*estimated slope
#' * **r.squared** square of the correlation coefficient
#' * **adj.r.squared** adjusted square of the correlation coefficient
#' * **thalf** elimination half-life
#' * **start_th** time of first sample included in the thalf estimation
#' * **end_th time** of last sample included in the thalf estimation
#' * **includeCmax** include results of regression including Cmax in selection? (y/n)
#' * **points_excluded** are time points excluded from the half-life estimation? (y/n)
#' @export
#' @examples
#' example(correct.loq)
#' x %<>% mutate(includeCmax = 'Y')
#' th <- x %>% est.thalf('subject')
#' th %>%  head

est.thalf <- function(
  x,
  by = character(0),
  timevar="time",
  depvar="dv",
  includeCmax="Y",
  exclvar=NA){
  supplied <- character(0)
  if(!missing(includeCmax)) supplied <- 'includeCmax'
  x <- group_by_at(x, vars(by))
  x <- do(
    .data = x,
    .est.thalf(
      .,
      timevar = timevar,
      depvar = depvar,
      includeCmax = includeCmax,
      exclvar = exclvar,
      supplied = supplied
    )
  )
  x <- ungroup(x)
  x
}
.est.thalf <- function(
  x,
  timevar,
  depvar,
  includeCmax,
  exclvar,
  supplied
){

  if('includeCmax' %in% names(x)){
    if('includeCmax' %in% supplied){
      warning('includeCmax supplied as column overrides like-named argument')
    }
    includeCmax <- unique(x$includeCmax)
    x$includeCmax <- NULL
  }
  if(length(includeCmax) > 1) {
    warning('includeCmax has length > 1; only first value will be used')
    includeCmax <- includeCmax[[1]]
  }


  data_in = x %>% mutate(timevar=x[[timevar]],
                         depvar=x[[depvar]],
                         exclvar=x[[exclvar]])

  if (!is.na(exclvar) & !(exclvar %in% names(x))) stop(paste("Exclusion variable",exclvar,"does not exist"), call.=F)

  if (!is.na(exclvar)) {

    anyexcl=0
    if (any(data_in$exclvar==1)) {anyexcl=1}

    data_in = data_in %>%
      filter(exclvar!=1|is.na(exclvar))  # remove samples to be excluded from the regression

  }

  data_in = data_in %>%
    filter(!is.na(depvar)&depvar>0) %>%
    filter(timevar>=first(timevar[depvar==max(depvar,na.rm=T)]))      # include Cmax

  if (tolower(includeCmax)=="n") {
    data_in=data_in %>% filter(timevar>first(timevar[depvar==max(depvar,na.rm=T)]))   # exclude Cmax
  }

  est=0
  i=length(data_in$timevar)-2
  result = data.frame(matrix(ncol=7,nrow = 1))

  if (i>=1) {
    est=1
    result = data.frame(matrix(ncol=7,nrow = i))
  }

  while (i>=1) {
    ## make subset of data frame for each number of data points
    xx = data_in[i:nrow(data_in), ]
    ## execute loglinear model fit
    lmcall=lm(log(xx$depvar)~xx$timevar)
    lmcall.estimates = as.numeric(summary(lmcall)$coef[,"Estimate"])
    ## save results of loglin lm fit
    result[i,1]=length(data_in$timevar)-i+1
    result[i,2]=exp(lmcall.estimates[1])  # exponentiate to see the actual intercept
    result[i,3]=lmcall.estimates[2]*-1
    result[i,4]=summary(lmcall)$r.squared
    result[i,5]=summary(lmcall)$adj.r.squared
    result[i,6]=data_in$timevar[i]
    result[i,7]=last(data_in$timevar)
    i=i-1
  }
  names(result) = c('no.points','intercept','lambda_z','r.squared','adj.r.squared','start_th','end_th')
  if (est==1) {
    result=result %>% mutate(sel=no.points[adj.r.squared==max(adj.r.squared)],
                             thalf=log(2)/lambda_z,
                             includeCmax=includeCmax
    ) %>%
      filter(sel==no.points) %>%
      select(-sel)
  }
  else {
    result=result %>% mutate(no.points=NA,intercept=NA,lambda_z=NA,r.squared=NA,adj.r.squared=NA,thalf=NA,
                             start_th=NA,end_th=NA)
  }

  if (!is.na(exclvar)) {

    if (anyexcl==1) {
      result = result %>%
        mutate(points_excluded="Y")
    }
    else {
      result = result %>%
        mutate(points_excluded="N")
    }
  }
  else {
    result = result %>%
      mutate(points_excluded="N")
  }

  return(result)
}
