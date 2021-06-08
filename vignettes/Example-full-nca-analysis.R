## ---- results='hide'----------------------------------------------------------
library(dplyr)
library(qpNCA)
library(knitr)

## -----------------------------------------------------------------------------

mutate_cond <- function (.data, condition, ..., envir = parent.frame()){
  condition <- eval(substitute(condition), .data, envir)
  if(!any(condition))return(.data) # do nothing if nothing to do
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
locf <- function(x){
  good <- !is.na(x)
  positions <- seq(length(x))
  good.positions <- good * positions
  last.good.position <- cummax(good.positions)
  last.good.position[last.good.position == 0] <- NA
  x[last.good.position]
}


## ---- results="markup", warnings=F--------------------------------------------

head(Theoph) %>% kable()

input.data <- Theoph

#we need nominal time variable for some tasks.
ntad <- data.frame(rn=c(1:11),ntad=c(0,0.25,0.5,1,2,4,5,7,9,12,24))

input.data %<>% 
           group_by(Subject) %>%
           mutate(subject=as.numeric(Subject),
                  rn=row_number(),
                  dose=Dose*Wt,
                  bloq=ifelse(conc==0,1,0),
                  loq=0.1,
                  excl_th=0
                  ) %>%
           left_join(ntad) %>%
           ungroup %>%
           arrange(subject,ntad) %>%
           select(subject,ntad,tad=Time,conc,dose,bloq,loq,excl_th)

input.data %<>%
           mutate_cond(condition=subject==2&ntad%in%c(24),conc=NA) %>%
           mutate_cond(condition=subject==4&ntad%in%c(9),conc=NA) %>%
           mutate_cond(condition=subject==3&ntad==9,excl_th=1) %>%
           mutate_cond(condition=subject==6&ntad==24,conc=0,bloq=1) %>%
           filter(!(subject==5&ntad==12))


## ---- results="markup", warnings=F, fig.width = 7-----------------------------

# Create a covariates file, containing at least the dose given

cov = input.data %>%
      distinct(subject,dose)

nca = qpNCA(
      input.data,
      by = "subject",
      nomtimevar = "ntad",
      timevar = "tad",
      depvar = "conc",
      bloqvar = "bloq",
      loqvar = "loq",
      loqrule = 1,
      includeCmax = "Y",
      exclvar = "excl_th",
      plotdir = NA,
      timelab = "Time (h)",
      deplab = "Conc (ng/mL)",
      tau = 24,
      tstart = 4,
      tend = 9,
      teval = 12,
      covariates = cov,
      dose = "dose",
      factor = 1,
      reg = "sd",
      ss = "n",
      route = "EV",
      method = 1
      )


## ---- results="markup", warnings=F, fig.width = 7-----------------------------

# Covariates:

nca$covariates %>% kable()

# Corrections applied:

nca$corrections %>% kable()

# half-life estimation:

nca$half_life %>% kable()

# PK parameters:

nca$pkpar %>% kable()

# Regression plots:

nca$plots


