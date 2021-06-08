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


## ---- results="markup", warnings=F--------------------------------------------

loqed.data = input.data %>%
  correct.loq(
  by = "subject",
  nomtimevar = "ntad",
  timevar = "tad",
  depvar = "conc",
  bloqvar = "bloq",
  loqvar = "loq",
  loqrule = 1
)

# Result:

loqed.data %>% filter(loqrule.nr!="") %>% select(subject,ntad,conc,loqrule.nr,loqrule.txt) %>% kable()


## ---- results="markup", warnings=F--------------------------------------------

th = loqed.data %>% 
     est.thalf(
     by = "subject",
     timevar = "tad",
     depvar = "conc",
     includeCmax = "Y",
     exclvar = "excl_th"
     )

# Result:

head(th) %>% kable()


## ---- results="markup", warnings=F, fig.width = 7-----------------------------

plot_reg(
  loqed.data,
  by = "subject",
  th = th,
  bloqvar = "bloq",
  timevar = "tad",
  depvar = "conc",
  timelab = "Time (h)",
  deplab = "Conc (ng/mL)",
  exclvar = "excl_th",
  plotdir = NA
)


## ---- results="markup", warnings=F--------------------------------------------

ctmax = input.data %>% calc.ctmax(
  by = "subject",
  timevar="tad",
  depvar="conc"
)

# Result:

head(ctmax) %>% kable()


## ----  results="markup", warnings=F-------------------------------------------

ct.data= loqed.data %>%
    correct.time(
    by="subject",
    nomtimevar="ntad",
    timevar="tad",
    depvar="conc",
    th=th,
    tau=24,
    tstart=4,
    tend=9,
    teval=12,
    reg="sd",
    method=1
  ) 

# Result:

ct.data %>% filter(subject<=5&(trule.nr!=""|create.nr!="")) %>% select(subject,ntad,conc,applies.to.time,trule.nr,trule.txt,create.txt) %>% kable()


## ----  results="markup", warnings=F-------------------------------------------

cc.ct.data = ct.data %>%
    correct.conc(
    by ="subject",
    nomtimevar="ntad",
    tau=24,
    tstart=4,
    tend=9,
    teval=12,
    th=th,
    reg="sd",
    ss="n",
    route="EV",
    method=1
  )

# Result:

cc.ct.data %>% filter(crule.nr!="") %>% select(subject,ntad,conc,applies.to.conc,crule.nr,crule.txt) %>% kable()


## ----  results="markup", warnings=F-------------------------------------------

tab_corr = cc.ct.data %>% 
     tab.corr(
       by="subject",
       nomtimevar="ntad"
     )

# Result:

tab_corr %>% kable()


## ----  results="markup", warnings=F-------------------------------------------

par = cc.ct.data %>%
      calc.par(by = 'subject',
               tau=24,
               teval=12,
               tstart=4,
               tend=9,
               route="EV",
               method=1)

# Result:

head(par) %>% kable()


## ----  results="markup", warnings=F-------------------------------------------

# Create a covariates file, containing at least the dose given

cov = input.data %>%
      distinct(subject,dose)

par <- par %>% 
  calc.par.th(
    by="subject",
    th=th ,
    covariates=cov,
    dose="dose",
    factor=1,  
    reg="sd",
    ss="n"
    )

# Result:

head(par) %>% kable()


## ----  results="markup", warnings=F-------------------------------------------

par_all = left_join(par, ctmax)

# Result:

head(par_all) %>% kable()


