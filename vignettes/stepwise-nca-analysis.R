## ---- results='hide'----------------------------------------------------------
library(dplyr)
library(qpNCA)

## ---- results="markup", warnings=F--------------------------------------------
head(Theoph)

#we need nominal time variable for some tasks.
NTAD <- c(0,0.25,0.5,1,2,4,5,7,9,12,24)

# or use metrumrg::snap
for(i in 1:nrow(Theoph)){
  time  <- Theoph$Time[[i]]
  delta <- abs(NTAD - time)
  best  <- min(delta)
  index <- match(best, delta)
  nom   <- NTAD[[index]]
  Theoph$NTAD[[i]] <- nom
}

Theoph1 <- Theoph %>% 
  # mutate(NTAD=metrumrg::snap(Time, NTAD)) %>%
  mutate(Subject=as.numeric(Subject)) %>%
  mutate(BLQ = ifelse(conc<0.20, 1, 0),
         LOQ = 0.20) %>%
  mutate(conc = ifelse(BLQ==1, 0, conc)) %>%
  as.data.frame() %>%
  arrange(Subject, Time)
              

## ---- results="markup", warnings=F--------------------------------------------
ctmax <- Theoph1 %>% calc.ctmax(
  by = 'Subject',
  timevar="Time",
  depvar="conc"
)
head(ctmax)

## ----  results="markup", warnings=F-------------------------------------------
th = Theoph1 %>% est.thalf(
  by = 'Subject',
  timevar="Time",
  depvar="conc",
  includeCmax="Y"
)
head(th)

## ----  results="markup", warnings=F-------------------------------------------

## 3. Correct deviations
tc = Theoph1 %>%
  correct.loq(
    by='Subject',
    nomtimevar="NTAD",
    timevar="Time",
    depvar="conc",
    bloqvar="BLQ",
    loqvar="LOQ",
    loqrule=1
  ) %>%
  correct.time(
    by='Subject',
    nomtimevar="NTAD",
    timevar="Time",
    depvar="conc",
    th=th,
    reg="sd"
  ) %>%  
  correct.conc(
    by ='Subject',
    nomtimevar="NTAD",
    teval=12,
    th=th,
    reg="sd",
    ss="n",
  )

head(tc)


## ----  results="markup", warnings=F-------------------------------------------
par <- tc %>% calc.par(by = 'Subject', teval=12, route="EV")  
  #This wil get us both AUC0-12 and AUC0-24 as sampling ends at 24
head(par)

## ----  results="markup", warnings=F-------------------------------------------
cov <- data.frame(Subject=as.numeric(Theoph1$Subject), DOSE=Theoph$Dose) %>% 
  distinct(.,.keep_all = T)

par <- par %>% 
  calc.par.th(
    th=th ,
    cov=cov,
    dose="DOSE",
    factor=1,  
    reg="sd",
    ss="n",
    by='Subject'
  )
head(par)

## ----  results="markup", warnings=F-------------------------------------------
par_all = left_join(par, ctmax)
head(par_all)
dir <- tempdir()
write.csv(par_all,file.path(dir,"nca_analysis_stepbystep_results.csv"), row.names = F)

