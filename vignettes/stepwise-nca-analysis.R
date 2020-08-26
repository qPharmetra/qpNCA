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
ctmax <- Theoph1 %>%
 group_by(Subject) %>% 
  #Subject in Theoph are ordered factor, but factors are not in ascending order
 do(calc.ctmax(.,timevar="Time",depvar="conc")) %>%
 ungroup()
head(ctmax)

## ----  results="markup", warnings=F-------------------------------------------
th = Theoph1 %>% 
  group_by(Subject) %>%
  do(est.thalf(.,timevar="Time",depvar="conc",includeCmax="Y")) %>%
  ungroup()
head(th)

## ----  results="markup", warnings=F-------------------------------------------

## 3. Correct deviations
tc = Theoph1 %>%
  group_by(Subject) %>%
  do(correct.loq(.,nomtimevar="NTAD",timevar="Time",depvar="conc",
                   bloqvar="BLQ",loqvar="LOQ",loqrule=1)) %>%
  do(correct.time(.,nomtimevar="NTAD",timevar="Time",depvar="conc",
                 tau=,tstart=,tend=,teval=,th=th,reg="sd",by='Subject')) %>%  
  do(correct.conc(.,nomtimevar="NTAD",tau=,tstart=,tend=,teval=12,
                   th=th,reg="sd",ss="n",by='Subject')) %>%
  ungroup()

head(tc)


## ----  results="markup", warnings=F-------------------------------------------
par <- tc %>%
  group_by(Subject) %>%
  do(calc.par(.,tau=,tstart=,tend=,teval=12, route="EV")) %>%  
  #This wil get us both AUC0-12 and AUC0-24 as sampling ends at 24
  ungroup()
head(par)

## ----  results="markup", warnings=F-------------------------------------------
cov <- data.frame(Subject=as.numeric(Theoph1$Subject), DOSE=Theoph$Dose) %>% 
  distinct(.,.keep_all = T)

par = calc.par.th(x=par,th=th ,cov=cov,
                  dose="DOSE",factor=1,  
                  reg="sd",ss="n",by='Subject') 
head(par)

## ----  results="markup", warnings=F-------------------------------------------
par_all = left_join(par, ctmax)
head(par_all)
# write.csv(par_all, "nca_analysis_stepbystep_results.csv", row.names = F)
