---
title: "stepwise-nca-analysis.rmd"
author: "Krina Mehta"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Stepwise NCA Analysis}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::knitr}
editor_options: 
  chunk_output_type: console
---

#Load libraries

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(qpNCA)
```

#Check out Theoph data and prepare data for NCA

```r
head(Theoph)
```

```
##   Subject   Wt Dose Time  conc
## 1       1 79.6 4.02 0.00  0.74
## 2       1 79.6 4.02 0.25  2.84
## 3       1 79.6 4.02 0.57  6.57
## 4       1 79.6 4.02 1.12 10.50
## 5       1 79.6 4.02 2.02  9.66
## 6       1 79.6 4.02 3.82  8.58
```

```r
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
```

#Calculate Cmax and Tmax

```r
ctmax <- Theoph1 %>% calc.ctmax(
  by = 'Subject',
  timevar="Time",
  depvar="conc"
)
head(ctmax)
```

```
## # A tibble: 6 x 3
##   Subject  cmax  tmax
##     <dbl> <dbl> <dbl>
## 1       1  6.44  1.15
## 2       2  7.09  3.48
## 3       3  7.56  2.02
## 4       4  8     0.98
## 5       5  8.2   1.02
## 6       6  8.33  1.92
```

#Calculate half-life

```r
th = Theoph1 %>% est.thalf(
  by = 'Subject',
  timevar="Time",
  depvar="conc",
  includeCmax="Y"
)
head(th)
```

```
## # A tibble: 6 x 11
##   Subject no.points intercept lambda_z r.squared adj.r.squared start_th end_th
##     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>         <dbl>    <dbl>  <dbl>
## 1       1         3      8.21   0.0916     0.999         0.998     9.22   23.8
## 2       2         4      9.86   0.0883     0.999         0.998     6.98   24.2
## 3       3         7      8.81   0.0818     0.992         0.991     2.02   24.1
## 4       4         3      8.56   0.0955     1.00          1.00      9.03   24.1
## 5       5         3     12.5    0.102      0.999         0.999     9      24.2
## 6       6         4     11.1    0.104      0.997         0.996     7.03   24.3
## # ... with 3 more variables: thalf <dbl>, includeCmax <chr>,
## #   points_excluded <chr>
```

#Corrections of missing time and concentrations

```r
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
```

```
## # A tibble: 6 x 70
##   Subject    Wt  Dose   BLQ   LOQ  conc  Time  NTAD bloqvar1 loqvar1 loqrule.nr
##     <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>   <dbl> <chr>     
## 1       1    80     4     1   0.2  0    0      0           1     0.2 "LOQ1"    
## 2       1    80     4     0   0.2  1.29 0.27   0.25        0     0.2 ""        
## 3       1    80     4     0   0.2  3.08 0.580  0.5         0     0.2 ""        
## 4       1    80     4     0   0.2  6.44 1.15   1           0     0.2 ""        
## 5       1    80     4     0   0.2  6.32 2.03   2           0     0.2 ""        
## 6       1    80     4     0   0.2  5.53 3.57   4           0     0.2 ""        
## # ... with 59 more variables: loqrule.txt <chr>, firstmeast <dbl>,
## #   consecutive <dbl>, lambda_z.x <dbl>, start_th.x <dbl>, end_th.x <dbl>,
## #   includeCmax.x <chr>, points_excluded.x <chr>, create.nr <chr>,
## #   create.txt <chr>, trule.nr <chr>, trule.txt <chr>, applies.to.time <chr>,
## #   t0.flag <dbl>, tau.flag <dbl>, tstart.flag <dbl>, tend.flag <dbl>,
## #   teval.flag <dbl>, missflag <dbl>, misstime <lgl>, diff <dbl>, flag <dbl>,
## #   leaddv <dbl>, lagdv <dbl>, leadtime <dbl>, lagtime <dbl>, conc.tau <dbl>,
## #   time.tau <dbl>, conc.teval <dbl>, time.teval <dbl>, conc.part <dbl>,
## #   time.part <dbl>, conc.lastall <dbl>, time.lastall <dbl>, newdepvar <lgl>,
## #   lambda_z.y <dbl>, start_th.y <dbl>, end_th.y <dbl>, includeCmax.y <chr>,
## #   points_excluded.y <chr>, ptime <dbl>, crule.nr <chr>, crule.txt <chr>,
## #   applies.to.conc <chr>, lambda_z <lgl>, depvar <dbl>, timevar <dbl>,
## #   lead.cteval <dbl>, lag.cteval <dbl>, lead.tteval <dbl>, lag.tteval <dbl>,
## #   t0val <dbl>, tauval <lgl>, back_extrap <dbl>, lc1 <dbl>, lc2 <dbl>,
## #   lt1 <dbl>, lt2 <dbl>, firstmeasc <dbl>
```

#calculate PK parameters NOT based on lambda_z 

```r
par <- tc %>% calc.par(by = 'Subject', teval=12, route="EV")  
  #This wil get us both AUC0-12 and AUC0-24 as sampling ends at 24
head(par)
```

```
## # A tibble: 6 x 26
##   Subject route method tlast clast.obs tlast.ok t0.ok aucall auclast aumcall
##     <dbl> <chr>  <dbl> <dbl>     <dbl>    <dbl> <dbl>  <dbl>   <dbl>   <dbl>
## 1       1 EV         1  23.8      0.92        1     1   73.8    73.8    609.
## 2       2 EV         1  24.2      1.15        1     1   90.7    90.7    782.
## 3       3 EV         1  24.1      1.25        1     1   88.6    88.6    740.
## 4       4 EV         1  24.1      0.86        1     1   80.1    80.1    617.
## 5       5 EV         1  24.2      1.05        1     1   99.3    99.3    803.
## 6       6 EV         1  24.3      0.9         1     1   91.5    91.5    707.
## # ... with 16 more variables: aumclast <dbl>, mrtall <dbl>, mrtlast <dbl>,
## #   calc.tau <dbl>, auctau <lgl>, aumctau <lgl>, tau <lgl>, calc.teval <dbl>,
## #   auc12 <dbl>, teval <dbl>, calc.part <dbl>, aucpart <lgl>, tstart <lgl>,
## #   tend <lgl>, c0 <lgl>, area.back.extr <lgl>
```

#calculate PK parameters based on lambda_z 

```r
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
```

```
## # A tibble: 6 x 63
##   Subject method tlast clast.obs tlast.ok t0.ok aucall auclast aumcall aumclast
##     <dbl>  <dbl> <dbl>     <dbl>    <dbl> <dbl>  <dbl>   <dbl>   <dbl>    <dbl>
## 1       1      1  23.8      0.92        1     1   73.8    73.8    609.     609.
## 2       2      1  24.2      1.15        1     1   90.7    90.7    782.     782.
## 3       3      1  24.1      1.25        1     1   88.6    88.6    740.     740.
## 4       4      1  24.1      0.86        1     1   80.1    80.1    617.     617.
## 5       5      1  24.2      1.05        1     1   99.3    99.3    803.     803.
## 6       6      1  24.3      0.9         1     1   91.5    91.5    707.     707.
## # ... with 53 more variables: mrtall <dbl>, mrtlast <dbl>, calc.tau <dbl>,
## #   auctau <lgl>, aumctau <lgl>, tau <lgl>, calc.teval <dbl>, auc12 <dbl>,
## #   teval <dbl>, calc.part <dbl>, aucpart <lgl>, tstart <lgl>, tend <lgl>,
## #   c0 <lgl>, area.back.extr <lgl>, no.points <dbl>, intercept <dbl>,
## #   lambda_z <dbl>, r.squared <dbl>, adj.r.squared <dbl>, start_th <dbl>,
## #   end_th <dbl>, thalf <dbl>, includeCmax <chr>, points_excluded <chr>,
## #   DOSE <dbl>, factor <dbl>, reg <chr>, ss <chr>, route <chr>,
## #   clast.pred <dbl>, aucinf.obs <dbl>, aucinf.pred <dbl>, aumcinf.obs <dbl>,
## #   aumcinf.pred <dbl>, cl.f.obs <dbl>, cl.f.pred <dbl>, mrtinf.obs <dbl>,
## #   mrtinf.pred <dbl>, vz.f.obs <dbl>, vz.f.pred <dbl>, vss.obs <lgl>,
## #   vss.pred <lgl>, pctextr.obs <dbl>, pctextr.pred <dbl>, pctback.obs <dbl>,
## #   pctback.pred <dbl>, cl.obs <lgl>, cl.pred <lgl>, vz.obs <lgl>,
## #   vz.pred <lgl>, cl.f.ss <lgl>, cl.ss <lgl>
```

#All NCA para in one df

```r
par_all = left_join(par, ctmax)
```

```
## Joining, by = "Subject"
```

```r
head(par_all)
```

```
## # A tibble: 6 x 65
##   Subject method tlast clast.obs tlast.ok t0.ok aucall auclast aumcall aumclast
##     <dbl>  <dbl> <dbl>     <dbl>    <dbl> <dbl>  <dbl>   <dbl>   <dbl>    <dbl>
## 1       1      1  23.8      0.92        1     1   73.8    73.8    609.     609.
## 2       2      1  24.2      1.15        1     1   90.7    90.7    782.     782.
## 3       3      1  24.1      1.25        1     1   88.6    88.6    740.     740.
## 4       4      1  24.1      0.86        1     1   80.1    80.1    617.     617.
## 5       5      1  24.2      1.05        1     1   99.3    99.3    803.     803.
## 6       6      1  24.3      0.9         1     1   91.5    91.5    707.     707.
## # ... with 55 more variables: mrtall <dbl>, mrtlast <dbl>, calc.tau <dbl>,
## #   auctau <lgl>, aumctau <lgl>, tau <lgl>, calc.teval <dbl>, auc12 <dbl>,
## #   teval <dbl>, calc.part <dbl>, aucpart <lgl>, tstart <lgl>, tend <lgl>,
## #   c0 <lgl>, area.back.extr <lgl>, no.points <dbl>, intercept <dbl>,
## #   lambda_z <dbl>, r.squared <dbl>, adj.r.squared <dbl>, start_th <dbl>,
## #   end_th <dbl>, thalf <dbl>, includeCmax <chr>, points_excluded <chr>,
## #   DOSE <dbl>, factor <dbl>, reg <chr>, ss <chr>, route <chr>,
## #   clast.pred <dbl>, aucinf.obs <dbl>, aucinf.pred <dbl>, aumcinf.obs <dbl>,
## #   aumcinf.pred <dbl>, cl.f.obs <dbl>, cl.f.pred <dbl>, mrtinf.obs <dbl>,
## #   mrtinf.pred <dbl>, vz.f.obs <dbl>, vz.f.pred <dbl>, vss.obs <lgl>,
## #   vss.pred <lgl>, pctextr.obs <dbl>, pctextr.pred <dbl>, pctback.obs <dbl>,
## #   pctback.pred <dbl>, cl.obs <lgl>, cl.pred <lgl>, vz.obs <lgl>,
## #   vz.pred <lgl>, cl.f.ss <lgl>, cl.ss <lgl>, cmax <dbl>, tmax <dbl>
```

```r
dir <- tempdir()
write.csv(par_all,file.path(dir,"nca_analysis_stepbystep_results.csv"), row.names = F)
```

