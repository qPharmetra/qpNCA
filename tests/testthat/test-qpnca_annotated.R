library(dplyr)
library(magrittr)
library(testthat)
library(tidyr)
library(qpNCA)
#library(knitr)
as_csv <- function (
  x, as.is = TRUE,
  na.strings = c("", "\\s",".", "NA"),     #JH specific
  strip.white = TRUE,
  check.names = FALSE,
  ...
){
  args <- list(...)
  y <- do.call(
    utils::read.csv,
    c(
      list(
        x,
        as.is = as.is,
        na.strings = na.strings,
        strip.white = strip.white,
        check.names = check.names
      ),
      args
    )
  )
  y
}
my_qpNCA <- function(     # create function with constant arguments filled in
  x,
  by = 'id',
  nomtimevar = 'ntad',
  timevar = 'tad',
  depvar = 'dv',
  bloqvar = 'bloq',
  loqvar = 'loq',
  loqrule = loqrule,
  includeCmax = 'Y',
  exclvar = 'excl_th',
  plotdir = NULL,
  timelab = 'Time (h)',
  deplab = 'Concentration (ng/mL)',
  tau = 24,
  tstart = 4,
  tend = 12,
  teval = 18,
  covariates = x %>% distinct(id) %>% mutate(dose=1) %>% select(id,dose),
  dose = 'dose',
  factor = 1000,
  reg = reg,
  ss = ss,
  route = route,
  method = method
)qpNCA(
  x,
  by = by,
  nomtimevar = nomtimevar,
  timevar = timevar,
  depvar = depvar,
  bloqvar = bloqvar,
  loqvar = loqvar,
  loqrule = loqrule,
  includeCmax = includeCmax,
  exclvar = exclvar,
  plotdir = plotdir,
  timelab = timelab,
  deplab = deplab,
  tau = tau,
  tstart = tstart,
  tend = tend,
  teval = teval,
  covariates = covariates,
  dose = dose,
  factor = factor,
  reg = reg,
  ss = ss,
  route = route,
  method = method
)

test_results <- function(x, reg, ss, route, method, loqrule){      # function that actually runs my_qpNCA
  y <- x %>% distinct(id) %>% mutate(dose=1) %>% select(id,dose)
  z <- my_qpNCA(
    x,
    by = 'id',
    nomtimevar = 'ntad',
    timevar = 'tad',
    depvar = 'dv',
    bloqvar = 'bloq',
    loqvar = 'loq',
    loqrule = loqrule,
    includeCmax = 'Y',
    exclvar = 'excl_th',
    plotdir = NULL,
    timelab = 'Time (h)',
    deplab = 'Concentration (ng/mL)',
    tau = 24,
    tstart = 4,
    tend = 12,
    teval = 18,
    covariates = y,
    dose = 'dose',
    factor = 1000,
    reg = reg,
    ss = ss,
    route = route,
    method = method
  ) %$% pkpar
  z %<>%
    select(
    -area.back.extr,-r.squared,-calc.part,
    -calc.teval,-calc.tau,-t0.ok,-tlast.ok,
    -factor, -loqrule
  ) %>%
  gather(
    "parameter","value_test",-id,-dose,-includeCmax,
    -method,-reg,-route,-ss,-tau,-teval,-tstart,-tend
  ) %>%
  arrange(id,parameter) %>%
  mutate(id=as.numeric(id))
  suppressWarnings(z$value_test %<>% as.numeric %>% round(4))
  z
}

reference_results <- function(x, route){   # function that loads the reference results (Excel)
  y <- x %>% mutate(route = route)
  y %<>% select(id,parameter,value_reference)
  suppressWarnings(y$value_reference %<>% as.numeric %>% round(4))
  y
}

merged_results <- function(x, y){            # function that merges actual results and (Excel) reference results
  x %<>% left_join(y)
  x %<>% mutate(identical = as.integer(value_test == value_reference))
  x %<>% mutate(identical = ifelse(is.na(value_reference) & is.na(value_test),1,identical))
  x %<>% mutate(identical = ifelse(is.na(identical),0,identical))
  x
}

test_that('PO SD results are stable',{       # tests each "testing curves group"
  test <- test_results(
    'profiles.csv' %>%
      as_csv %>%
      filter(id%in%c(-1,1:6)),
    reg="SD",
    ss="N",
    route="EV",
    loqrule=2,
    method=1
  )
  refr <- reference_results(
    'reference.csv' %>%
      as_csv %>%
      filter(id%in%c(-1,1:6)),
    route = 'ev'
  )
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing",ncurve,"extravascular, single dose curves\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profiles\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('IV Bolus SD results are stable',{
  test <- test_results(
    as_csv('profiles.csv') %>%
      filter(id%in%c(-2,7:10)),
    reg="SD",
    ss="N",
    route="IVB",
    loqrule=2,
    method=1
  )
  refr <- reference_results(as_csv('reference.csv')%>% filter(id%in%c(-2,7:10)), route = 'ivb')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing",ncurve,"intravascular bolus, single dose curves\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profiles\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('IV Infusion SD results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(-3)),
                       reg="SD",ss="N",route="IVI",loqrule=2,method=1 )
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(-3)), route = 'ivi')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing",ncurve,"intravascular infusion, single dose curve\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('PO MD non-steady state results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(-4,11:17)),
                       reg="MD",ss="N",route="EV",loqrule=2,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(-4,11:17)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing",ncurve,"extravascular, multiple dose, non-steady state curves\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profiles\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('PO MD steady state results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(18:20)),
                       reg="MD",ss="Y",route="EV",loqrule=2,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(18:20)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing",ncurve,"extravascular, multiple dose, steady state curves\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profiles\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('IV Bolus MD results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(-5,32:35)),
                       reg="MD",ss="N",route="IVB",loqrule=2,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(-5,32:35)), route = 'ivb')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing",ncurve,"intravascular bolus, multiple dose curves\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profiles\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('IV Infusion MD results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(-6)) ,reg="MD",ss="N",route="IVI",loqrule=2,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(-6)), route = 'ivi')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing",ncurve,"intravascular infusion, multiple dose curve\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 1 results after EV adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(21)) ,reg="SD",ss="N",route="EV",loqrule=1,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(21)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 1 (",ncurve,"extravascular, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 2 results after EV adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(22)) ,reg="SD",ss="N",route="EV",loqrule=2,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(22)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 2 (",ncurve,"extravascular, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 3 results after EV adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(23)) ,reg="SD",ss="N",route="EV",loqrule=3,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(23)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 3 (",ncurve,"extravascular, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 4 results after EV adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(24)) ,reg="SD",ss="N",route="EV",loqrule=4,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(24)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 4 (",ncurve,"extravascular, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 1 results after IVB adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(25)) ,reg="SD",ss="N",route="IVB",loqrule=1,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(25)), route = 'ivb')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 1 (",ncurve,"intravascular bolus, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 2 results after IVB adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(26)) ,reg="SD",ss="N",route="IVB",loqrule=2,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(26)), route = 'ivb')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 2 (",ncurve,"intravascular bolus, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 3 results after IVB adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(27)) ,reg="SD",ss="N",route="IVB",loqrule=3,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(27)), route = 'ivb')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 3 (",ncurve,"intravascular bolus, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('LOQ Rule 4 results after IVB adm. are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(28)) ,reg="SD",ss="N",route="IVB",loqrule=4,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(28)), route = 'ivb')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing LOQ Rule 4 (",ncurve,"intravascular bolus, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('Exclude datapoints from lamba_z calculation results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(29)) ,reg="SD",ss="N",route="EV",loqrule=2,method=1)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(29)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing exclusion of data points from lambda_z calculation (",ncurve,"extravascular, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('Trapezoidal rule 2 results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(30)) ,reg="SD",ss="N",route="EV",loqrule=2,method=2)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(30)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing Trapezoidal Rule 2 (",ncurve,"extravascular, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('Trapezoidal rule 3 results are stable',{
  test <- test_results(as_csv('profiles.csv') %>% filter(id%in%c(31)) ,reg="SD",ss="N",route="EV",loqrule=2,method=3)
  refr <- reference_results(as_csv('reference.csv') %>% filter(id%in%c(31)), route = 'ev')
  comp <- merged_results(test,refr)
  ncurve=length(unique(comp$id))
  npar=length(unique(comp$parameter))
  ndiff=length(comp$identical[comp$identical==0])
  cat("\n")
  cat(paste("Testing Trapezoidal Rule 3 (",ncurve,"extravascular, single dose curve )\n"))
  cat(paste("Calculated",npar,"PK parameters for",ncurve,"concentration-time profile\n"))
  cat(paste("Compared these results with PK parameter values calculated using Microsoft Excel\n"))
  cat(paste("The comparison resulted in",ndiff,"differences\n"))
  expect_identical(comp$value_test, comp$value_reference)
})



