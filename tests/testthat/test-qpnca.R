library(dplyr)
library(magrittr)
library(testthat)
library(tidyr)
library(qpNCA)
#library(knitr)
as_csv <- function (
  x, as.is = TRUE,
  na.strings = c("", "\\s",".", "NA"),
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
my_qpNCA <- function(
  x,
  by = 'id',
  nomtimevar = 'ntad',
  timevar = 'tad',
  depvar = 'dv',
  bloqvar = 'bloq',
  loqvar = 'loq',
  loqrule = 2,
  includeCmax = 'Y',
  exclvar = 'excl_th',
  plotdir = NULL,
  pdfdir = NA,
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
  method = 1
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
  pdfdir = pdfdir,
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
  method = 1
)
test_results <- function(x, reg, ss, route){
  y <- x %>% distinct(id) %>% mutate(dose=1) %>% select(id,dose)
  z <- my_qpNCA(
    x,
    by = 'id',
    nomtimevar = 'ntad',
    timevar = 'tad',
    depvar = 'dv',
    bloqvar = 'bloq',
    loqvar = 'loq',
    loqrule = 2,
    includeCmax = 'Y',
    exclvar = 'excl_th',
    plotdir = NULL,
    pdfdir = NA,
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
    method = 1
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

reference_results <- function(x, id, from, to, route){
  y <- x %>% filter(id == id|between(id, from, to))
  y %<>% mutate(route = route)
  y %<>% select(id,parameter,value_reference)
  suppressWarnings(y$value_reference %<>% as.numeric %>% round(4))
  y
}

merged_results <- function(x, y){
  x %<>% left_join(y)
  x %<>% mutate(identical = as.integer(value_test == value_reference))
  x %<>% mutate(identical = ifelse(is.na(value_reference) & is.na(value_test),1,identical))
  x %<>% mutate(identical = ifelse(is.na(identical),0,identical))
  x
}

test_that('PO SD results are stable',{
  test <- test_results(
    as_csv('posd.csv') %>% filter(id == 6),
    reg="SD",ss="N",route="EV")
 refr <- reference_results(as_csv('voucher.csv'), id = -1, from = 1, to = 6, route = 'ev')
 # refr <- reference_results(as_csv('voucher.csv'), id = 5, from = 5, to = 5, route = 'ev')
  comp <- merged_results(test,refr)
  diff <- comp %>% filter(identical == 0) %>% data.frame
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('IV SD results are stable',{
  test <- test_results(as_csv('ivsd.csv'),reg="SD",ss="N",route="IVB" )
  refr <- reference_results(as_csv('voucher.csv'), id = -2, from = 7, to = 10, route = 'ivb')
  comp <- merged_results(test,refr)
  diff <- comp %>% filter(identical == 0) %>% data.frame
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('PO MD non-steady state results are stable',{
  test <- test_results(as_csv('pomdnoss.csv'),reg="MD",ss="N",route="EV")
  refr <- reference_results(as_csv('voucher.csv'), id = -3, from = 11, to = 17, route = 'ev')
  comp <- merged_results(test,refr)
  diff <- comp %>% filter(identical == 0) %>% data.frame
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('PO MD steady state results are stable',{
  test <- test_results(as_csv('pomdss.csv'),reg="MD",ss="Y",route="EV")
  refr <- reference_results(as_csv('voucher.csv'), id = Inf, from = 18, to = 20, route = 'ev')
  comp <- merged_results(test,refr)
  diff <- comp %>% filter(identical == 0) %>% data.frame
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('IV MD results are stable',{
  test <- test_results(as_csv('ivmd.csv'),reg="MD",ss="N",route="IVB")
  refr <- reference_results(as_csv('voucher.csv'), id = -5, from = Inf, to = -Inf, route = 'ivb')
  comp <- merged_results(test,refr)
  diff <- comp %>% filter(identical == 0) %>% data.frame
  expect_identical(comp$value_test, comp$value_reference)
})

test_that('all results are stable',{
  skip_on_cran()
  test <- as_csv('profiles.csv')

  #test %>% itemize(rule, id, desc) %>% data.frame

  #undebug(qpNCA:::check.input)
#   test %>%
#     filter(id == 25) %>%
#     my_qpNCA(reg = 'sd', ss = 'n', route = 'IVB') %$%
#     pkpar %>% as.list
#
# head(test)

test %<>% mutate(
  reg = case_when(
    grepl('SD', rule) ~ 'SD',
    grepl('MD', rule) ~ 'MD'
  )
)

test %<>% mutate(ss = 'n')

test %<>% mutate(
  route = case_when(
    grepl('INF', rule) ~ 'IVI',
    grepl('IV', rule) ~ 'IVB',
    grepl('PO', rule) ~ 'EV',
    grepl('SD', rule) ~ 'EV',
    grepl('MD', rule) ~ 'EV'
  )
)

test %<>% mutate(
  loqrule = case_when(
    grepl('LOQ_1',rule) ~ 1,
    grepl('LOQ_2',rule) ~ 2,
    grepl('LOQ_3',rule) ~ 3,
    grepl('LOQ_4',rule) ~ 4,
    TRUE ~ 2
  )
)

test %<>% mutate(method = 1)

test %>%
  select(id, rule, desc, reg, route, loqrule, method) %>%
  unique %>%
  data.frame %>% write.csv(row.names = FALSE, 'scenarios.csv')

out <- qpNCA(
  test,
  by = 'id',
  nomtimevar = 'ntad',
  timevar = 'tad',
  depvar = 'dv',
  bloqvar = 'bloq',
  loqvar = 'loq',
  #loqrule = 2,
  includeCmax = 'Y',
  exclvar = 'excl_th',
  plotdir = NULL,
  pdfdir = NA,
  timelab = 'Time (h)',
  deplab = 'Concentration (ng/mL)',
  tau = 24,
  tstart = 4,
  tend = 12,
  teval = 18,
  covariates = test %>% distinct(id) %>% mutate(dose=1) %>% select(id,dose),
  factor = 1000,
  # reg = reg,
  # ss = ss,
  # route = route,
  # method = 1,
  dose = 'dose'
)

expect_equal_to_reference(file = '001.rds',out)

})

