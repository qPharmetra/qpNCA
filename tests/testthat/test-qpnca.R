library(dplyr)
library(magrittr)
library(ggplot2)
library(qpToolkit)
library(tidyr)
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
test_results <- function(x, reg, ss, route){
  y <- x %>% distinct(id) %>% mutate(dose=1) %>% select(id,dose)
  z <- qpNCA(
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
    pdfdir =,
    timelab = 'Time (h)',
    deplab = 'Concentration (ng/mL)',
    tau = 24,
    tstart = 4,
    tend = 12,
    teval = 18,
    covfile = y,
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
    -calc.teval,-calc.tau,-t0.ok,-tlast.ok
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
  test <- test_results(as_csv('posd.csv'),reg="SD",ss="N",route="EV")
  refr <- reference_results(as_csv('voucher.csv'), id = -1, from = 1, to = 6, route = 'ev')
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
