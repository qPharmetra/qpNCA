library(readxl)
library(dplyr)
library(magrittr)
library(csv)
library(tidyr)
sheets=excel_sheets('input_rules_1.0.12.xlsx')
sheet_list_org=lapply(sheets, function(x) {read_excel('input_rules_1.0.12.xlsx',sheet=x)})
names(sheet_list_org)=sheets
list2env(sheet_list_org, envir=.GlobalEnv)
render <- function(..., path)bind_rows(...) %>%
  mutate(dv = as.numeric(dv)) %>%
  mutate(loq=as.numeric(loq), bloq=ifelse(is.na(loq),0,1),excl_th=0) %>%
  as.csv(path)

render(path = '../../tests/testthat/posd.csv',`START_POSD`,`SDT_1`,`SDT_2`,`SDT_3`,`SDC_1`,`SDC_2`,`SDC_3`)
render(path = '../../tests/testthat/ivsd.csv',`START_IVSD`,`SDC_4a`,`SDC_4b`,`SDC_4c`,`SDC_4d`)
render(path = '../../tests/testthat/pomdnoss.csv',`START_POMD`,`MDT_1`,`MDT_2`,`MDT_3`,`MDT_3a`,`MDC_1_noss`,`MDC_2_noss`,`MDC_3_noss`)
render(path = '../../tests/testthat/pomdss.csv',`MDC_1_ss`,`MDC_2_ss`,`MDC_3_ss`)
render(path = '../../tests/testthat/ivmd.csv',`MDC_4a`,`MDC_4b`,`MDC_4c`,`MDC_4d`)


sum_excel <-read_excel(
  'validation_qPNCA_1.0.12.xlsx',
  sheet='summ',
  guess_max = 100
)

sum_excel %<>%
  select(-2) %>% # remove empty column
  filter(!is.na(Parameter)) %>%
  gather('rule','value',2:length(.)) %>%
  spread(Parameter,value) %>%
  mutate(id=as.numeric(id)) %>%
  arrange(id) %>%
  select(id,rule,everything()) %>%
  gather('parameter','value_reference',-id,-rule,-dose,-factor,-includeCmax,-method,-reg,-route,-ss,-tau,-teval,-tstart,-tend) %>%
  arrange(id,rule,parameter) %>%
  mutate_at(vars(dose,factor,method,tau,teval,tstart,tend),as.numeric) %>%
  mutate_at(vars(reg,ss),tolower) %>%
  mutate(value_reference=ifelse(parameter=='intercept',exp(as.numeric(value_reference)),as.numeric(value_reference)))
# exponentiate intercept before comparing to qPNCA result


# rename some parameters to match v1.0.22 terminology
mutate_cond <-
function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
sum_excel %<>%
  mutate_cond(parameter=='cl.f.obs'&route=='IV'&ss=='n',parameter='cl.obs') %>%
  mutate_cond(parameter=='cl.f.obs'&route=='IV'&ss=='y',parameter='cl.ss') %>%
  mutate_cond(parameter=='cl.f.pred'&route=='IV'&ss=='n',parameter='cl.pred') %>%
  mutate_cond(parameter=='cl.f.pred'&route=='PO'&ss=='y',parameter='cl.f.ss') %>%
  mutate_cond(parameter=='cl.f.obs'&route=='PO'&ss=='y',value_reference=NA) %>%
  mutate_cond(parameter=='vz.f.obs'&route=='IV'&ss=='n',parameter='vz.obs') %>%
  mutate_cond(parameter=='vz.f.obs'&route=='IV'&ss=='y',parameter='vss.obs') %>%
  mutate_cond(parameter=='vz.f.pred'&route=='IV'&ss=='n',parameter='vz.pred') %>%
  mutate_cond(parameter=='vz.f.pred'&route=='IV'&ss=='y',parameter='vss.pred') %>%
  mutate_cond(parameter=='vss.obs'&route=='PO',value_reference=NA) %>%
  mutate_cond(parameter=='vss.pred'&route=='PO',value_reference=NA)

sum_excel %>% as.csv('../../tests/testthat/voucher.csv')
