library(readxl)
library(dplyr)
library(magrittr)
library(csv)
library(tidyr)
library(wrangle)
prev_sheets=excel_sheets('input_rules_1.0.12.xlsx')
sheets=excel_sheets('input_rules_1.0.22.xlsx')
stopifnot(setequal(sheets, prev_sheets))
sheet_list_org=lapply(sheets, function(x) {read_excel('input_rules_1.0.22.xlsx',sheet=x)})
names(sheet_list_org)=sheets

sheet_list_org[[6]]['dv']
sheet_list_org[[6]] %<>% mutate(dv = as.numeric(dv))
sheet_list_org[[3]]['dv']
sheet_list_org[[3]] %<>% mutate(dv = as.numeric(dv))
sheet_list_org %>% do.call(bind_rows,.) %>% write.csv('../../tests/testthat/profiles.csv')

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


used <- c(
  'START_POSD','SDT_1','SDT_2','SDT_3','SDC_1','SDC_2','SDC_3',
  'START_IVSD','SDC_4a','SDC_4b','SDC_4c','SDC_4d',
  'START_POMD','MDT_1','MDT_2','MDT_3','MDT_3a','MDC_1_noss','MDC_2_noss','MDC_3_noss',
  'MDC_1_ss','MDC_2_ss','MDC_3_ss','MDC_4a','MDC_4b','MDC_4c','MDC_4d'
)

setdiff(sheets, used)

START_INFSD$dv %<>% as.numeric
START_INFMD$dv %<>% as.numeric

render(
  path = '../../tests/testthat/other.csv',
  `START_INFSD`, # problematic
  `START_IVMD`,  # ok
  `START_INFMD`, # problematic
  `POSDLOQ_1`,
  `POSDLOQ_2`,   `POSDLOQ_3`,   `POSDLOQ_4`,   `IVSDLOQ_1`,
  `IVSDLOQ_2`,   `IVSDLOQ_3`,   `IVSDLOQ_4`,   `POSDEXCL`,
  `POSDTRAP_2`,  `POSDTRAP_3`
)

old_sum_excel <-read_excel(
  'validation_qPNCA_1.0.12.xlsx',
  sheet='summ',
  guess_max = 100
)
sum_excel <-read_excel(
  'validation_qPNCA_1.0.22.xlsx',
  sheet='summ',
  guess_max = 100
)

dim(old_sum_excel)
dim(sum_excel)

setdiff(names(old_sum_excel), names(sum_excel))
setdiff(names(sum_excel), names(old_sum_excel))

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
