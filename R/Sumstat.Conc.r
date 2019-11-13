#' Creates summary statistics of concenctrations by time variable and grouping variable
#'
#' @param ds
#' @param dvVar
#' @param timeVar
#' @param by
#' @param nsig
#' @param na.rm
#'
#' @examples
#'library(dplyr)
#'Theoph1 <- Theoph %>% mutate(FORM=ifelse(Dose<5, "LowDose", "HighDose"))
#'
#'NTAD <- c(0,0.3,0.5,1,2,4,5,7,9,12,24)
#'Theoph1 <- Theoph1 %>%
#'  mutate(NTAD=metrumrg::snap(Time, NTAD)) %>%
#'  mutate(Subject=as.numeric(as.character(Subject)),     #converting from factor to numeric
#'         BQL = ifelse(conc<=0.25, 1, 0),                #just adding few BLQs to demonstrate functionality
#'         conc= ifelse(conc<=0.25, NA, round(conc, 2)))  #just adding digits to demonstrate nsig functionality
#'
#'test <- Sumstat.Conc(Theoph1, timeVar="NTAD", dvVar="conc", nsig=3,  na.rm=TRUE, by=NULL)
#'
#' @export
Sumstat.Conc = function(ds, dvVar, timeVar, by, nsig, na.rm){

  by1= rlang::syms(by)
  timeVar=rlang::sym(timeVar)
  by=c(by1, timeVar)
  dvVar= rlang::sym(dvVar)

  out = ds %>% group_by(!!!by) %>%
    summarize(N= n(),
              Mean=signif(mean(!!dvVar, na.rm=na.rm),nsig),
              SD=signif(sd(!!dvVar, na.rm=na.rm), nsig),
              Median=signif(median(!!dvVar, na.rm=na.rm), nsig),
              Min=signif(min(!!dvVar, na.rm=na.rm), nsig),
              Max=signif(max(!!dvVar, na.rm=na.rm), nsig),
              L95CI=signif((qt(0.025, N) * (sqrt(var(!!dvVar, na.rm=na.rm))/sqrt(N))) + Mean, nsig),
              U95CI=signif((qt(0.975, N) * (sqrt(var(!!dvVar, na.rm=na.rm))/sqrt(N))) + Mean, nsig),
              GeoMean=signif(qpToolkit::geomean(!!dvVar, na.rm=na.rm), nsig))

  out <- as.data.frame(lapply(out, function(y) gsub("NaN", "NC", y)))
  out <- as.data.frame(lapply(out, function(y) gsub("-Inf", "NC", y)))
  out <- as.data.frame(lapply(out, function(y) gsub("Inf", "NC", y)))
  out <- as.data.frame(lapply(out, function(y) ifelse(is.na(y), "NC", y)))

  out <- out %>%
    tidyr::gather(c("N", "Mean", "SD", "Median", "Min", "Max", "L95CI", "U95CI", "GeoMean"),
                  key="key", value="value") %>%
    tidyr::spread(!!timeVar, "value")

  out$key <- factor(out$key,
                    c("N", "Mean", "SD", "Median", "Min", "Max", "L95CI", "U95CI", "GeoMean"),
                    levels=c("N", "Mean", "SD", "Median", "Min", "Max", "L95CI", "U95CI", "GeoMean"))

  out <- out %>% arrange(!!!by1, key) %>% select(!!!by1, Statistics=key, everything())

  return(out)
}
