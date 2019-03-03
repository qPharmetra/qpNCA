#' Generate summary statistics of concentrations
#'
#' @param x
#' @param LOQ
#' @param na.rm
#' @param ns
#' @param bloqRule
#' @param oddCode
#' @param bloqCode
#' @param ...
#'
#' @return
#' @export
nca.sumstat.conc = function(x, LOQ, na.rm = T, ns = 3, bloqRule = "BLOQ=0", oddCode = Cs(M,NS), bloqCode="BLOQ", ...){
  # Summarize a vector concentration measurements for output to a report
  #   x     = character vector containing numeric values to be summarized
  #   LOQ   = limit of quantitation for x
  #   na.rm = flag for removal of NA values
  #   ns    = number of significant figures for summary values
  #   bloqRule  = Rule for treatment of BLOQ values
  #   oddCode   = alphabetical codes for concentrations
  #   bloqCode  = alphabetic code for BLOQ values

  # prepopulate the output assuming nothing needs to be summarized (i.e. all(x)=="NA")
  N = 0
  Mean = geoMean = SD = Min = Median = Max = CV = UCI = LCI = "NC"
  cvec = c(N, Mean, geoMean, SD, Min, Median, Max, CV, UCI, LCI)

  ## First apply BLOQ rule to the input data
  ## Keep in mind that the function expects character data
  ## One option is to set to 0
  if(bloqRule == "BLOQ=0") x[x == bloqCode] = "0"

  # Option2 is BLOQ <- LOQ/2
  if(bloqRule == "BLOQ=LOQ/2") x[x == bloqCode] = as.character(LOQ/2)

  ## make character entries NA (not "NA") before turning character into numeric
  if(length(whichNumeric(x))>0)
    x[-whichNumeric(x)] = NA

  ## check if we have evaluable numbers
  nn = length(x[whichNumeric(x)])
  if(all(x[!is.na(x)]=="0" & bloqRule=="BLOQ=0")) nn = -1
  if(all(2*as.numeric(x[!is.na(x)])==LOQ & bloqRule=="BLOQ=LOQ/2")) nn = -1

  # Set proper number of significant figures
  ns2 = ns+1
  myfmt1 = paste("%#.", ns, "g", sep="")
  myfmt2 = paste("%#.", ns2, "g", sep="")

  if(nn == -1){
    # all BLOQ
    cvec = c(as.character(length(x[whichNumeric(x)])),
             c(bloqCode, bloqCode, "NC", rep(bloqCode, 3), "NC", "NC", "NC"))
  } else if(nn == 0) { # No observations or all NA
    cvec = rep("NC", 10)
  } else if(mean(as.numeric(x), na.rm=T) < LOQ){
    # mean BLOQ
    x = as.numeric(x[whichNumeric(x)])
    cvec = c(N = length(x[whichNumeric(x)]),
             Mean = bloqCode,
             geoMean = bloqCode,
             SD = "NC",
             Min = bloqCode,
             Median = bloqCode,
             Max = sprintf(myfmt1, signif(max(as.numeric(x[whichNumeric(x)]), na.rm=na.rm),ns)),
             CV = "NC",
             UCI = "NC",
             LCI = "NC")
  } else if(nn == 0) {
    # No observations
    cvec = rep("NC", 10)
  } else if (nn <= 2){
    # Few observations
    N = sprintf("%.0f", nn)
    x = as.numeric(x[whichNumeric(x)])
    Mean = sprintf(myfmt2, signif(mean(x, na.rm = na.rm), ns2))
      if(grepl("e", Mean)) Mean = sprintf("%g", signif(mean(x, na.rm = na.rm), ns2))
    geoMean = sprintf(myfmt2, signif(geomean(x[x>0], na.rm = na.rm), ns2))
      if(grepl("e", geoMean)) geoMean = sprintf("%g", signif(geomean(x, na.rm = na.rm), ns2))
    SD = "NC"
    Min = sprintf(myfmt1, signif(min(x, na.rm=na.rm),ns))
    if(Min < LOQ) Min = bloqCode
      if(grepl("e", Min)) Min = sprintf("%g", signif(min(x, na.rm=T), ns)) # correct for scientific notation
    Median = sprintf(myfmt2, signif(median(x, na.rm=na.rm), ns2))
      if(Median < LOQ) Median = bloqCode
      if(grepl("e", Median)) Median = sprintf("%g", signif(mean(x, na.rm = na.rm), ns2))
    Max = sprintf(myfmt1, signif(max(x, na.rm=na.rm), ns))
      if(grepl("e", Max)) Max = sprintf("%g", signif(max(x, na.rm=T), ns)) # correct for scientific notation
    CV = "NC"
    UCI = "NC"
    LCI = "NC"
    ## Format the summary stats
    cvec = c(N, Mean, geoMean, SD, Min, Median, Max, CV, UCI, LCI)
    invalids = c(grep("NA", cvec),
                 grep("NaN", cvec),
                 grep("Inf", cvec))

    cvec[invalids] = "NC"
    ## get rid of trailing decimal dots
    cvec = sub("(.*)\\.$","\\1",cvec)
  } else {
    # normal case
    N = sprintf("%.0f", nn)
    x = as.numeric(x[whichNumeric(x)])
    Mean = sprintf(myfmt2, signif(mean(x, na.rm = na.rm), ns2))
      if(grepl("e", Mean)) Mean = sprintf("%g", signif(mean(x, na.rm = na.rm), ns2))
    if(all(x[!is.na(x)]>0)) {
      geoMean = sprintf(myfmt2, signif(geomean(x, na.rm = na.rm),ns2))
      if(grepl("e", geoMean)) geoMean = sprintf("%g", signif(geomean(x, na.rm = na.rm), ns2))
    } else geoMean = "NC"
    SD = sprintf(myfmt2, signif(sd(x, na.rm=na.rm),ns2))
      if(grepl("e", SD)) SD = sprintf("%g", signif(sd(x, na.rm=na.rm), ns)) # correct for scientific notation
    Min = sprintf(myfmt1, signif(min(x, na.rm=na.rm), ns))
      if(as.numeric(Min) < LOQ) Min = bloqCode
      if(grepl("e", Min)) Min = sprintf("%g", signif(min(x, na.rm=na.rm), ns)) # correct for scientific notation
    Median = sprintf(myfmt2, signif(median(x, na.rm=na.rm),ns2))
      if(as.numeric(Median) < LOQ) Median = bloqCode
      if(grepl("e", Median)) Median = sprintf("%g", signif(median(x, na.rm=na.rm), ns2)) # correct for scientific notation
    Max = sprintf(myfmt1, signif(max(x, na.rm=na.rm), ns))
      if(grepl("e", Max)) Max = sprintf("%g", signif(max(x, na.rm=na.rm), ns)) # correct for scientific notation
    CV = sprintf("%.1f", sd(x, na.rm=na.rm)/mean(x, na.rm = na.rm)*100)
    se = sd(x, na.rm=na.rm)/sqrt(nn)
    uci = mean(x, na.rm = na.rm) + qt(p=.975, df=nn-1)*se
    lci = mean(x, na.rm = na.rm) - qt(p=.975, df=nn-1)*se
    UCI = sprintf(myfmt2, signif(uci, ns2))
      if(grepl("e", UCI)) UCI = sprintf("%g", signif(uci, ns2)) # correct for scientific notation
    LCI = sprintf(myfmt2, signif(lci, ns2))
      if(grepl("e", LCI)) LCI = sprintf("%g", signif(lci, ns2)) # correct for scientific notation

    # correct trailing decimals
    if(substring(LCI, nchar(LCI)) == ".") LCI = paste(LCI, "0", sep="")
    if(substring(UCI, nchar(UCI)) == ".") UCI = paste(UCI, "0", sep="")

    ## Format the summary stats
    cvec = c(N, Mean, geoMean, SD, Min, Median, Max, CV, LCI, UCI)
    invalids = c(grep("NA", cvec),
                 grep("NaN", cvec),
                 grep("Inf", cvec))

    cvec[invalids] = "NC"
    ## get rid of trailing decimal dots
    cvec = sub("(.*)\\.$","\\1",cvec)
  }

  names(cvec) = c(Cs(N, Mean, "Geometric Mean", SD, Min, Median, Max), "CV(\\%)", "Lower 95\\% CI", "Upper 95\\% CI")
  return(cvec)
}
