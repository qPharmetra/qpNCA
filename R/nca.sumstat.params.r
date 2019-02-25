nca.sumstat.params = function(x, ns = 3, na.rm=T, ...){
  # Summarize a vector PK Parameters for output to a report
  # For parameter values that have 0's reported, e.g. AUC when all concentratinos are 0, we 
  #   assign NC
  # function parameters
  #   x     = character vector containing numeric values to be summarized
  #   na.rm = flag for removal of NA values
  #   ns    = number of significant figures for summary values
  #   bloqRule  = Rule for treatment of BLOQ values
  #   oddCode   = alphabetical codes for concentrations
  #   bloqCode  = alphabetic code for BLOQ values
  
  # prepopulate the output assuming nothing needs to be summarized (i.e. all(x)=="NA")
  N = 0
  Mean = geoMean = SD = Min = Median = Max = CV = LCI = UCI = "NC"
  cvec = c(N, Mean, geoMean, SD, Min, Median, Max, CV, LCI, UCI)
  
  ## make character entries NA (not "NA") before turning character into numeric
  if(length(whichNumeric(x))>0)
    x[-whichNumeric(x)] = NA
  
  ## check if we have any evaluable numbers
  x[as.numeric(x)==0] = "NC"
  x = x[whichNumeric(x)]
  #x = x[as.numeric(x)>0]
  nn = length(x)
  if(nn==0) nn = -1
  
  # Set proper number of significant figures
  ns2 = ns+1
  myfmt1 = paste("%#.", ns, "g", sep="")
  myfmt2 = paste("%#.", ns2, "g", sep="")
  
  if(nn == -1){
    # all 0's
    cvec = c("0",rep("NC", 9))
  } else if (nn <= 2){ 
    # Few observations
    N = sprintf("%.0f", nn)
    x = as.numeric(x[as.numeric(x)>0])
    Mean = sprintf(myfmt2, signif(mean(x, na.rm = T),ns2))
      if(grepl("e", Mean)) Mean = sprintf("%g", signif(mean(x, na.rm=T), ns2)) # correct for scientific notation
    geoMean = sprintf(myfmt2, signif(geomean(x, na.rm = T),ns2))
      if(grepl("e", geoMean)) geoMean = sprintf("%g", signif(geomean(x, na.rm=T), ns2)) # correct for scientific notation
    SD = "NC"
    Min = sprintf(myfmt1, signif(min(x, na.rm=T),ns))
      if(grepl("e", Min)) Min = sprintf("%g", signif(min(x, na.rm=T), ns)) # correct for scientific notation
    Median = sprintf(myfmt2, signif(median(x, na.rm=T),ns2))
      if(grepl("e", Median)) Median = sprintf("%g", signif(median(x, na.rm=T), ns)) # correct for scientific notation
    Max = sprintf(myfmt1, signif(max(x, na.rm=T), ns))
      if(grepl("e", Max)) Max = sprintf("%g", signif(max(x, na.rm=T), ns)) # correct for scientific notation
    CV = "NC"
    LCI = "NC"    
    UCI = "NC"    
    ## Format the summary stats
    cvec = c(N, Mean, geoMean, SD, Min, Median, Max, CV, LCI, UCI)
    invalids = c(grep("NA", cvec),
                 grep("NaN", cvec),
                 grep("Inf", cvec))
    
    cvec[invalids] = "NC"
    ## get rid of trailing decimal dots
    cvec = sub("(.*)\\.$","\\1",cvec)
  } else {
    # Most common
    N = sprintf("%.0f", nn)
    x = as.numeric(x[as.numeric(x)>0][whichNumeric(x)])
    Mean = sprintf(myfmt2, signif(mean(x, na.rm = T), ns2))
      if(grepl("e", Mean)) Mean = sprintf("%g", signif(mean(x, na.rm=T), ns2)) # correct for scientific notation
    geoMean = sprintf(myfmt2, signif(geomean(x, na.rm = T),ns2))
      if(grepl("e", geoMean)) geoMean = sprintf("%g", signif(geomean(x, na.rm=T), ns2)) # correct for scientific notation
    SD = sprintf(myfmt2, signif(sd(x, na.rm=T), ns2))
      if(grepl("e", SD)) SD = sprintf("%g", signif(sd(x, na.rm=T), ns2)) # correct for scientific notation
    Min = sprintf(myfmt1, signif(min(x, na.rm=T), ns))
      if(grepl("e", Min)) Min = sprintf("%g", signif(min(x, na.rm=T), ns)) # correct for scientific notation
    Median = sprintf(myfmt2, signif(median(x, na.rm=T),ns2))
      if(grepl("e", Median)) Median = sprintf("%g", signif(median(x, na.rm=T), ns)) # correct for scientific notation
    Max = sprintf(myfmt1, signif(max(x, na.rm=T), ns))
      if(grepl("e", Max)) Max = sprintf("%g", signif(max(x, na.rm=T), ns)) # correct for scientific notation
    CV = sprintf("%.1f", sd(x, na.rm=T)/mean(x, na.rm = T)*100)
    se = sd(x, na.rm=T)/sqrt(nn)
    uci = mean(x, na.rm = na.rm) + qt(p=.975, df=nn-1)*se
    lci = mean(x, na.rm = na.rm) - qt(p=.975, df=nn-1)*se
    UCI = sprintf(myfmt2, signif(uci, ns2)) 
      if(grepl("e", UCI)) UCI = sprintf("%g", signif(uci, ns2)) # correct for scientific notation
    LCI = sprintf(myfmt2, signif(lci, ns2)) 
      if(grepl("e", LCI)) LCI = sprintf("%g", signif(lci, ns2)) # correct for scientific notation

        ## get rid of trailing decimal dots
    LCI = sub("(.*)\\.$","\\1",LCI)
    UCI = sub("(.*)\\.$","\\1",UCI)

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
