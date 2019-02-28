### _____________________________________________________________________________________________
### CLIENT:		qPharmetra 
### PROJECT:	qP PROJECT
###
### PROGRAM:	table.summarize.r
### LANGUAGE:	R 2.11.1
### PURPOSE:	creates standard (demographics) tables
### AUTHOR:		K. Prins
### DATE:     July, Aug 2010
### NOTES:		
### Copyright 2010, qPharmetra LLC, all rights reserved 
### _____________________________________________________________________________________________

examples = F

## mean (SE)
conDataFun1 = function(y, nSignif) paste(signif(mean(y),nSignif),
  " (",signif(sqrt(var(y)),nSignif-1),")", sep="")                                                    
## median (range)
conDataFun2 = function(y, nSignif) paste(signif(median(y),nSignif),
  " (", signif(min(y),nSignif-1),
  " - ",signif(max(y),nSignif-1), ")",sep="")
## mean (95% CI)
conDataFun3 = function(y, nSignif)
{
  se = sqrt(var(y))/sqrt(length(y)) 
  ci = qt(c(0.975, 0.025), length(y)) * se + mean(y)
  ci = ci[order(ci)]
  
  return(
    paste(signif(mean(y),nSignif),  " (", 
          signif(ci[1],nSignif-1),  " - ",
          signif(ci[2],nSignif-1),  ")",
          sep="")
        )
}

## count(%)
catDataFun = function(y){paste(table(y), " (",round(100*table(y)/sum(table(y))),"%)", sep="")}                   

## full.names reverts to standard full names based on qP standardized names
## it is some sort of inverse of the "standardize.names" functionality
full.names = function(x)
{
  exchange.values(x, 
    Cs(wt,ht,bmi,age,sex,crcl,race),
    c("Weight (kg)","Height (cm)","BMI (kg/m2)","Age (yr)","Gender","Creatinine Clearance (L/min)","Race")
    )      
}
      
## tabStats determines is variable is continuous or categorical and returns the statistics
tabStats = function(x, BY, nSignif = 3, conFunc1 = conDataFun1, conFunc2 = conDataFun2, catFunc = catDataFun, parName)
{
  if(missing(parName))
  {
    parName = deparse((match.call()[2]))
    parName = substring(parName,1, (nchar(parName)-2))
  }
  x = c(x, x)
  BY = lapply(BY, function(by) c(as.character(by), rep("All", length(by))))
  if(is.factor(x) | is.character(x))
  {
    if(is.character(x)) x = as.factor(x)
    tmp = data.frame((t(aggregate(x, by = BY, FUN = catFunc))))
    row.names(tmp)[(length(names(BY))+1) : nrow(tmp)] = levels(x)
  }
  if(is.numeric(x))
  {
    tmp1 = data.frame(t(aggregate(x, by = BY, FUN = function(y, nSignif) 
      conFunc1(y, nSignif = nSignif), nSignif = nSignif)))
    tmp2 = data.frame(t(aggregate(x, by = BY, FUN = function(y, nSignif) 
      conFunc2(y, nSignif = nSignif), nSignif = nSignif)))
    tmp = rbind(tmp1, tmp2[2,])
  }
  names(tmp) = levels(as.factor(unlist(BY)))
  tmp2 = as.data.frame(tmp[,1])
  names(tmp2) = "parameter"
  tmp2$parameter = row.names(tmp)
  tmp2$parameter[1] = parName
  if(is.numeric(x)) tmp2$parameter = c(parName, "Mean (SD)", "Median (range)")
  tmp = cbind(tmp2, tmp)
  row.names(tmp) = 1:nrow(tmp)
  for(i in 1:ncol(tmp)) tmp[,i] = as.character(tmp[,i])
  tmp[1, ] = c(parName, rep("", length(2:ncol(tmp))))
  return(tmp)
}

## test it out
examples = F
if(examples){
ok = duplicated(pkpdData$id) == F
tmp = tabStats(x=pkpdData$race[ok], BY=list(dose = pkpdData$dose[ok]))
tabStats(x=pkpdData$race[ok], BY=list(dose = pkpdData$dose[ok]), parName = "Race")
tabStats(pkpdData$sex[ok], list(dose = pkpdData$dose[ok]))
tabStats(x=pkpdData$wt, BY = list(dose = pkpdData$dose))
tabStats(x=pkpdData$bmi, BY = list(dose = pkpdData$dose))
}

## tabSummarize does tabStats using a formula and the target data set
tabSummarize = function(formula, data, nSignif = 3, extra.blank.line = FALSE)
{
  allX = all.vars(getResponseFormula(formula)[[2]])
  allY = all.vars(getCovariateFormula(formula)[[2]])
  BY = lapply(1:length(allX), function(x, allX, data) eval(as.name(allX[[x]]), data), 
    data = data, allX = allX)
  names(BY) = allX

  YYY = lapply(allY, function(yyy, data) eval(as.name(yyy), data), data = data)
  names(YYY) = allY
  theData = do.call("rbind",lapply(1:length(YYY), 
    function(z, YYY, BY, extra.blank.line,nSignif) 
    {                                          
      stats = tabStats(x = YYY[[z]], BY = BY, parName = names(YYY)[z], nSignif = nSignif)
      if(extra.blank.line == TRUE) 
      {
        EBL = stats[1,]
        EBL[1,] = rep("", ncol(stats))
        stats = rbind(EBL, stats)
      }
      return(stats)
    }, YYY = YYY, BY = BY, extra.blank.line = extra.blank.line, nSignif = nSignif))
  row.names(theData) = 1 : nrow(theData)
  
  theData$parameter = full.names(theData$parameter)
  names.order = as.character(unique(eval(as.name(allX[1]), data))) ## sorting properly 
  theData = theData[, c("parameter", names.order, "All")]

   ## insert (N=xxx)
  ndf = theData[1,]
  NNN = tapply(YYY[[1]], BY, length)
  NNN = NNN[names.order]
  ndf[1,] = c("", paste("(N=",c(as.numeric(NNN), length(YYY[[1]])), ")", sep = ""))#, length(YYY[[1]]))
  theData = rbind(ndf, theData)
  return(theData)
}
#exportData(demoTable)

if(examples)
{
## validate that the right stats are spat out by tabSummarize()
## you will need the pkpdData object
options(width = 150)
ok = duplicated(pkpdData$id) == F
tabSummarize(formula = dose ~ race + wt + bmi + sex, data = pkpdData[ok, ], nSignif = 3)
round(table(pkpdData$race[ok], pkpdData$dose[ok]) / apply(table(pkpdData$race[ok], pkpdData$dose[ok]), 2, sum) * 100)
## OK

## demographics by treatment 
tabSummarize(formula = trt ~ race + wt + bmi + sex, data = pkpdData[ok, ], nSignif = 3)
## demographics by treatment, when it is an ordered factor 
pkpdData$trt.o = ordered(pkpdData$trt, levels = unique(pkpdData$trt))
levels(pkpdData$trt.o)
demoTable = tabSummarize(formula = trt.o  ~ race + wt + bmi + sex, data = pkpdData[ok,], nSignif = 3)
print(demoTable)
demoTable = tabSummarize(formula = trt.o  ~ race + wt + bmi + sex, data = pkpdData[ok,], nSignif = 3, extra.blank.line = TRUE)
print(demoTable)
pkpdData$trt.o = NULL

} # end if(examples)

###
### === end of script ============
###
