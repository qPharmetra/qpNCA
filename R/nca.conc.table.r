nca.conc.table = function(obs,                      # source dataframe
                        LOQ = .0001,              # LOQ set to a small value by default
                        sumVar = "conc",          # variable to be summarized
                        subjVar = "subjid",       # subject identifier
                        timeVar = "nom.time",     # sort variable
                        omits = NULL,             # vector of omitted subjects
                        times = NULL,             # vector of all times to be considered
                        carryAlong = NULL,        # Vector of variable names to be sorted and included in the table
                        # must be unique within subject
                        nsig = 3,                 # number of significant digits for the table
                        bloqCode="BLOQ"           # Code for BLOQ values
)
{ # obs = obs.tems  
  # Take a long thin dataset, and reshape to short and wide
  # Round numeric values of sumVar to nsig significant digits
  # Summarize the values in sumVar
  myfmt = paste("%#.", nsig, "g", sep="")
  nsubj = lunique(obs[,subjVar] )
  
  ind = whichNumeric(obs[, sumVar])
#  obs[, sumVar][ind] = sprintf("%g", signifString(as.numeric(obs[, sumVar][ind]),nsig))
  obs[, sumVar][ind] = signifString(as.numeric(obs[, sumVar][ind]),nsig)
  if(is.null(times)) times = obs[, timeVar] 
#  obs[, timeVar] = ordered(obs[, timeVar], levels=sunique(obs[,timeVar]))
  #names(obs)[grep(subjVar, names(obs))] = "subjid"  # rename subjVar to subjid for convenience
  
  # Expand observation table to include all times
  allTimes = expand.grid(sunique(obs[, subjVar]), sunique(times))
  names(allTimes) = c(subjVar, timeVar)
  
  if(is.null(carryAlong)) {
    obs = merge(obs, allTimes, by=intersect(names(obs), names(allTimes)), all=T)
  } else {
    subj = obs[!duplicated(obs[,subjVar]), c(subjVar, carryAlong)]
    allTimes = merge(allTimes[, names(allTimes) %nin% carryAlong], subj)
    obsNms = names(obs)[names(obs) %nin% carryAlong]
    obs = merge(obs[, obsNms], allTimes, by=intersect(obsNms, names(allTimes)), all=T)
  }  
  
  obs[, sumVar][isMissing(obs[,sumVar])] = "M"
  
  ## form the table
  obs = obs[order(obs[, subjVar], obs[,timeVar]),] 
  if(is.null(carryAlong)){
    tabl = reshape(obs[, c(subjVar, timeVar, sumVar)], direction="wide",
                   idvar=subjVar, timevar=timeVar, v.names=sumVar)
    tabl = tabl[order(tabl[,subjVar]),]
    if(nsubj>1) tabl = data.frame(apply(tabl, 2, function(x) sub("(.*)\\.$","\\1",x)))
    
  } else{
    tabl = reshape(obs[, c(subjVar, timeVar, sumVar, carryAlong)], direction="wide"
                   , idvar=c(subjVar, carryAlong), timevar=timeVar, v.names=sumVar)
    tabl = tabl[do.call(order, tabl[,c(carryAlong, subjVar)]),]
    if(nsubj>1) tabl = data.frame(apply(tabl, 2, function(x) sub("(.*)\\.$","\\1",x)))
  }
  
  stats = data.frame(
    lapply(tabl[tabl[,subjVar] %nin% omits, (2+length(carryAlong)):length(tabl)],
           FUN = nca.sumstat.conc, LOQ=LOQ, ns=nsig, bloqCode=bloqCode))
  if(is.null(carryAlong)){
    stats = cbind(rownames(stats), stats)
    names(stats)[1] = subjVar
  } else{
    blanks = data.frame(matrix(rep("", length(carryAlong)*nrow(stats), ncol=length(carryAlong)), nrow=nrow(stats)) )
    names(blanks) = carryAlong
    stats = cbind(rownames(stats), blanks, stats)
    names(stats)[1] = subjVar
  }
  
  ret = data.frame(rbind(as.matrix(tabl), as.matrix(stats)))
  names(ret) = names(tabl)
  return(ret)
}

