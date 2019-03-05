#' Format raw individual concentrations table for report
#'
#' @param obs dataset
#' @param LOQ value of LOQ
#' @param sumVar concentration/DV variable
#' @param subjVar subject or id variable
#' @param timeVar time or visit variable
#' @param omits any subjects to be ommitted
#' @param times nominal time
#' @param carryAlong any variables to be carried in the table
#' @param nsig number of significant digits
#' @param bloqCode variable name for BLQs are flagged
#'
#' @return report ready individual concentrations table with summary statistics
#' @example
#' library(dplyr)
#' get dataset ready
#'NTAD <- c(0,0.3,0.5,1,2,4,5,7,9,12,24)
#'Theoph1 <- Theoph %>%
#'  mutate(NTAD=metrumrg::snap(Time, NTAD)) %>%
#'  mutate(Subject=as.numeric(as.character(Subject)), #converting from factor to numeric
#'         BQL = ifelse(conc<=0.25, 1, 0),                #just adding few BLQs to demonstrate functionality
#'        conc= ifelse(conc<=0.25, NA, conc))            #just adding few BLQs to demonstrate functionality
#'
#' #Tabulate concentrations vs time for each subject and calculate summary stats of concentrations per timepoint (can use visit instead of time)
#'ConcTab = nca.conc.table(Theoph1,
#'                         sumVar = "conc",
#'                         subjVar = "Subject",
#'                         timeVar = "NTAD",
#'                         LOQ = 0.250,               #just adding to demonstrate functionality
#'                         bloqCode = "BQL",
#'                         nsig=3)
#'
#' @export
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
  obs[, sumVar][ind] = formatted.signif(as.numeric(obs[, sumVar][ind]),nsig)
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
  obs[, sumVar][obs[,bloqCode]==1] = bloqCode

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

