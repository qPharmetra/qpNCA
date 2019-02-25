ncaTabulate = function(ds, ns = nsignif, loq=.1, units=parUnits, omits = NULL){
  ncaTab = nca.summary.table(ds, subjVar="subjectc", nsig = ns,
                              LOQ=loq, omits=omits)
  un = units$unit[match(tolower(names(ncaTab)), units$param)]
  names(ncaTab) = units$dispName[match(tolower(names(ncaTab)), units$param)]
  ncaSum = ncaTab[(nrow(ncaTab)-9):nrow(ncaTab),]
  names(ncaTab)[1] = "Patient"
  names(ncaSum)[1] = "Statistic"
  return(list(ncaTab=ncaTab, ncaSum=ncaSum, units=un))
}