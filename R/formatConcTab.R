#' Format raw concentrations table to latex suitable for qP NCA reports
#'
#' @param ds
#' @param dsTimes
#' @param nsig
#' @param omits
#'
#' @return
#' @export
formatConcTab = function(ds, dsTimes=NULL, nsig=3, omits = NULL){
  lloq = 0.075
  ds$dv[ds$dv==0] = "BQL"
  if(is.null(dsTimes)) dsTimes = sunique(ds$ntad)

  temp =   nca.conc.table(obs = ds, subjVar = "subjectc"
                          , timeVar="ntad", sumVar = "dv", LOQ=lloq
                          , bloqCode="BQL", nsig=nsig
                          , omits=omits
                          , times = dsTimes
  )
  # remove geometric mean from summary table
  temp = temp[-grep("Geometric", temp$subjectc),]

  # Adjust names
  names(temp) = c( "Patient", "Pre-dose", paste(dsTimes[2:length(dsTimes)], "h"))
  return(temp)

}

formatConcTabOth = function(ds, dsTimes=NULL, nsig=3, omits = NULL){
  lloq = 0.075
  ds$dv[ds$dv==0] = "BQL"
  if(is.null(dsTimes)) dsTimes = sunique(ds$ntad)

  temp =   nca.conc.table(obs = ds, subjVar = "subjectc"
                          , timeVar="timepoint", sumVar = "dv", LOQ=lloq
                          , bloqCode="BQL", nsig=nsig
                          , omits=omits
                          , times = dsTimes
  )
  # remove geometric mean from summary table
  temp = temp[-grep("Geometric", temp$subjectc),]

  # Adjust names
  names(temp) = c( "Patient", paste(dsTimes[1:length(dsTimes)], ""))
  return(temp)

}
