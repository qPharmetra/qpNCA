#' Summarize NCA parameters table
#'
#' @param ds
#' @param subjVar
#' @param carryAlong
#' @param omits
#' @param nsig
#' @param LOQ
#'
#' @return
#' @export
#'
#' @examples
nca.summary.table = function(ds, subjVar = "subjid", carryAlong=NULL, omits=NULL
                       , nsig=rep(3,ncol(ds)-1-length(carryAlong)), LOQ=.0001){
  ## Summarize  nca parameter table
  ##   ds = formatted nca parameter table (output from formatNCAtab)
  ##   subjVar = variable name of subject identifier
  ##   carryAlong = character vector of variable names to be included in output table; it expected that these
  ##                will be placed before teh parameter portion of the table.
  ##   omits = vector of subject identifiers (subjid) to be omitted from calculation of summary statistics
  ##   nsig = vector of default number of significant figures to report parameter values
  ##   LOQ = limit of quantitation for concentration variables
  ##
  names(ds)[grep(subjVar, names(ds))] = "subjid"  # rename subjVar to subjid for convenience

  #
  ds$subjid = ordered(ds$subjid, levels = unique(ds$subjid))

  stats = data.frame(sapply(names(ds)[names(ds) %nin% c("subjid", carryAlong)],
                            function(nm, dss=ds,oms=omits, nss=nsig, loq=LOQ){
                              ok = names(dss)==nm
                              xx=dss[dss$subjid %nin% oms, ok]
                              ns = nss[ok[(2+length(carryAlong)):ncol(dss)]]
                              # interogate nm to find whether it is a concentration or another parmeter
                              if(tolower(nm) %in% Cs(cmax, cmin, cavg)) {
                                return(nca.sumstat.conc(xx, LOQ=loq,ns))
                              } else {# other parameter, e.g., AUC, thalf, CL, etc.
                                return(nca.sumstat.params(xx, ns))
                              }
                            },
                            USE.NAMES=T
  ))

  if(is.null(carryAlong)){
    stats = cbind(subjid = rownames(stats), stats)
    names(stats) = names(ds)    # Some names get changed in lapply
    ds = ds[order(ds$subjid),]
    ds$subjid = as.character(ds$subjid)
  } else{
    blanks = data.frame(matrix(rep("", length(carryAlong)*nrow(stats)), nrow=nrow(stats)) )
    names(blanks) = carryAlong
    stats = cbind(subjid=rownames(stats), blanks, stats)
    names(stats) = names(ds)    # Some names get changed in lapply
    ds = ds[do.call(order, ds[,c(carryAlong, "subjid")]),]
    ds$subjid = as.character(ds$subjid)
  }
  return(rbind(ds, stats))
}
