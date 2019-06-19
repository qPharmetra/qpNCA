#' Format NCA parameters table
#' @importFrom gdata trim
#' @param nca
#' @param carryAlong
#' @param parms
#' @param ns
#' @param uAdjust
#'
#' @return
#' @export
nca.format.table <- function(nca, carryAlong = NULL,
                        parms = Cs(Cmax,  Tmax, AUCall, AUCINF.pred, Cl.F.pred, HL.Lambda.z),
                        ns =    rep(3, length(parms)),
                        uAdjust = rep(1, length(parms)) )
{
  # process nca table without affecting the order of the observations
  # edit names to remove underscores
  names(nca) = gsub("_", ".", names(nca))
  parms = gsub("_", ".", parms)
  tmp = nca[, parms]
  tmp = data.frame(lapply(tmp, function(X) as.numeric(as.character(X))))
  nr = nrow(tmp)
  # format all of the variables in keepers
  tmp = data.frame(matrix(sapply(1:ncol(tmp), function(i, ds, adju, fmt,pm)  {
    ds[,i] = as.numeric(as.character(ds[,i]))*adju[i]
    # out = signifString(ds[,i], ns[i])
    out = sprintf("%g", signif(ds[,i], ns[i]))
    # BLOQ for 0's in concentrations; NA for other 0's
    if(pm[i] %in% Cs(Cmax,Cmin,Cavg)) {out[asNumeric(out)==0] = "BLOQ"
    } else {out[asNumeric(out)==0 | isMissing(out)] = "NC"}
    return(swap(out, "NA", "NC"))
  },
  ds=tmp, fmt=fmts, adju=uAdjust, pm=parms), nrow = nr))
  names(tmp) = parms
  if(!is.null(carryAlong))
  {
    nca = cbind(nca[, carryAlong], tmp)
    names(nca) = c(carryAlong, parms)
  } else{
    nca = tmp
    names(nca) = parms
  }

  # get rid of trailing .'s
  nca = data.frame(apply(nca, 2,
                         function(x) {
                           x = gdata::trim(as.character(x))
                           x = sub("(.*)\\.$","\\1",x)
                         }))
  return(nca)
}
