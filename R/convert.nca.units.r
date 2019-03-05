convert.nca.units = function(ncaList, parUnits = parUnits){
  nca = casefold(nca, F)
  parameters = casefold(parameters, F)
  nca[nca %nin% parameters] = "NA"
  un = swap(gsub("_",".",nca), gsub("_",".",parameters), units)
  gsub("NA","",un)
}

end.parameter.and.units = function(nca =  ncaSumTab, units = NCA.units)
{
  txs = nca
  ok = units != ""
  names(txs)[ok] = paste(rep(" ($",ncol(nca[,ok])),
         sub("\\*","\\\\cdot ",units[ok]),
        rep("$)}",ncol(nca[,ok])),sep = "")
  names(txs)[units == ""] = rep("}", sum(units == ""))
  return(names(txs))
}
#end.parameter.and.units(ncaSumTab,NCA.units)

convert.units = function(nca =  ncaSumTab, units = NCA.units)
{
  txs = nca
  ok = units != ""
  names(txs)[ok] = paste(rep(" ($",ncol(nca[,ok])),
         sub("\\*","\\\\cdot ",units[ok]),
        rep("$)",ncol(nca[,ok])),sep = "")
  names(txs)[units == ""] = rep("", sum(units == ""))
  return(names(txs))
}



