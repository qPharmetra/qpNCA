# ds = ncaTab
# un=c("",un[-1])
# nams = names(ds)
formatTab = function(myLst){
  # Format the subject-level tables and return Latex-ready tables
  ds = myLst$ncaTab
  un = c("", myLst$units[-1])
  nams = names(ds)
  ok = 2:length(un)
  un[ok] = paste("$(",sub("\\*","\\\\cdot ",  un[ok]),")$", sep = "")
  un[1] = Cs(Patient) ## put statistic at bottom of cell
  un[grepl("Accumulation", names(ds))] = "Ratio"
  un[names(ds)=="R\\SP{2}"] = "R\\SP{2}"
  un[names(ds)=="R\\SP{2}"] = "R\\SP{2}"
  ds = insert.row(ds,un,1)
  names(ds) = nams
  names(ds)[1] = "~"
  # Split "Accumulation Ratio" over two lines
  names(ds)[grep("Accumulation", names(ds))] = "Accumulation"
  return(ds)
}


