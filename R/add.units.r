#' Formatting for parameter tables
#'
#' @param tablst
#' @param rawtab
#'
#' @return
#' @export
#'
#' @examples
add.units = function(tablst, rawtab){
  # Some formatting for parameter tables
  un = c(rep("",2), tablst$units[-1])
  nams = names(rawtab)
  # Adjust the header for parameter units
  ok = un != ""
  un[ok] = paste("$(",sub("\\*","\\\\cdot ",  un[ok]),")$", sep = "")
  un[!ok] = nams[!ok] # put names of columns with no units at the bottom of the cell
  un[grepl("Accumulation", names(rawtab))] = "Ratio"
  nams[!ok] = "~"
  rawtab = insert.row(rawtab,un,1)
  names(rawtab) = nams
  return(rawtab)
}
