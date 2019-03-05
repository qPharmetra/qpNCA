#' Insert blank rows
#'
#' @param rdata
#' @param what
#' @param pos
#'
#' @return
#' @export
#'
#' @examples
insert.row = function(rdata, what, pos = 2 )
{
  part1 = if(pos==1) rdata[0,] else rdata[seq(1:(pos-1)), ]
  part2 = rdata[pos:nrow(rdata), ]
  ndata = data.frame(rbind(part1, what,part2))
  ndata
}

