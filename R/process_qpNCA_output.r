#' Processes list output from qpNCA function and provides an output with selected varibles in a data frame.
#'
#' @param data
#' @param covs
#' @param by
#' @param keepers
#'
#' @return
#' @export
process_qpNCA_output = function(data = ncaresults
                                , covs = c(DOSE, FORM)
                                , by = ID
                                , keepers = c(cmax, tmax, thalf, clast.pred,
                                              aucinf.pred, cl.f.pred, vz.f.pred)
)
{
  by= enquo(by)
  keepers = enquo(keepers)
  covs = enquo(covs)
  pkdat = data$pkpar %>% group_by(!!by) %>% select(!!keepers)
  covdat = data$covariates %>% group_by(!!by) %>% select(!!covs)
  out = left_join(pkdat, covdat)
  return(out)
}
