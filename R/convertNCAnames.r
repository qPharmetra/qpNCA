#' Format NCA parameter names for output in latex nca table
#'
#' @param ds
#'
#' @return
#' @export
convertNCAnames = function(ds){
  NCAParams = Cs(C0, Tlag, Tmax, Cmax, Cmax.D, Tlast, CmaxSS
                 , Cend, AUClast, AUC24, AUCtau, AUCtauSS, AUMClast, MRTlast
                 , Lambda.z, HL.Lambda.z, hl.eff
                 , AUCINF.obs, AUCINF.D.obs, Vz, Vz.F, Vz.obs, Vz.F.obs, CL.obs, CL.F.obs
                 , AUCINF.pred, AUCINF.D.pred, Vz.pred, Vz.F.pred, CL.pred, CL.F.pred
                 , Vss.obs, Vss.pred,Tau, Tmin, Cmin, CLss, CLss.F
                 , Accumulation.Index, ARCALC, ARCMAX, ARCTROUG
                 , R2, LAMZLL, LAMZNPT)
  LatexParams = c("C\\SB{0}", "t\\SB{lag}", "T\\SB{max}", "C\\SB{max}", "C\\SB{max}/D", "T\\SB{last}", "C\\SB{max,SS}"
                  , "C\\SB{24}", "AUC\\SB{last}", "AUC\\SB{24}", "AUC\\SB{tau}", "AUC\\SB{tau,SS}", "AUMC\\SB{last}", "MRT\\SB{last}"
                  , "$\\lambda$\\SB{z}", "t\\SB{1/2}", "t\\SB{1/2,eff}"
                  , "AUC\\SB{$\\infty$}", "AUC\\SB{$\\infty$}/D", "V\\SB{z}", "V\\SB{z}/F", "V\\SB{d}", "V\\SB{d}/F", "CL", "CL/F"
                  , "AUC\\SB{$\\infty$}", "AUC\\SB{$\\infty$}/D", "V\\SB{z}", "V\\SB{z}/F", "CL", "CL/F"
                  , "V\\SB{SS}", "V\\SB{SS}", "Tau", "t\\SB{min}", "C\\SB{trough}", "CL\\SB{SS}/F", "CL\\SB{SS}/F"
                  , "Accumulation Index", "R\\SB{AUC}", "R\\SB{max}", "R\\SB{trough}"
                  , "R\\SP{2}", "$\\lambda$\\SB{Z} T1", "$\\lambda$\\SB{Z} n\\SB{pt}")
  names(ds) = swap(names(ds), NCAParams, LatexParams)
  return(ds)
}
