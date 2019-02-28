calc.par.th <- function(x=par,th=th,cov=cov,dosevar="dose",factor=1, reg="sd", ss="n") {

  result=left_join(x,th) %>%
    left_join(cov) %>%
    mutate(dosevar=.[[dosevar]],
           reg=tolower(reg),
           ss=tolower(ss),
           clast.pred=exp(intercept-lambda_z*tlast),
           aucinf.obs=auclast+clast.obs/lambda_z,
           aucinf.pred=auclast+clast.pred/lambda_z,
           aumcinf.obs=aumclast+tlast*clast.obs/lambda_z + clast.obs/lambda_z^2,
           aumcinf.pred=aumclast+tlast*clast.pred/lambda_z + clast.pred/lambda_z^2,
           cl.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/aucinf.obs, dosevar*factor/auctau),
           cl.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/aucinf.pred, dosevar*factor/auctau),
           mrt.obs= ifelse(tolower(ss)=="n", aumcinf.obs/aucinf.obs,
                           (aumctau + tau*(aucinf.obs-auctau))/auctau),
           mrt.pred= ifelse(tolower(ss)=="n", aumcinf.pred/aucinf.pred,
                            (aumctau + tau*(aucinf.pred-auctau))/auctau),
           vz.f.obs=  ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.obs),
                             dosevar*factor/(lambda_z*auctau)),
           vz.f.pred= ifelse(tolower(ss)=="n",dosevar*factor/(lambda_z*aucinf.pred),NA),
           vss.obs= mrt.obs*cl.f.obs,
           vss.pred= mrt.pred*cl.f.pred,
           pctextr.obs=(clast.obs/lambda_z)/aucinf.obs*100,
           pctextr.pred=(clast.pred/lambda_z)/aucinf.pred*100,
           pctback.obs=area.back.extr/aucinf.obs*100,
           pctback.pred=area.back.extr/aucinf.pred*100
    ) %>%
    select(-dosevar)

  return(result)

}
