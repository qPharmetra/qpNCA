#' A log-linear plot will be made for each curve. If elimination half-life was estimated for that curve,
#' the following will be indicated in the plot:
#'
#' - Cmax (Yellow, even if no half-life was estimated)
#' - points used in regression and resulting regression line (green)
#' - points excluded from regression (red crossed)
#' - estimate of elimination half-life and adjusted R-squared
#'
#' Input dataset:
#' - uncorrected dataset, used for half-life estimation
#' - dataset containing results of the half-life estimation
#'
#' USAGE:
#' plot.reg(x,by=c("id","period"),th=th,bloqvar="bloq",timevar="tad",depvar="dv",excl="excl_th",savedir=)
#'
#' @param x input dataset name (if called within dplyr: .)
#' @param by by-variable(s), e.g. c("subject","day")
#' @param th file name of file with half-life estimation information for each curve
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param timevar variable name containing the sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param timelab X-axis label (default: "timevar")
#' @param deplab Y-axis label (default: "depvar")
#' @param excl variable name containing information about points to be excluded (these should have <excl>=1)
#' @param savedir folder where individual plot files will be saved
#'
#' @return If the attribute 'savedir'is empty, plots will be generated in standard output, otherwise plots will be saved as
#' PNG file in the designated folder
#' PLOT.REG: Plot regression curves
#'
#' @export
plot.reg <- function(x,by=c("id"),th=th,bloqvar="bloq",timevar="tad",depvar="dv",timelab="timevar",deplab="depvar",excl=NA,savedir=NA) {

  #define function to print the plot title based on by-variables
  titlefun <- function(df,by) {
    plottitle=""
    for (i in 1:length(by)) {
      plottitle <- paste0(plottitle,by[i],": ",unique(df[[by[i]]])," ")
    }
    plottitle=substr(plottitle,1,(nchar(plottitle)-1))
    return(plottitle)
  }
  # end definition

  # define function to create filename based on by variables
  filenamefun <- function(df,by) {
    filename=""
    for (i in 1:length(by)) {
      filename <- paste0(filename,by[i],"_",unique(df[[by[i]]]),"_")
    }
    filename=substr(filename,1,(nchar(filename)-1))
    return(filename)
  }
  # end definition

  data_in = x %>% mutate(timevar=x[[timevar]],
                         depvar=x[[depvar]],
                         bloqvar=x[[bloqvar]])

  if   (!is.na(excl)) { data_in = data_in %>% mutate(excl=x[[excl]]) }
  else                { data_in = data_in %>% mutate(excl=0) }


  plot = left_join(data_in,th) %>%
    filter(bloqvar==0) %>%
    mutate(start_conc=exp(intercept-lambda_z*start_th),
           end_conc=exp(intercept-lambda_z*end_th),
           thalf_txt=ifelse(!is.na(thalf),
                            paste0("Half-life: ",round(thalf,2),"  "),"Half-life not calculated  "),
           radj_txt=ifelse(!is.na(thalf),
                           paste0("Adj. R-squared: ",round(adj.r.squared,4),"  ")," ")) %>%
    group_by_(by) %>%
    mutate(max_conc=10**(ceiling(log10(max(depvar,na.rm=T)))),
           min_conc=10**(floor(log10(min(depvar,na.rm=T)))),
           cmax=max(depvar),
           tmax=first(timevar[depvar==cmax])
    ) %>%
    ungroup() %>%
    filter(!is.na(depvar), !is.na(cmax))

  plots = plot %>% group_by_(by) %>%
    do(  plots=ggplot(data=.) +

           geom_line(data=., mapping=aes(x=timevar, y=depvar), color=qp.blue) +

           geom_point(data=., mapping=aes(x=tmax, y=cmax), color="gold3", size=6) +

           geom_point(data=., mapping=aes(x=timevar, y=depvar), color=qp.blue, size=4) +

           geom_segment(data=.,x=.$start_th,xend=.$end_th,
                        y=log10(.$start_conc), yend=log10(.$end_conc),
                        color=qp.green, linetype="solid", size=1) +

           geom_point(data=.[.$timevar>=.$start_th&.$timevar<=.$end_th,],
                      mapping=aes(x=timevar, y=depvar), size=4, color=qp.green) +

           geom_point(data=.[.$excl==1,],
                      mapping=aes(x=timevar, y=depvar), size=4, color=qp.blue) +

           geom_point(data=.[.$excl==1,],
                      mapping=aes(x=timevar, y=depvar), shape=4, size=6, color="red") +

           annotate(geom="text",label=unique(.$thalf_txt), x=Inf, y=Inf, hjust=1, vjust=2) +

           annotate(geom="text",label=unique(.$radj_txt), x=Inf, y=Inf, hjust=1, vjust=4) +

           scale_y_log10() +
           xlab(timelab) +
           ylab(deplab) +
           theme(legend.direction="horizontal",legend.position="bottom", aspect.ratio=0.4,
                 legend.text = element_text(size=15),
                 plot.title =element_text(hjust = 0, size=15),
                 panel.background = element_rect(fill = 'white', colour="black"),
                 panel.grid.major = element_line(colour = "grey"),
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
                 axis.text.x=element_text(colour="black", size=15),
                 axis.text.y=element_text(colour="black", size=15)) +
           ggtitle(titlefun(.,by))

    ) %>% ungroup

  #  plots$plots

  filename=paste0(filenamefun(plot,by),".png")
  if (!is.na(savedir))  {
    if (file.exists(savedir)) {
      mapply(ggsave, file = file.path(savedir,filename), plot=plots$plots)
    }
    else {
      dir.create(savedir)
      mapply(ggsave, file = file.path(savedir,filename), plot=plots$plots)
    }
  }
  else {
    print(plots$plots)
  }

}
