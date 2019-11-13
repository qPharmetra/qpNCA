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
plot.reg <- function(x,by=c("id"),th=NA,bloqvar="bloq",timevar="tad",depvar="dv",timelab="timevar",deplab="depvar",excl=NA,plotdir=NA) {

  data_in = x %>% mutate(timevar=x[[timevar]],
                         depvar=x[[depvar]],
                         bloqvar=x[[bloqvar]],
                         exclvar=x[[excl]])

  if (!is.na(excl) & !(excl %in% names(x))) stop(paste("Exclusion variable",excl,"does not exist"), call.=F)

  data_in = data_in %>% mutate(excl=0)
  if(!is.na(excl)) { data_in = data_in %>% mutate(excl=exclvar) }  # if exclvar exists, set excl to exclvar, else set to 0

  plot = left_join(data_in,th) %>%
    filter(bloqvar==0) %>%
    mutate(start_conc=exp(intercept-lambda_z*start_th),
           end_conc=exp(intercept-lambda_z*end_th),
           thalf_txt=ifelse(!is.na(thalf),
                            paste0("Half-life: ",round(thalf,2),"  "),"Half-life not calculated  "),
           radj_txt=ifelse(!is.na(thalf),
                           paste0("Adj. R-squared: ",round(adj.r.squared,4),"  ")," ")) %>%
    group_by_at(by) %>%
    mutate(max_conc=10**(ceiling(log10(max(depvar,na.rm=T)))),
           min_conc=10**(floor(log10(min(depvar,na.rm=T)))),
           cmax=max(depvar,na.rm=T),
           tmax=first(timevar[depvar==cmax])
    ) %>%
    ungroup() %>%
    filter(!is.na(depvar), !is.na(cmax))

  plots = plot %>% group_by_at(by) %>%
    do(  plots=ggplot(data=.) +

           geom_line(data=., mapping=aes(x=timevar, y=depvar), color=qp.blue, size=0.5) +

           geom_point(data=., mapping=aes(x=tmax, y=cmax), color="gold3", size=4) +        # Cmax

           geom_point(data=.[.$timevar>=.$start_th&.$timevar<=.$end_th&.$excl!=1,],
                      mapping=aes(x=timevar, y=depvar), size=4, color=qp.green) +          # points used in regression

           geom_point(data=., mapping=aes(x=timevar, y=depvar), color=qp.blue, size=2.5) +   # curve

           geom_segment(data=.,x=.$start_th,xend=.$end_th,
                        y=log10(.$start_conc), yend=log10(.$end_conc),                     # regression line
                        color=qp.green, linetype="solid", size=0.5) +

           #           geom_point(data=.[.$excl==1,],
           #                      mapping=aes(x=timevar, y=depvar), size=4, color=qp.blue) +

           geom_point(data=.[.$excl==1,],
                      mapping=aes(x=timevar, y=depvar), shape=4, size=4, color="red") +    # exclusions

           annotate(geom="text",label=unique(.$thalf_txt), x=Inf, y=Inf, hjust=1, vjust=1.5, size=3) +

           annotate(geom="text",label=unique(.$radj_txt), x=Inf, y=Inf, hjust=1, vjust=3, size=3) +

           annotate(geom="text", label="Cmax",x=unique(.$tmax), y=unique(.$cmax), hjust=0.5, vjust=-1, color="black", size=3) +

           annotate(geom="text", label="  Green data points were included in the elimination half-life estimation, red-crossed data points were excluded", x=-Inf, y=Inf, hjust=0, vjust=1.5, color="black", size=2.5)  +

           scale_y_log10(limits=c(unique(.$min_conc),unique(.$max_conc)),

                         breaks=(10**(seq(log10(unique(.$min_conc)),log10(unique(.$max_conc))))),
                         labels=(10**(seq(log10(unique(.$min_conc)),log10(unique(.$max_conc))))),
                         expand=expand_scale(mult = c(0.05, .25))
           ) +
           xlab(timelab) +
           ylab(deplab) +
           theme(plot.margin=unit(c(5,5,0,0),"mm"),
                 #                 legend.direction="horizontal",legend.position="bottom", aspect.ratio=0.4,
                 #                 legend.text = element_text(size=15),
                 plot.title =element_text(hjust = 0, size=10),
                 panel.background = element_rect(fill = 'white', colour="black"),
                 panel.grid.major = element_line(colour = "grey"),
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),
                 axis.text.x=element_text(colour="black", size=8),
                 axis.text.y=element_text(colour="black", size=8)) +
           ggtitle(titlefun(.,by))
    )   %>%ungroup

  plots=plots %>% mutate(filename=paste0(filenamefun(.,by),".png"))

  if (!is.na(plotdir))  {
    if (file.exists(plotdir)) {
      mapply(ggsave, file = file.path(plotdir,plots$filename), plot=plots$plots, width=18, height=8, units="cm")
    }
    else {
      dir.create(plotdir)
      mapply(ggsave, file = file.path(plotdir,plots$filename), plot=plots$plots, width=18, height=8, units="cm")
    }
  }
  else {
    print(plots$plots)
  }


}
