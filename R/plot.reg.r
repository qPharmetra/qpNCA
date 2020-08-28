#' Plot Regression Curves
#'
#' Plots regression curves for each set of records defined using \code{by}.
#' A log-linear plot will be made for each curve.
#'
#' If elimination half-life was estimated for that curve,
#' the following will be indicated in the plot:
#' * Cmax (Yellow, even if no half-life was estimated)
#' * points used in regression and resulting regression line (green)
#' * points excluded from regression (red crossed)
#' * estimate of elimination half-life and adjusted R-squared
#'
#' Input dataset:
#' * uncorrected dataset, used for half-life estimation
#' * dataset containing results of the half-life estimation
#'
#'
#' @param x input dataset name (if called within dplyr: .)
#' @param by column names in x indicating grouping variables
#' @param th file name of file with half-life estimation information for each curve
#' @param bloqvar variable name containing the BLOQ flag (0: no, 1: yes)
#' @param timevar variable name containing the sampling time
#' @param depvar variable name containing the dependent variable (e.g., concentration)
#' @param timelab X-axis label (default: "timevar")
#' @param deplab Y-axis label (default: "depvar")
#' @param exclvar variable name containing information about points to be excluded (these should have exclvar = 1)
#' @param plotdir directory where individual plot files will be saved
#' @param ... ignored
#' @import ggplot2
#'
#' @return (invisible) plotdir.  If the attribute 'plotdir' is empty, plots will be generated in standard output, otherwise plots will be saved as
#' PNG file in the designated directory.
#' @export
#' @importFrom dplyr left_join group_by_at first ungroup
#' @examples
#' example(est.thalf)
#' x %>% plot_reg(by = 'subject', th = th)
#'
plot_reg <- function(
   x,
   by = character(0),
   th = NA,
   bloqvar = "bloq",
   timevar = "tad",
   depvar = "dv",
   timelab = "timevar",
   deplab = "depvar",
   exclvar = NA,
   plotdir = NA,
   ...
  ){
  data_in = x %>% mutate(
    timevar = x[[timevar]],
    depvar = x[[depvar]],
    bloqvar = x[[bloqvar]],
    exclvar = x[[exclvar]]
  )

  if (!is.na(exclvar) &
      !(exclvar %in% names(x)))
    stop(paste("Exclusion variable", exclvar, "does not exist"), call. = F)

  data_in = data_in %>% mutate(excl = 0)
  if (!is.na(exclvar)) {
    data_in = data_in %>% mutate(excl = exclvar)
  }  # if exclvar exists, set excl to exclvar, else set to 0

  plot = left_join(data_in, th, by = by) %>%
    filter(bloqvar == 0) %>%
    mutate(
      start_conc = exp(log(intercept) - lambda_z * start_th),
      # intercept was already exponentiated in calc_thalf
      end_conc = exp(log(intercept) - lambda_z * end_th),
      # intercept was already exponentiated in calc_thalf
      thalf_txt = ifelse(
        !is.na(thalf),
        paste0("Half-life: ", round(thalf, 2), "  "),
        "Half-life not calculated  "
      ),
      radj_txt = ifelse(!is.na(thalf),
                        paste0(
                          "Adj. R-squared: ", round(adj.r.squared, 4), "  "
                        ), " ")
    ) %>%
    group_by_at(by) %>%
    mutate(
      max_conc = 10 ** (ceiling(log10(
        max(depvar, na.rm = T)
      ))),
      min_conc = 10 ** (floor(log10(
        min(depvar, na.rm = T)
      ))),
      cmax = max(depvar, na.rm = T),
      tmax = first(timevar[depvar == cmax])
    ) %>%
    ungroup() %>%
    filter(!is.na(depvar),!is.na(cmax))

  plots = plot %>%
    filter(depvar > 0) %>%
    group_by_at(by) %>%
    do(
      plots = ggplot(data = .) +

        geom_line(
          data = .,
          mapping = aes(x = timevar, y = depvar),
          color = "#144A90",
          size = 0.5
        ) +

        geom_point(
          data = .,
          mapping = aes(x = tmax, y = cmax),
          color = "gold3",
          size = 4
        ) +        # Cmax

        geom_point(
          data = .[.$timevar >= .$start_th & .$timevar <= .$end_th &
                     .$excl != 1, ],
          mapping = aes(x = timevar, y = depvar),
          size = 4,
          color = "#4CB54F"
        ) +          # points used in regression

        geom_point(
          data = .,
          mapping = aes(x = timevar, y = depvar),
          color = "#144A90",
          size = 2.5
        ) +   # curve

        geom_segment(
          data = .,
          x = .$start_th,
          xend = .$end_th,
          y = log10(.$start_conc),
          yend = log10(.$end_conc),
          # regression line
          color = "#4CB54F",
          linetype = "solid",
          size = 0.5
        ) +

        #           geom_point(data=.[.$excl==1,],
        #                      mapping=aes(x=timevar, y=depvar), size=4, color="#144A90") +

        geom_point(
          data = .[.$excl == 1, ],
          mapping = aes(x = timevar, y = depvar),
          shape = 4,
          size = 4,
          color = "red"
        ) +    # exclusions

        annotate(
          geom = "text",
          label = unique(.$thalf_txt),
          x = Inf,
          y = Inf,
          hjust = 1,
          vjust = 1.5,
          size = 3
        ) +

        annotate(
          geom = "text",
          label = unique(.$radj_txt),
          x = Inf,
          y = Inf,
          hjust = 1,
          vjust = 3,
          size = 3
        ) +

        annotate(
          geom = "text",
          label = "Cmax",
          x = unique(.$tmax),
          y = unique(.$cmax),
          hjust = 0.5,
          vjust = -1,
          color = "black",
          size = 3
        ) +

        annotate(
          geom = "text",
          label = "  Green data points were included in the elimination half-life estimation, red-crossed data points were excluded",
          x = -Inf,
          y = Inf,
          hjust = 0,
          vjust = 1.5,
          color = "black",
          size = 2.5
        )  +

        scale_y_log10(
          limits = c(unique(.$min_conc), unique(.$max_conc)),

          breaks = (10 ** (seq(
            log10(unique(.$min_conc)), log10(unique(.$max_conc))
          ))),
          labels = (10 ** (seq(
            log10(unique(.$min_conc)), log10(unique(.$max_conc))
          ))),
          expand = expand_scale(mult = c(0.05, .25))
        ) +
        xlab(timelab) +
        ylab(deplab) +
        theme(
          plot.margin = unit(c(5, 5, 0, 0), "mm"),
          #                 legend.direction="horizontal",legend.position="bottom", aspect.ratio=0.4,
          #                 legend.text = element_text(size=15),
          plot.title = element_text(hjust = 0, size = 10),
          panel.background = element_rect(fill = 'white', colour =
                                            "black"),
          panel.grid.major = element_line(colour = "grey"),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(colour = "black", size = 8),
          axis.text.y = element_text(colour = "black", size = 8)
        ) +
        ggtitle(titlefun(., by))
    )   %>% ungroup

  plots = plots %>% mutate(filename = paste0(filenamefun(., by), ".png"))

  if (!is.na(plotdir))  {
    if (file.exists(plotdir)) {
      mapply(
        ggsave,
        file = file.path(plotdir, plots$filename),
        plot = plots$plots,
        width = 18,
        height = 8,
        units = "cm"
      )
    } else {
      dir.create(plotdir)
      mapply(
        ggsave,
        file = file.path(plotdir, plots$filename),
        plot = plots$plots,
        width = 18,
        height = 8,
        units = "cm"
      )
    }
  }
  else {
    print(plots$plots)
  }
  invisible(plots$plots)
}

#' Create Title for Regression Plots
#'
#' Creates title for regression plots in plot_reg() using by-variables.
#'
#' @param x dataset containing concentration-time information of the current curve
#' @param by column names in x indicating grouping variables
#' @return character
titlefun <- function(x,by) {
  plottitle=""
  for (i in 1:length(by)) {
    plottitle <- paste0(plottitle,by[i],": ",unique(x[[by[i]]])," ")
  }
  plottitle=substr(plottitle,1,(nchar(plottitle)-1))
  return(plottitle)
}

#' Create File Name for Regression Plots.

#' Creates file name for regression plots (*.png) from by-variables in plot_reg function.
#' Also used by nca.sum().
#' @param x data.frame
#' @param by column names in x indicating grouping variables
#' @return character

filenamefun <- function(x,by) {
  filename=""
  for (i in 1:length(by)) {
    filename <- paste0(filename,by[i],"_",x[[by[i]]],"_")
  }
  filename=substr(filename,1,(nchar(filename)-1))
  return(filename)
}

