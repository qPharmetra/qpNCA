calc.ctmax <- function(x,timevar="time",depvar="dv") {
              result=x %>% mutate(depvar=x[[depvar]],        # calculated dependent variable           (internal)
                                  timevar=x[[timevar]]       # calculated time variable                (internal)
              ) %>%
                filter(!is.na(depvar))            # na.rm=T in the max function does not work?

              if (dim(result)[1]==0) {            # if all concentrations are NA, set Cmax and Tmax to NA
                result=result%>%summarise(cmax=NA,
                                          tmax=NA)
              }
              else {
                result=result%>% summarise(cmax=max(depvar),
                                           tmax=first(timevar[depvar==cmax]) # there might be more than 1 timepoint with Cmax
              )
              }

              return(result)
            }
