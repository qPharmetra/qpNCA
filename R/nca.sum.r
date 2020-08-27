#' Create Summary for Each Curve
#'
#' Creates a PDF summary for each curve.
#'
#' @param x data.frame
#' @param corrfile a correction table like the output of \code{\link{tab.corr}}
#' @param by column names in x indicating grouping variables
#' @param pdfdir path to directory for pdf ouput
#' @return (invisible) pdfdir
#' @importFrom rmarkdown render
#' @export
#' @importFrom dplyr group_by_at ungroup
#' @examples
#' example(qpNCA)
nca.sum <- function(x,corrfile, by=c("subject"),pdfdir=NA) {

  if (is.na(pdfdir)) invisible(pdfdir)
  if (!file.exists(pdfdir)) { dir.create(pdfdir) }

  x %>% # this step is needed to remove old pdf files, otherwise it will kind of append the content
    group_by_at(by) %>%
    do(x=if (file.exists(paste0(pdfdir,"/",filenamefun(.,by),".pdf"))) file.remove(paste0(pdfdir,"/",filenamefun(.,by),".pdf"))
    ) %>%
    ungroup

  x %>%
    group_by_at(by) %>%
    do(
      x=rmarkdown::render(
       input = "./create.summary.nca.pages.pdf.Rmd",
       output_format = "pdf_document",
       output_file = paste0(filenamefun(.,by),".pdf"),
       output_dir = pdfdir,
       quiet=T,
       clean=TRUE
    )
  ) %>%
  ungroup
  invisible(pdfdir)
}
