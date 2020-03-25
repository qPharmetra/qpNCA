#' Create summarising PDF for each curve
#'
#' @param x
#' @param corrfile
#' @param by
#' @param plotdir
#' @param pdfdir
#'
#' @return
#' @export
#'
#' @examples
nca.sum <- function(x,corrfile=corrtab,by=c("subject"),plotdir=NA,pdfdir=NA) {
  
  if (is.na(pdfdir)) return()
  if (!file.exists(pdfdir)) { dir.create(pdfdir) }
  
  x %>% # this step is needed to remove old pdf files, otherwise it will kind of append the content
    group_by_at(by) %>%
    do(x=if (file.exists(paste0(pdfdir,"/",filenamefun(.,by),".pdf"))) file.remove(paste0(pdfdir,"/",filenamefun(.,by),".pdf"))
    ) %>%
    ungroup
  
  x %>% 
    group_by_at(by) %>%
    do(x=rmarkdown::render(input = "./create.summary.nca.pages.pdf.Rmd", 
                           output_format = "pdf_document",
                           output_file = paste0(filenamefun(.,by),".pdf"),
                           output_dir = pdfdir,
                           quiet=T,
                           clean=TRUE)
    ) %>%
    ungroup
  
}
