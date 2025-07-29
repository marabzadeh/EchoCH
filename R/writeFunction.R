#' WriteFunction
#'
#' Function to write the calculated analysis in a df form in the given path
#'
#' @param
#' @return
#' @examples
#' WriteFunction (file, path)

WriteFunction <- function(file, path){

   writexl::write_xlsx(file,paste0(path, file, ".xlsx"))

}
