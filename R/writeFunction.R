#' WriteFunction
#'
#' Function to write the calculated analysis in a df form in the given path
#'
#' @param file file to be written
#' @param path path of the file
#' @param name name of the file
#' @return a xlsx file written
#' @examples
#' WriteFunction (file, path, name)
WriteFunction <- function(input.df, path, name){

  writexl::write_xlsx(input.df, path = paste0(path, name, ".xlsx"))

}
