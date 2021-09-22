#' Clean up DLL
#'
#' As recommended in the "Compiled Code" chapter of the book
#' "R packages" (Version 2) by Hadley Wickham and Jenny Bryan.
#' @param libpath Argument
.onUnload <- function (libpath) {
  library.dynam.unload("igcop", libpath)
}
