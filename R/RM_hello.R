# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' @title hello
#' @description Print a hello statement
#' @param x The name of the person to address
#' @return The output from \code{\link{print}}
#' @export
#' @examples
#' hello("John")

hello <- function(x) {
  print(paste("Hello ", x))
}
