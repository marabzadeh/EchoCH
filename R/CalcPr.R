#' calcPr
#'
#' Function to compute Fitness likelihood
#'
#' @param N predicted Neff
#' @param n the number of mutated allele counts - after
#' @param m the number of mutated allele counts - before
#' @return  Fitness likelihood
#' @examples
#' Pr <- calcPr(N,n,m)

calcPr <- function(N,n,m) {

  n <- ceiling(n*1000)
  m <- ceiling(m*1000)
  N <- ceiling(N)

  n1 <- N*n/(n+m)
  m1 <- N*m/(n+m)

  n2 <- ceiling(n1*100/N)
  m2 <- ceiling(m1*100/N)
  N1 <- 100

  f1 <- factorial(N1) / (factorial(m2)* factorial(N1-m2) )
  f2 <- (n2/N1)^m2
  f3 <- (1- n2/N1)^(N1-m2)

  p <- f1*f2*f3

  return (p)
}
