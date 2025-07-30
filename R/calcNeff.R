#' calcNeff
#'
#' Function to calculate bottleneck effect size for each patient
#'
#' @param s number of mutations
#' @param mvec the list of AF - after
#' @param nvec the list of AF - before
#' @return  N-effective (bottleneck effect size)
#' @examples
#' Neff <- calcNeff(s, mvec, nvec)
calcNeff <- function(s, mvec, nvec){
  #N <- s/(2* sumi(KL(m|n) ) )
  N <- NULL
  KL_i_m_n <- 0
  for (nm in 1:s){
    mve <- mvec[nm] #ceiling(mvec[i] *1000)
    nve <- nvec[nm] #ceiling(nvec[i] *1000)

    f <- mve/nve
    if ( is.nan(f) ){      # both zero
      ff <- 0
    }else if ( is.infinite(f) ){ # before zero
      ff <- 1
    }
    else if (f == 0){
      ff <- 1                  # after zero
    } else{
      ff <- log(f)
    }
    #print(ff)
    KL_i_m_n <- KL_i_m_n + mve * ff
    #print(KL_i_m_n)
  }
  N <- s/ (2*KL_i_m_n)
  #N
  return(N)
}
