#' L2 distance
#'
#' Compute the L2 or euclidean distance between two numeric vectors.
#'
#' @param x a vector.
#' @param y a vector.
#' @returns a real number.
#' @examples
#' L2(c(-1,-1), c(1,2))
L2 <- function(x,y) {
  sqrt(sum((x-y)^2))
}



#' Threshold indicator function
#'
#' Indicator function that gives 1 if x <= y, 0 if not.
#'
#' @param x a numeric value.
#' @param y a numeric value.
#' @returns 1 if x <= y, 0 if not.
thresh_ind <- function(x,y) {
  if(x<=y){
    x=1
  }
  else{
    x=0
  }
  return(x)
}


#' Sinusoidal embedding
#'
#' Create a sinusoidal embedding from a vector of time points and a period parameter.
#'
#' @param x a vector of time points.
#' @param period period of the sin function.
#' @returns a numeric vector with the sin embedding.
sin_embedding<-function(x,period){
  sin(x/period*2*pi)
}


#' Importance
#'
#' Compute covariates importance based on CRPS.
#'
#' @param CRPS_FULL CRPS computed from the full model.
#' @param CRPS_shuffled CRPS computed from the shuffled model.
#' @returns Importance numeric value.
Importance <- function(CRPS_FULL,CRPS_shuffled) {
    return(mean(CRPS_shuffled-CRPS_FULL))
   }
