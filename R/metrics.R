#' Approximate CRPS
#'
#' Compute the approximate CRPS from a sample of prediction based on the generalized quantile function.
#'
#' @param y_sample sample of predictions for the observation.
#' @param y_obs observation value.
#' @returns numeric CRPS value.
CRPS_approx<- function(y_sample,y_obs) {
m=length(y_sample)
y_ordre=sort(y_sample,decreasing = FALSE)
somme=0
for (i in 1:m) {
  somme=somme+(y_ordre[i]-y_obs)*(m*(y_obs<y_ordre[i])-i+0.5)
}
valeur=(2*somme)/(m^2)
return(valeur)
}

#' Mean squared error (MSE)
#'
#' Compute the mean squared error (MSE) or the mean absolute error (MAE).
#'
#' @param y_obs observation value.
#' @param y_estimated estimated value for y_obs.
#' @param absolute compute the mean absolute error instead (TRUE = mean absolute error, FALSE = mean squared error).
#' @returns numeric average value
MSE<- function(y_estimated,y_obs,absolute=FALSE){
  somme=0

  for (i in 1:length(y_obs)) {

    if(absolute==TRUE){
      somme=somme+abs(y_obs[i]-y_estimated[i])
    }else{
      somme=somme+(y_obs[i]-y_estimated[i])^2
    }
  }
  return(somme/length(y_obs))
}

#' Multivariate student's t density
#'
#' Evaluate the multivariate student's density at one point given the distribution parameters.
#'
#' @param x a numeric vector.
#' @param v degrees of freedom ( v > 0).
#' @param mu mean parameter vector.
#' @param Sig scale parameter matrix.
#' @param log_d compute the density on the log-scale if log_d is TRUE.
#' @returns  numeric density value.
multivariate_t_density <- function(x,v,mu,Sig,log_d=TRUE) {
  p=dim(Sig)[1]
  dens=lgamma((v+p)/2)-((v+p)/2)*log(1+(t(x-mu)%*%solve(Sig)%*%(x-mu))/v)-(lgamma(v/2)+(p/2)*log(v*pi)+(1/2)*log(det(Sig)))
  if(log_d==FALSE){
    dens=exp(dens)
  }
  return(dens)
}

#' Empirical Variance
#'
#' Compute the empirical variance from vectors of observations.
#'
#' @param X data.frame or matrix such that rows=observations and columns=variables.
#' @returns covariance matrix.
Empirical_Variance <- function(X) {
  somme=0
  X_mean=rowMeans(X)
  for (q in 1:nrow(X)) {
    somme=somme+t(t(X[q,]-X_mean))%*%t(X[q,]-X_mean)
  }
  return(somme/nrow(X))
}

#' Clusters overlap indexes
#'
#' Function that computes the homogeneity, the separation and the similarity index in clustered data.
#'
#' @param X matrix or data.frame with columns = variables and rows=observations.
#' @param c vector of clusters indicator values.
#' @returns a data.frame containing the metrics.
overlap_index <- function(X,c) {

  nb=length(unique(c))
  n=dim(X)[1];p=dim(X)[2]
  mean_X=matrix(data=NA,nrow = nb,ncol=dim(X)[2])
  taille<-NULL;y<-NULL

  for (i in 1:nb) {
    if(sum(c==i)==1){
      mean_X[i,]=mean(X[which(c==i),])
    }
    else{
      mean_X[i,]=colMeans(X[which(c==i),] )
    }
    taille[i]=sum(c==i)
  }

  somme=0
  for (i in 1:n) {
    somme=somme+L2(X[i,],mean_X[c[i],])
  }

  Homogeneity=somme/n
  somme=0
  somme2=0
  for (i in 1:(nb-1)) {
    for (j in (i+1):nb) {
      somme=somme+taille[i]*taille[j]*L2(mean_X[i,],mean_X[j,])
      somme2=somme2+taille[i]*taille[j]
    }
  }
  Separation=somme/(somme2)
  ratio=1-(n/(n-1))*(Homogeneity/(Homogeneity+Separation))
  return( data.frame(Homogeneity=Homogeneity,Separation=Separation,ratio=ratio)
  )
}

#' Exponential correlation matrix
#'
#' Produce an exponential correlation matrix based on a distance matrix.
#'
#' @param distance An JxJ matrix of distances between the J stations.
#' @param lambda length-scale or frequency parameter.
#' @param jitter add a jitter to the correlation matrix diagonal if jitter = TRUE.
#' @param jitter_value jitter value.
#' @returns a real JxJ correlation matrix.
#' @examples
#' #Generate a random set of 100 stations coordinates for the example's purpose
#' coord = data.frame(latitude = runif(100,-25,25), longitude = runif(100, -25, 25))
#' #Compute the exponential correlation matrix with fixed lambda and added jitter
#' corr_exp(coord, 0.5, jitter = FALSE, 0.001)
corr_exp=function(distance, lambda, jitter=FALSE, jitter_value=0.001){
  mat=matrix(data = NA,nrow=nrow(distance),ncol=ncol(distance))
  for (i in 1:nrow(distance)) {
    for (j in 1:nrow(distance)) {
      mat[i,j]=exp(-(distance[i,j]^2)/lambda^2)
    }
    if(jitter==TRUE){
      mat[i,i]=mat[i,i]+jitter_value
    }
  }
  return(mat)
}
