

#' Global Bayesian regression model
#'
#' Run the global Bayesian regression model and get predicted values with mean average error and exact CRPS.
#'
#' @param y_test vector of observations responses from the testing data.
#' @param X_test matrix of observations covariates from the testing data used for regression (rows = observations, columns = covariates, first column = unit vector).
#' @param X_train matrix of observations covariates from the training data and used for regression (rows = observations, columns = covariates, first column = unit vector).
#' @param y_train vector of observations responses from the training data.
#' @param beta_0 mean of the Gaussian prior.
#' @param V_0 covariance matrix of the Gaussian prior.
#' @param a_0 shape parameter of the inverse-gamma prior.
#' @param b_0 rate parameter of the inverse-gamma prior.
#' @param approx if approx is TRUE, approximate the CRPS instead of using the exact formula to compute it. Recommended to always use approx = FALSE, since
#' approximate method is for examples only.
#' @param nb_draws if approx is TRUE, number of draw from the estimated predictive distribution used to approximate the CRPS.
#' @param plot_hist_PIT if plot_hist_PIT is TRUE, plot the PIT histogram.
#' @returns data.frame with the pair of observed/predicted values with associated AE, CRPS and PIT metrics.
global_Bayes_reg=function(y_test,X_test,X_train,y_train,beta_0,V_0,a_0,b_0,approx=FALSE,nb_draws=200,plot_hist_PIT=FALSE){

  #Compute posterior distribution parameters
  n=length(y_train)
  V_n=solve(solve(V_0)+t(X_train)%*%X_train)
  beta_n=V_n%*%(solve(V_0)%*%beta_0+t(X_train)%*%y_train )
  a_n=a_0+n/2
  b_n=b_0+(t(y_train)%*%y_train +t(beta_0)%*%solve(V_0)%*%beta_0-t(beta_n)%*%solve(V_n)%*%beta_n)/2

  # If approx is TRUE, approximate the CRPS with approximate predictive distribution based ob Monte-Carlo. If FALSE (recommended), use the exact
  # student t predictive distribution to compute the CRPS.
  if(approx==TRUE){
    #Sample from posterior
    sig_sample=rep(0,nb_draws)
    beta_sample=matrix(nrow=p_reg+1,ncol=100)

    for (k in 1:100) {
      sig_sample[k]=rinvgamma(1,shape=a_n,rate=b_n)
      beta_sample[,k]=rmvnorm(1,beta_n,sig_sample[k]*V_n)
    }
    AE_vec=abs(y_test-rowMeans(X_test%*%beta_sample))


    predi_mat=matrix(nrow=length(y_test),ncol=nb_draws)
    CRPS_vec=rep(0,length(y_test))


    # Draw nb_draws density parameters to approximate the predictive distribution CRPS
    for (i in 1:length(y_test)) {
      for(j in 1:nb_draws){
        tir=sample(1:100,1)
        predi_mat[i,j]=rnorm(1,mean= beta_sample[,tir]%*%X_test[i,],sd=sqrt(sig_sample[tir]))
      }
      CRPS_vec[i]=CRPS_approx(predi_mat[i,],y_test[i])
      if(i%%1000==0){
        print(i)
      }
    }
    PIT_vec=pit_sample(y_test,predi_mat)
  }
  else{
    #Student's t predictive distribution CRPS
    CRPS_vec=y_pred=AE_vec=rep(0,length(y_test))
    for (i in 1:length(y_test)) {
      nu_t=2*a_n
      mu_t= X_test[i,]%*%beta_n
      Sigma_t=(b_n/a_n)*(1+X_test[i,]%*%V_n%*%X_test[i,])
      Sigma_t=sqrt(Sigma_t)
      CRPS_vec[i]=crps_t(y_test[i],df=nu_t,location = mu_t,scale =Sigma_t)
      AE_vec[i]=abs(y_test[i]-mu_t)
      if(i%%1000==0){
        print(i)
      }
    }

    # Compute the PIT values
    PIT_vec=pst(y_test, mu=mu_t,sigma=Sigma_t,nu=nu_t,log=FALSE,lower.tail = TRUE)
  }

  if(plot_hist_PIT==TRUE){
    hist(PIT_vec,main=paste("Global bayesian regression model PIT Histogram with nb_draws= ",nb_draws,sep = ""))
  }

  #Print and return the results
  print(c("MAE","mean_CRPS"));print(c(mean(AE_vec),mean(CRPS_vec)))
  return(data.frame(observed = y_test, predicted = X_test%*%beta_n, AE = AE_vec, CRPS = CRPS_vec, PIT = PIT_vec))

}



#' Local Bayesian regression model
#'
#' Run the local Bayesian regression model and get predicted values with mean average error and exact CRPS.
#'
#' @param y_test vector of observations responses from the testing data.
#' @param X_test matrix of observations covariates from the testing data used for regression (rows = observations, columns = covariates, first column = unit vector).
#' @param X_train matrix of observations covariates from the training data and used for regression (rows = observations, columns = covariates, first column = unit vector).
#' @param y_train vector of observations responses from the training data.
#' @param grid_train spacetime index data.frame/matrix with rows = training observations, first column = time points and second column = stations.
#' @param grid_test spacetime index data.frame/matrix with rows = testing observations, first column = time points and second column = stations.
#' @param beta_0 mean of the Gaussian prior.
#' @param V_0 covariance matrix of the Gaussian prior.
#' @param a_0 shape parameter of the inverse-gamma prior.
#' @param b_0 rate parameter of the inverse-gamma prior.
#' @param approx if approx is TRUE, approximate the CRPS instead of using the exact formula to compute it. Recommended to always use approx = FALSE, since
#' approximate method is for examples only.
#' @param nb_draws if approx is TRUE, number of draw from the estimated predictive distribution used to approximate the CRPS.
#' @param plot_hist_PIT if plot_hist_PIT is TRUE, plot the PIT histogram.
#' @param local_station If local_station is TRUE, run the space local model. Else, run the time local model.
#' @returns data.frame with pairs of observed/predicted values with associated AE, CRPS and PIT metrics.
local_Bayes_reg=function(y_test,X_test,X_train,y_train,grid_train,grid_test,beta_0,V_0,a_0,b_0,approx = FALSE,nb_draws=200,plot_hist_PIT = FALSE,local_station=TRUE){

  # Initialized the results vectors
  CRPS_vec = AE_vec = predi_vec = PIT_vec = rep(0, length(y_test))
  p_reg = ncol(X_train)

  if(local_station==TRUE){

    # At each station, compute the CRPS, MAE and PIT values based on the local bayes regression model
    for(j in sort(unique(grid_train[,2]))){

      # Indexes indicating observations coming from station j
      idx_j_train=which(grid_train[,2]==j)
      idx_j_test=which(grid_test[,2]==j)
      print(c("station",j))

      # Run locally the bayesian regression model and getting predictive values
      modele_reg_bayes_j=global_Bayes_reg(y_test[idx_j_test],X_test[idx_j_test,],X_train[idx_j_train,],
                                             y_train[idx_j_train],rep(0,p_reg),diag(p_reg), v, s, approx = approx, nb_draws = nb_draws,plot_hist_PIT = FALSE)
      CRPS_vec[idx_j_test]=modele_reg_bayes_j$CRPS
      AE_vec[idx_j_test]=modele_reg_bayes_j$AE
      predi_vec[idx_j_test]=modele_reg_bayes_j$predicted
      PIT_vec[idx_j_test]=modele_reg_bayes_j$PIT
    }
  }
  else{

    # At each time, compute the CRPS, MAE and PIT values based on the local bayes regression model
    for(j in sort(unique(grid_train[,1]))){

      # Indexes indicating observations coming from station j
      idx_j_train=which(grid_train[,1]==j)
      idx_j_test=which(grid_test[,1]==j)
      print(c("time",j))

      # Run locally the Bayesian regression model and getting predictive values
      modele_reg_bayes_j=global_Bayes_reg(y_test[idx_j_test],X_test[idx_j_test,],X_train[idx_j_train,],
                                          y_train[idx_j_train],rep(0,p_reg),diag(p_reg), v, s, approx = approx, nb_draws = nb_draws, plot_hist_PIT = FALSE)
      CRPS_vec[idx_j_test]=modele_reg_bayes_j$CRPS
      AE_vec[idx_j_test]=modele_reg_bayes_j$AE
      predi_vec[idx_j_test]=modele_reg_bayes_j$predicted
      PIT_vec[idx_j_test]=modele_reg_bayes_j$PIT
    }
  }

  Resultats=data.frame(observed = y_test, predicted = predi_vec, AE = AE_vec, CRPS_vec = CRPS_vec, PIT = PIT_vec)

  if(plot_hist_PIT==TRUE){
    hist(PIT_vec,main=paste("Local bayesian regression model PIT Histogram PIT with nb_draws = ",nb_draws,sep = ""))
  }

  #Print and return the results
  print(c("MAE","mean CRPS"));print(c(mean(AE_vec),mean(CRPS_vec)))
  return(Resultats)

}



