#' Mean CRPS function of the EMOS model.
#'
#' Computes the mean CRPS for the EMOS model based on numerical weather predictions (NWP), using the exact CRPS formula for the normal distribution.
#'
#'@param par Optimized parameters obtained on training data with BFGS method, for example.
#'@param x testing data NWP means.
#'@param s testing data NWP variances.
#'@param y testing data observed values.
#'@returns mean CRPS value.
CRPS_emos<- function(par,x,s,y) {
a<-par[1]
b<-par[2]
c<-par[3]
d<-par[4]
mean(sqrt(c+d*s)*(((y-(a+b*x))/(sqrt(c+d*s)))*(2*pnorm((y-(a+b*x))/(sqrt(c+d*s)))-1)+
                    2*dnorm((y-(a+b*x))/(sqrt(c+d*s)))-1/sqrt(pi)))
}

#' MAE function of the EMOS model.
#'
#' Computes the mean absolute error (MAE) for the EMOS model based on numerical weather predictions (NWP).
#'
#'@param par Optimized parameters obtained on training data with BFGS method, for example.
#'@param x testing data NWP means.
#'@param y testing data observed values.
#'@returns MAE value.
MAE_emos<-function(par,x,y){
mean(abs(y-(par[1]+par[2]*x)))
}



#' CRPS function of the EMOS model.
#'
#' Computes each observation CRPS with the mean CRPS for the EMOS model based on numerical weather predictions (NWP),
#'  using the exact CRPS formula for the normal distribution.
#'
#'@param par Optimized parameters obtained on training data with BFGS method, for example.
#'@param x testing data NWP means.
#'@param y testing data observed values.
#'@returns List with vector of CRPS and mean CRPS value.
CRPS_emos_list<- function(par,x,s,y) {
  a<-par[1]
  b<-par[2]
  c<-par[3]
  d<-par[4]
  CRPS_list=vector("list",2)
  CRPS_list[[1]]=sqrt(c+d*s)*(((y-(a+b*x))/(sqrt(c+d*s)))*(2*pnorm((y-(a+b*x))/(sqrt(c+d*s)))-1)+
                                2*dnorm((y-(a+b*x))/(sqrt(c+d*s)))-1/sqrt(pi))
  CRPS_list[[2]]=mean(sqrt(c+d*s)*(((y-(a+b*x))/(sqrt(c+d*s)))*(2*pnorm((y-(a+b*x))/(sqrt(c+d*s)))-1)+
                                     2*dnorm((y-(a+b*x))/(sqrt(c+d*s)))-1/sqrt(pi)))

  print(mean(sqrt(c+d*s)*(((y-(a+b*x))/(sqrt(c+d*s)))*(2*pnorm((y-(a+b*x))/(sqrt(c+d*s)))-1)+
                            2*dnorm((y-(a+b*x))/(sqrt(c+d*s)))-1/sqrt(pi))))
  return(CRPS_list)
}


#' AE function of the EMOS model.
#'
#' Computes each observation average error (AE) with the mean absolute error (MAE) for the EMOS model based on numerical weather predictions (NWP).
#'
#'@param par Optimized parameters obtained on training data with BFGS method, for example.
#'@param x testing data NWP means.
#'@param y testing data observed values.
#'@returns list with vector of AE and MAE value.
MAE_emos_list<-function(par,x,y){
  print(mean(abs(y-(par[1]+par[2]*x))))
  MAE_list=vector("list",2)
  MAE_list[[1]]=abs(y-(par[1]+par[2]*x))
  MAE_list[[2]]=mean(abs(y-(par[1]+par[2]*x)))
  return(MAE_list)

}


#' Global mean-EMOS model.
#'
#' Run the global mean-EMOS model with the BFGS method.
#'
#'@param df_train training data.frame/matrix with following columns: space and time indexes as columns "time" and "station", "x" and "s" columns as NWP means and variances
#' and "y" as observed training values.
#'@param df_test testing dataframe/matrix with the same columns as training data.
#'@returns list containing testing observations/predictions pairs with their associated CRPS/AE and PIT metrics and the mean CRPS and MAE of all testing data.
global_mean_EMOS<- function(df_train,df_test) {

  liste=vector("list",2)

  #Optimized parameter obtained with the  BFGS method (optim).
  param_opti=optim(par=c(1,1,1,1),fn=CRPS_emos,gr=NULL,method="L-BFGS-B",lower = c(-10000,-10000,0,0),x=df_train$x,y=df_train$y,s=df_train$s)
  a<-param_opti$par[1]
  b<-param_opti$par[2]
  c<-param_opti$par[3]
  d<-param_opti$par[4]

  x = df_test$x
  s = df_test$s
  y = df_test$y

  #Compute CRPS/MAE.
  CRPS_liste=sqrt(c+d*s)*(((y-(a+b*x))/(sqrt(c+d*s)))*(2*pnorm((y-(a+b*x))/(sqrt(c+d*s)))-1)+
                                2*dnorm((y-(a+b*x))/(sqrt(c+d*s)))-1/sqrt(pi))
  AE_liste=abs(y-(par[1]+par[2]*x))

  PIT_liste = pnorm(y, mean = a+b*x, sd=sqrt(c+d*s))

  liste[[1]]=cbind(df_test[,c("time","station")],observed = y, predicted = par[1]+par[2]*x, AE = AE_liste, CRPS = CRPS_liste, PIT = PIT_liste)
  liste[[2]]=c(MAE = mean(liste[[1]]$AE),mean_CRPS = mean(liste[[1]]$CRPS))

  names(liste) = c("table_metrics", "mean_metrics")

  #Print mean CRPS/MAE
  print(liste[[2]])
  return(liste)
}

#'  Local mean-EMOS model.
#'
#' Run the Local mean-EMOS model with the BFGS method. Here, local refers to having an individual models for each stations, meaning that all time points
#' associated with a station share the same model, opposed to the global model where all observations share the same model regardless of spacetime.
#'
#'@param df_train training dataframe/matrix with following columns: spacetime indexes as columns "time" and "station", "x" and "s" columns as NWP means and variances
#' and "y" as observed training values.
#'@param df_test testing dataframe/matrix with the same columns as training data.
#'@returns list containing testing observations/predictions pairs with their associated CRPS/AE and PIT metrics, space points mean CRPS/MAE over time points
#'and the mean CRPS and MAE of all testing data.
local_mean_EMOS<- function(df_train,df_test) {

  AE_CRPS_local_EMOS=matrix(nrow = length(unique(df_test$station)), ncol = 2)
  AE_CRPS_local_EMOS_long=matrix(nrow=nrow(df_test) ,ncol=2)
  PIT_EMOS_local=rep(0,length(unique(df_test$station)))

  liste=vector("list", 3)


  for(j in unique(df_train$station)){
    idx_j=which(df_train$station == j)
    idx_j_test=which(df_test$station == j)
    param_opti_emos=optim(par=c(1,1,1,1),fn=CRPS_emos,gr=NULL,method="L-BFGS-B",lower = c(-10000,-10000,0,0),
                          x=df_train$x[idx_j_train],y=df_train$y[idx_j_train],s=df_train$s[idx_j_train])

    AE_EMOS_j=MAE_emos_list(param_opti_emos$par,x=df_test$x[idx_j_test],y=df_test$y[idx_j_test])
    CRPS_EMOS_j=CRPS_emos_list(param_opti_emos$par,x=df_test$x[idx_j_test],s=df_test$s[idx_j_test],y=df_test$y[idx_j_test])

    AE_CRPS_local_EMOS[j,]=c(AE_EMOS_j[[2]],CRPS_EMOS_j[[2]])
    AE_CRPS_local_EMOS_long[idx_j_test,1]=AE_EMOS_j[[1]] ; AE_CRPS_local_EMOS_long[idx_j_test,2]=CRPS_EMOS_j[[2]]
    PIT_EMOS_local[idx_j_test]=pnorm(df_test$y[idx_j_test],mean=param_opti_emos$par[1]+param_opti_emos$par[2]*df_test$x[idx_j_test] ,sd=sqrt(param_opti_emos$par[3]+param_opti_emos$par[4]*df_test$s[idx_j_test]))
  }

  #Compute CRPS/MAE.
  liste[[1]]=cbind(df_test[,c("time","station")], AE = AE_CRPS_local_EMOS_long[,1], CRPS = AE_CRPS_local_EMOS_long[,2], PIT = PIT_EMOS_local)
  liste[[2]]=cbind(df_test[,'station'], AE = AE_CRPS_local_EMOS[,1], CRPS = AE_CRPS_local_EMOS[,2])
  liste[[3]]=c(MAE = mean(AE_CRPS_local_EMOS[,1]), mean_CRPS = mean(AE_CRPS_local_EMOS[,1]))

  names(liste) = c("full_table_metrics", "station_table_metrics", "mean_metrics")

  #Print mean CRPS/MAE
  print(liste[[3]])
  return(liste)
}



#' PIT histogram
#'
#' Creates a PIT histogram based on given PIT values and numbers of breaks wanted.
#'
#'@param PIT vector of PIT values.
#'@param nb_breaks number of histogram breaks desired.
#'
#'@returns PIT histogram plot.
hist_PIT=function(PIT,nb_breaks){
  h = hist(PIT,breaks=nb_breaks,plot=FALSE)
  h$density = h$counts/sum(h$counts)*100
  plot(h,freq=FALSE,main="",xlab="PIT",ylab="Frequency (%)")
}




