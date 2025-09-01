#' Simulate clustered data using MixSim.
#'
#' Simulate clustered data using MixSim by specifying the average and/or max overlap parameters. Given the chosen number of clusters,
#' stations, time points and dimensions for the covariates, simulate clustered data with covariate and response value, with chosen added spatiotemporal random effects.
#' Also generates random coordinates on an uniform square.
#'@param max_overlap maximum overlap used for the controlling overlap.
#'@param bar_overlap average overlap used for the controlling overlap.
#'@param lim_x specific boundaries on the space of each covariate dimension. By default (0,1).
#'@param K number of stations/space points.
#'@param N number of time points.
#'@param nb_clusters number of clusters.
#'@param p covariates dimension size.
#'@param tau_time variance parameter for the random time effects.
#'@param varphi correlation parameter for the random time effects.
#'@param tau_space variance parameter for the spatial random effects.
#'@param station_by_time how to structure the long data grid. If TRUE, list each station for time 1, then each station for time 2, etc...
#' If FALSE, do the opposite.
#'@param distribution "N" for normally distributed response, "LN" for a log-normal distribution instead.
#'@param sigma_2 response variance parameter. Used to control noise in the response value.
#'@param add_space if TRUE add spatial random effects.
#'@param add_time if TRUE, add time random effects.
#'@returns a list containing simulated training and test data, generated coordinates, distance matrix W and parameters used to generate the observations.
Simulation_MixSim <- function(max_overlap,bar_overlap,lim_x,K,N,nb_clusters,p,tau_time,varphi,tau_space,station_by_time=TRUE,
                              distribution,sigma_2,add_space=FALSE,add_time=TRUE){

  # Choose how we order spacetime date (station by time or time by station)
  if(station_by_time==TRUE){
    time_grid=expand.grid(station=1:K,time=1:N)
    time_grid=time_grid[,c(2,1)]

  } else{
    time_grid=expand.grid(time=1:N,station=1:K)
  }

  # Generate clusters indexes for all time point in each station with random discrete uniform sampling.
  c_train<-rep(0,K*N)
  c_test<-rep(0,K*N)
  for (k in 1:K) {
    c_train[which(time_grid$station==k)]=sample(1:nb_clusters,N,replace = TRUE)
    c_test[which(time_grid$station==k)]=sample(1:nb_clusters,N,replace = TRUE)
  }

  # MixSim to get X covariates parameters with the same proportion, then generate those X covariate according to those parameters and the c partition generated above
  Melange=MixSim(MaxOmega =max_overlap,BarOmega = bar_overlap,K=nb_clusters,p=p,int=lim_x,resN=1000)
  Data_train=Data_test=matrix(data=NA,nrow = nrow(time_grid)  ,ncol = p)
  for (i in 1:nrow(time_grid)) {
    Data_train[i,]=rmvnorm(1,mean = Melange$Mu[c_train[i],] ,sigma = Melange$S[,,c_train[i]])
    Data_test[i,]=rmvnorm(1,mean = Melange$Mu[c_test[i],] ,sigma = Melange$S[,,c_test[i]])
  }

  # Fuse training and testing data together for future use
  Data=rbind(Data_train,Data_test)
  c=c(c_train,c_test)
  grille=rbind(time_grid,time_grid)
  x_train=Data_train
  x_test=Data_test

  #If only 2 covariates, plot the data clusters
  if(p==2){
    plot(x_train,col=c_train,pch=16,main="covariables de l'annee (train)")
    plot(x_test,col=c_test,pch=16,main="covariables de l'annee (test)")
  }

  # Print Similarity ratios in simulated clusters
  print(overlap_index(scale(x_train),c_train))
  print(overlap_index(scale(x_test),c_test))


  # Generate random space effects if add_space == TRUE
  phi=rep(0,K)
  if(add_space==TRUE){

    # Generate coordinates on an uniform square
    coord=scale(cbind(runif(K,50,54.5),runif(K,11.3,15)))
    plot(coord,main="stations coordinates",pch=16)
    text(coord[,1],coord[,2],1:K)

    # Define spatial distance matrix W based on generated coordinates and compute fixed lambda based on specific formula
    W=as.matrix(dist(coord,diag = TRUE,upper = TRUE))
    lambda=sqrt(2*mean((W[W>0])^2))
    H_lambda=corr_exp(W,lambda,jitter=TRUE)
    H_lambda.inv=solve(H_lambda)

    # Random spatial effects vector
    phi=rmvnorm(n=1,mean=rep(0,K),sigma =tau_space*H_lambda)
  } else{
    lambda="NA"
    coord="NA"
    W="NA"
  }
  phi.long <- rep(phi, N)

  # Generate random time effects if add_time==TRUE, using and AR(1) process
  delta=rep(0,N)
  if(add_time==TRUE){

    delta[1]=rnorm(1,0,sqrt(tau_time))
    for (t in 2:N) {
      delta[t]=rnorm(1,varphi*delta[t-1],sqrt(tau_time))
    }
  }
  delta.long=kronecker(delta, rep(1, K))

  # generate, in each cluster, beta parameters from a standard normal distribution and the sigma variance parameters drawn uniformly around sigma_2 (function arg).
  sim_sigma<-NULL
  sim_beta=matrix(data=rnorm((p+1)*nb_clusters),ncol=nb_clusters,nrow=p+1,byrow = FALSE)
  sim_sigma=runif(nb_clusters,sigma_2-0.25*sigma_2,sigma_2+0.25*sigma_2)

  # Generate the response vectors y based on generated covariates X, parameters beta/sigma, random effects and chosen distribution ("N" or "LN")
  y_train<-NULL
  y_test<-NULL
  for (i in 1:(K*N)) {
    if(distribution=="N"){
      y_train[i]=rnorm(1,mean=sim_beta[1,c_train[i]]+x_train[i,]%*%sim_beta[2:(p+1),c_train[i]]+phi[time_grid$station[i]]+
                         delta[time_grid$time[i] ],sd=sqrt(sim_sigma[c_train[i]]))

      y_test[i]=rnorm(1,mean=sim_beta[1,c_test[i]]+x_test[i,]%*%sim_beta[2:(p+1),c_test[i]]+phi[time_grid$station[i]]+
                        delta[time_grid$time[i]],sd=sqrt(sim_sigma[c_test[i]]))
    }

    if(distribution=="LN"){
      mu_train=sim_beta[1,c_train[i]]+x_train[i,]%*%sim_beta[2:(p+1),c_train[i]]+phi[time_grid$station[i]]+delta[time_grid$time[i]]
      var_train=sim_sigma[c_train[i]]
      y_train[i]=rlnorm(1,meanlog=log(mu_train^2/sqrt(mu_train^2+var_train)) ,sdlog=sqrt(log(1+var_train/(mu_train^2))))

      mu_test=sim_beta[1,c_test[i]]+x_test[i,]%*%sim_beta[2:(p+1),c_test[i]]+phi[time_grid$station[i]]+delta[time_grid$time[i]]
      var_test=sim_sigma[c_test[i]]
      y_test[i]=rlnorm(1,meanlog=log(mu_test^2/sqrt(mu_test^2+var_test)) ,sdlog=sqrt(log(1+var_test/(mu_test^2))))

    }
  }

  # Keep in memory the parameters that generated the data
  sim_param=vector("list", 5 )
  sim_param[[1]]=sim_beta
  sim_param[[2]]=sim_sigma
  sim_param[[3]]=delta
  sim_param[[4]]=phi
  sim_param[[5]]=lambda

  # Keep all the results in a list
  results=vector("list", 5 )
  results[[1]]=data.frame(c_train,time_grid,x_train ,y_train)
  results[[2]]=data.frame(c_test,time_grid,x_test ,y_test)
  results[[3]]=W
  results[[4]]=coord
  results[[5]]=sim_param

  names(results)=c("train_data","test_data","W_matrix","coord","sim_param")

  return(results)
}



#' Simulate clustered data based on an existing dataset.
#'@param clean_data Clean dataset from which you want to generate clustered simulated datasets. Needs to have the following columns: "date", "station",
#' "time", "obs" and "t2m_mean".
#'Rest of the columns are other covariates.
#'@param idx_train vector of training data rows indexes.
#'@param idx_test vector of testing data rows indexes.
#'@param idx_covariate columns indexing which covariates in clean_data to be sampled for simulation, other than "t2m_mean", "lat" and "lon" which are
#' automatically included.
#'@param station_info data.frame containing stations id's, latitudes and longitudes as the following columns respectively: "stations_id", "lat" and "lon".
#'@param nb_clusters number of clusters in simulated dataset.
#'@param p number of total covariates used. Since we use "t2m_mean", "lat", "lon", this determines hot many p-3 covariates are sampled based on correlation with obs.
#'@param K number of stations wanted in the simulated dataset (sampled uniformly from all stations).
#'@param idx_time vector of selected time wanted in the simulated dataset.
#'@param distribution simulated response distribution ("N" for normal distribution or "LN" for log-normal)
#'@param same_stations if TRUE, keep times with a lot of stations and then keep stations present accross all selected times (used to make sure each station has
#' the same time points)
#'@param add_time if TRUE, add time effects (randomly distributed around a sin function).
#'@param time_sorted if TRUE, rearrange the N time points to a 1:N sequence, otherwise keep the original indexing.
#'@param mapping.daymonth mapping matrix to generate the random time effects.
#'@param disp_beta dispersion parameter for the random walk generation of the beta parameters.
#'@param nb_rep number of steps in the random walk.
#'@param sigma_time variance parameters for the simulated time effects.
#'@param factor_sigma variance parameter for the randomly generated betas parameters.
#'@param plot_ALL if TRUE, plot generated random time effects and noises.
#'@param jitter_coord if TRUE, add a jitter to the coordinates.
#'@param jitter_cov if TRUE, add a jitter to the simulated covariates covariance matrix.
#'@param plot.coord if TRUE, plot the selected stations coordinates.
#'@return a list containing the simulated dataset, the scaled coordinates, the clusters indicators vectors and the generated time effects vector.
Simulation_real_data = function(clean_data,idx_train,idx_test,station_info, nb_clusters, p, K, idx_time, distribution,same_stations=TRUE,
                                     add_time = TRUE,time_sorted=FALSE, mapping.daymonth,disp_beta,nb_rep,sigma_time,factor_sigma,
                                     plot_ALL=FALSE,jitter_coord=FALSE,jitter_cov=0.001,plot.coord=FALSE){

  # Draw K stations from all available stations
  stations_tirage=sample(unique(clean_data$station),K)

  # From all clean data, kep only previously drawn stations and time points selected with idx_time
  data_sim <- clean_data[clean_data$station %in% stations_tirage & clean_data$time %in% idx_time ,]

  #Keep only time points common to both training and testing data
  cond=identical( sort(unique(data_sim[idx_train,]$time)),sort(unique(data_sim[idx_test,]$time))
  )
  print(c("Same days in train and test? ",cond))
  if(cond==FALSE){
    idx_fixed_test=sort(unique(data_sim[idx_test,]$time))[sort(unique(data_sim[idx_test,]$time)) %in%sort(unique(data_sim[idx_train,]$time))]
    idx_fixed_train=sort(unique(data_sim[idx_train,]$time))[sort(unique(data_sim[idx_train,]$time)) %in%idx_fixed_test]
    print(c("same time?",identical( idx_fixed_test,idx_fixed_train)))
    data_sim<-data_sim[data_sim$time %in% idx_fixed_test ,]
  }

  #Now, keep stations common to both training and testing data (after keeping common times points)
  cond_2=identical(sort(unique(data_sim[idx_train,]$station)),sort(unique(data_sim[idx_test,]$station)))
  print(c("Same stations at all time points?",cond_2))
  if(cond_2==FALSE){
    idx_fixed_test=sort(unique(data_sim[idx_test,]$station))[sort(unique(data_sim[idx_test,]$station)) %in%sort(unique(data_sim[idx_train,]$station))]
    idx_fixed_train=sort(unique(data_sim[idx_train,]$station))[sort(unique(data_sim[idx_train,]$station)) %in%idx_fixed_test]
    print(c("same station?",identical(idx_fixed_test,idx_fixed_train)))
    data_sim<-data_sim[data_sim$station %in% idx_fixed_test ,]
  }

  #regroup the data
  data_sim<-rbind(data_sim[idx_train,],data_sim[idx_test,])


  #Change stations indexes to 1:K and time indexes to 1:N if time_sorted==TRUE
  N=length(unique(data_sim$time))
  vec_station=sort(unique(data_sim$station))
  for (i in 1:length(vec_station)) {
    data_sim$station[which(data_sim$station==vec_station[i])]=i
  }
  if(time_sorted==TRUE){
    vec_time=sort(unique(data_sim$time))
    for (i in 1:length(vec_time)) {
      data_sim$time[which(data_sim$time==vec_time[i])]=i
    }
  }

  #Coordinates and stations/time points numbers
  coord=station_info[station_info$station %in% stations_tirage,c("lat","lon")]
  if(plot.coord==TRUE){
    plot(coord)
  }
  nb_station=length(unique(data_sim$station))
  nb_time=length(unique(data_sim$time))


  # Uniformly sample covariates with weigths proportionnal to covariates corelation with observed value
  y=data_sim$obs
  corr_x=cor(data_sim$obs,data_sim[,idx_covariate])
  poids=abs(corr_x)/sum(abs(corr_x))
  variables_sim = colnames(data_sim)[sample(idx_covariate,p-3,prob=poids)]

  #Keep selected covariates
  data_sim_select=data_sim[,c("date","time","station","t2m_mean",variables_sim)]

  #If same_stations==TRUE, keep time points with a lot of stations available and keep data with stations common to all time points
  if(same_stations==TRUE){
    #data_sim_select_train=data_sim_select[stri_detect_fixed(data_sim_select$date,"2015")==TRUE,]
    data_sim_select_train=data_sim_select[idx_train,]
    stations_uniques_train=vector("list",length(unique(data_sim_select_train$time)))
    id_boucle=1
    for (i in sort(unique(data_sim_select_train$time))) {
      stations_uniques_train[[id_boucle]]=sort(unique(data_sim_select_train$station[which(data_sim_select_train$time==i)]))
      id_boucle=id_boucle+1
    }
    times_uniques_grand_train=sort(unique(data_sim_select_train$time))[which(lapply(stations_uniques_train,length)==K)]
    stations_uniques_grand_train=stations_uniques_train[which(lapply(stations_uniques_train,length)==K)]

    data_sim_select_test=data_sim_select[idx_test,]
    stations_uniques_test=vector("list",length(unique(data_sim_select_test$time)))
    id_boucle=1
    for (i in sort(unique(data_sim_select_test$time))) {
      stations_uniques_test[[id_boucle]]=sort(unique(data_sim_select_test$station[which(data_sim_select_test$time==i)]))
      id_boucle=id_boucle+1
    }
    times_uniques_grand_test=sort(unique(data_sim_select_test$time))[which(lapply(stations_uniques_test,length)==K)]
    stations_uniques_grand_test=stations_uniques_test[which(lapply(stations_uniques_test,length)==K)]
    times_uniques_grand=intersect(times_uniques_grand_test,times_uniques_grand_train)
    stations_communes=Reduce(intersect,c(stations_uniques_grand_train,stations_uniques_grand_test))
    data_sim_select=data_sim_select[which(data_sim_select$time%in%times_uniques_grand ), ]
    data_sim_select=data_sim_select[which(data_sim_select$station%in%stations_communes), ]
    y=y[which(data_sim_select$time%in%times_uniques_grand )]
    y=y[which(data_sim_select$station%in%stations_communes)]
  }


  #Add coordinates as covariates and jitter to those if jitter_coord==TRUE
  coord_long=kronecker_cpp(as.matrix(rep(1,2*length(unique(data_sim_select$time)))),as.matrix(coord))
  if(jitter_coord==TRUE){
    coord_long=coord_long+rnorm(length(coord_long),0,0.1)
  }
  data_sim_select=cbind(data_sim_select,coord_long)
  for (i in 1:length(unique(data_sim_select$station))) {
    data_sim_select$station[which(data_sim_select$station==sort(unique(data_sim_select$station))[i])]=i
  }


  # Extract days and months from dates
  dates_sim_day = day(data_sim_select$date)
  dates_sim_month = month( data_sim_select$date)
  data_sim_cl = data_sim_select
  data_sim_cl$day = dates_sim_day
  data_sim_cl$month = dates_sim_month

  #Reorder columns
  data_sim_cl=data_sim_cl[,c("date","time","station","day","month","1","2","t2m_mean",variables_sim)]

  #Nb de stations/time points
  K=length(unique(data_sim_cl$station))
  nb_time=length(unique(data_sim_cl$time))
  print(paste("nb stations: ",K,"nb time: ",nb_time,sep=""))

  #Standardized and non-standardized data
  data_sim_std = apply(data_sim_cl[, -c(1,2,3,4,5)], 2, FUN= function(x) { ( x - mean(x))/ sd(x) } )
  data_sim_nstd = data_sim_cl[, -c(1,2,3,4,5)]


  #For each station, get sub-partitions by divided the time points into floor(sqrt(T)) clusters
  ntime.cl = floor(sqrt( nb_time ))
  cl.station = vector("list", K )
  mean.station = vector("list", K )

  for ( i in sort(unique(data_sim_select$station))){
    station_i=which(sort(unique(data_sim_select$station))==i)
    obs.station = which( data_sim_select$station == i )
    data.cl = data_sim_std[ obs.station, ]
    Kmeans_stat = kmeans( data.cl, ntime.cl )
    cl.station[[station_i]] = as.vector(Kmeans_stat$cluster)
    mean.station[[station_i]] = matrix(0, nrow= length( unique(cl.station[[station_i]]) ), ncol= p)
    for ( j in seq(1, length( unique(cl.station[[station_i]])))){
      mean.station[[station_i]][j,] = colMeans( as.matrix( data.cl[which( Kmeans_stat$cluster == j), seq(1, p)]))
    }
  }

  ## Time effects based on sin function
  if ( add_time){
    for ( i2 in seq(1, nrow(data_sim_cl))){
      date_day = data_sim_cl$day[i2]
      date_month = data_sim_cl$month[i2]
      if ( data_sim_cl$month[i2] == 2 & data_sim_cl$day[i2] == 29 ) date_day = 28
      utime = 2*sin( pi + 0.5+ (2*pi/12)*(12/365)*mapping.daymonth[date_month, date_day])
      data_sim_cl$t[i2] = rnorm(1, utime, sd = sqrt(sigma_time))
    }
  }

  #K-means with nb_clusters number of clusters wanted
  mean.station.data = NULL
  for ( i in seq(1, K)){
    mean.station.data = rbind( mean.station.data, mean.station[[i]])
  }

  if ( nb_clusters > nrow(mean.station.data) ) {
    cat("Number of mean data ", nrow(mean.station.data), " is smaller than cluster saught ", nb_clusters, ".\n")
    nb_clusters = floor( sqrt(nrow(mean.station.data) ) )
    cat("nb_clusters has been changed to ", nb_clusters, ".\n")
  }

  Kmeans_sim=kmeans(mean.station.data, nb_clusters) ##each row is a mean
  cl.stat =as.vector(Kmeans_sim$cluster)

  # Assign each observation to her right cluster
  cl = rep(0, nrow(data_sim_select))
  i1 = 1
  for ( i in seq(1, K)){
    obs.station = which( data_sim_select$station == i )
    lk = nrow(mean.station[[i]])
    for ( i2 in seq(1, lk)){
      obs.mean = which( cl.station[[i]] == i2 )
      cl.s = cl.stat[i1]
      cl[obs.station][obs.mean] = cl.s
      i1 = i1 +1
    }
  }
  # If a cluster has too low observations, fuse it with the closest cluster
  cl.sizes = as.vector(table( cl ))
  cl.ssq = Kmeans_sim$withinss
  for ( i in seq(1, nb_clusters)){
    if( i!=1){
      dist_min = sum(( Kmeans_sim$centers[i,] - Kmeans_sim$center[1,])^2)
      ij=1
    } else{
      dist_min = sum(( Kmeans_sim$centers[i,] - Kmeans_sim$center[2,])^2)
      ij=2
    }

    if ( cl.sizes[i] < p+1)
    {
      print(c("cluster too small:",i,"de taille",cl.sizes[i]))
      ##look for bigger cluster
      for ( j in seq(1, nb_clusters)){
        if ( j != i ){
          dist = sum( (Kmeans_sim$centers[j,] - Kmeans_sim$center[i,])^2 )
          if ( dist < dist_min )
          {
            dist_min = dist
            ij = j
          }
        }
      }

      obs.cl = which(cl == i)
      cl[ obs.cl ] = ij
      print(c("deplace vers le cluster:",ij))
      cl.sizes[j] = cl.sizes[j] + cl.sizes[i]
      cl.ssq[j] = cl.ssq[j] + cl.ssq[i]
      cl.sizes[i] = 0
    }
  }
  cl=as.numeric(replace_cpp(cl)+1)
  nb_clusters=length(unique(cl))
  cl.sizes=cl.sizes[which(cl.sizes>0)]


  print(c("R coefficient on all standardized data:",sqrt(summary(lm(y~as.matrix(data_sim_std)))$adj.r.squared)))
  print(c("R coefficient on all non-standardized data:",sqrt(summary(lm(y~as.matrix(data_sim_nstd)))$adj.r.squared)))

  #Find beta parameters in each cluster with linear regression and mean/covariance matrix in each cluster
  beta_cluster=matrix(data=NA,nrow=p+1,ncol = nb_clusters)
  sigma_cluster=rep(0,nb_clusters)
  mean_param=matrix(data=NA,nrow=p,ncol=nb_clusters)
  cov_param=array(data = NA,dim = c(p,p,nb_clusters))

  for (i in seq(1,nb_clusters) ) {
    if ( cl.sizes[i] > 0 ){
      if(distribution=="N"){
        reg=lm(y[which(cl==i)] ~ as.matrix(data_sim_std[which(cl==i), seq(1,p)]) )
        coeff_reg=reg$coefficients
        if(anyNA(reg$coefficients)==TRUE){
          coeff_reg[is.na(coeff_reg)]=0
          print(c("changed coeff to 0 in the cluster ",i))
        }
      }

      if(distribution=="LN"){
        reg=lm(log(y[which(cl==i)]+abs(min(y[which(cl==i)]))+1) ~as.matrix(data_sim_std[which(cl==i), seq(1,p)]) )
        coeff_reg=reg$coefficients
        if(anyNA(reg$coefficients)==TRUE){
          coeff_reg[is.na(coeff_reg)]=0
          print(c("changed coeff to 0 in the cluster ",i))
        }
      }

      #Random walk step to generate the beta parameter
      for (k in 1:length(coeff_reg)) {
        beta_cluster[k,i]= rnorm(1, coeff_reg[k],  min(1, disp_beta*abs(coeff_reg[k])))
        if(nb_rep>0){
          for ( j in seq(1, nb_rep)){
            beta_cluster[k,i]= rnorm(1, beta_cluster[k,i],  min(1, disp_beta*abs(coeff_reg[k])))
          }
        }
      }
      sigma_cluster[i]= runif( 1, factor_sigma*summary(reg)$sigma, 2*factor_sigma*(summary(reg)$sigma) )
      mean_param[,i]=colMeans(as.matrix(data_sim_std[which(cl==i),seq(1,p)] ) )
      cov_param[,,i]= cov(as.matrix(data_sim_std[which(cl==i), seq(1,p)] ) )

    }
  }

  # Generate data based of previously generated clusters parameters by using distribution conditionned on the coordinates to keep them fixed
  X_sim=matrix(data=NA,nrow=length(cl) ,ncol=p)
  y_sim<-NULL
  for (i in 1:length(cl)) {
    X_sim[i,c(1,2)]=as.matrix(data_sim_std[i,c(1,2)])

    Sigma_11=cov_param[-c(1,2),-c(1,2),cl[i]]
    Sigma_12=cov_param[-c(1,2),c(1,2),cl[i]]
    Sigma_21=t(Sigma_12)
    if(det(cov_param[c(1,2),c(1,2),cl[i]])==0 ){
      cov_param[c(1,2),c(1,2),cl[i]]=cov(as.matrix(data_sim_std[which(cl==cl[i]), seq(1,2)])+rnorm(length(data_sim_std[which(cl==cl[i]), seq(1,2)]),0,jitter_cov) )
    }
    Sigma_22.inv=solve_cpp(cov_param[c(1,2),c(1,2),cl[i]],FALSE)
    mu_1=mean_param[-c(1,2),cl[i]]
    mu_2=mean_param[c(1,2),cl[i]]

    mu_cond_cluster=mu_1+Sigma_12%*%Sigma_22.inv%*%(data_sim_std[i,c(1,2)]-mu_2)
    Sigma_cond_cluster=Sigma_11-Sigma_12%*%Sigma_22.inv%*%Sigma_21
    X_sim[i,seq(3,p)] = rmvnorm(1, mean=mu_cond_cluster, sigma = Sigma_cond_cluster)

    if ( add_time ){
      if(distribution=="N"){
        y_sim[i]=rnorm(1, mean=beta_cluster[1,cl[i]] + X_sim[i,]%*%beta_cluster[2:(p+1),cl[i]]+data_sim_cl$t[i] ,sd=sqrt(sigma_cluster[cl[i]]))
      }
      if(distribution=="LN"){
        mu_log=beta_cluster[1,cl[i]] + X_sim[i,]%*%beta_cluster[2:(p+1),cl[i]]+data_sim_cl$t[i]
        var_log=sigma_cluster[cl[i]]
        y_sim[i]=rlnorm(1,meanlog=log(mu_log^2/sqrt(mu_log^2+var_log)) ,sdlog=sqrt(log(1+var_log/(mu_log^2))))
      }
    } else{
      if(distribution=="N"){
        y_sim[i]=rnorm(1, mean=beta_cluster[1,cl[i]] + X_sim[i,]%*%beta_cluster[2:(p+1),cl[i]] ,sd=sqrt(sigma_cluster[cl[i]]))
      }
      if(distribution=="LN"){
        mu_log=beta_cluster[1,cl[i]] + X_sim[i,]%*%beta_cluster[2:(p+1),cl[i]]
        var_log=sigma_cluster[cl[i]]
        y_sim[i]=rlnorm(1,meanlog=log(mu_log^2/sqrt(mu_log^2+var_log)) ,sdlog=sqrt(log(1+var_log/(mu_log^2))))
      }
    }
  }

  for (i in sort(unique(cl))) {
    print(c("Adjusted R squared in cluster:",i,sqrt(max(0,summary(lm(y_sim[which(cl==i)] ~X_sim[which(cl==i),] ))$adj.r.squared))))
  }

  print(c("adjusted R squared for all data: ",sqrt(max(0,summary(lm(y_sim~X_sim))$adj.r.squared))))

  effets_exemple=rnorm(length(cl),mean=0,sd=sqrt(mean(sigma_cluster)))
  print(c("sigma_clusters",sigma_cluster))
  colnames(X_sim)=c("lat","lon","t2m_mean",variables_sim)
  t_mean=rep(0,length(unique(data_sim_cl[idx_train,]$time)))
  effets_exemple_mean=rep(0,length(unique(data_sim_cl[idx_train,]$time)))

  for (i in sort(unique(data_sim_cl[idx_train,]$time))){
    t_mean[which(sort(unique(data_sim_cl[idx_train,]$time))==i)]=mean( data_sim_cl[idx_train,]$t[which(data_sim_cl[idx_train,]$time==i)])
    effets_exemple_mean[which(sort(unique(data_sim_cl[idx_train,]$time))==i)]=mean(effets_exemple[which(data_sim_cl[idx_train,]$time==i)])
  }

  if(plot_ALL==TRUE){
    par(mfrow=c(2,2))
    plot(sort(unique(data_sim_cl[idx_train,]$time)),t_mean,main="space-average time effects at each time point for the training data")
    plot(sort(unique(data_sim_cl[idx_train,]$time)),effets_exemple_mean,main="space-average independent effects at each time point for the training data")
    hist(effets_exemple,main="independent effects")
    hist(data_sim_cl$t,main="time effects")
  }



  #Keep the results
  results=vector("list", 4)
  results[[1]]=data.frame(date=data_sim_select[,1] ,time=data_sim_select[,2],station=data_sim_select[,3],X_sim ,y_sim)
  results[[2]]=scale(coord)
  results[[3]]=cl
  results[[4]]=data_sim_cl$t

  names(results)=c("data", "scaled_coordinates", "partition", "random_effects")

  return(results)
}


