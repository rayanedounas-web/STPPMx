

#' Compute relevant matrices on a parameter grid
#'
#' Pre-compute two relevant matrices to be then used in the STPPMx model with mixed effects. Like explained in my paper, this allows to greatly speed up
#' mixed effects STPPMx model by pre-calculating fastidious matrix for a fixed grid of selected parameters.
#'
#'@param grid_varphi vector of time correlation parameter values such that we calculates the matrices for each.
#'@param lambda space correlation parameter.
#'@param W Distance matrix between stations coordinates.
#'@param grid data.frame or matrix with rows = observations, first column = time points and second columns = stations.
#'@param jitter_H if TRUE, add a jitter to the diagonal of the space correlation matrix.
#'@param jitter_varphi if TRUE, add a jitter to the diagonal of the time correlation matrix.
#'@returns List containing, for each varphi, matrices with 2 columns each containing matrices calculated for each time point.
grid_matrices <-function(grid_varphi,lambda,W,grid,jitter_H=FALSE,jitter_varphi=TRUE){
  list_matrices=vector("list",length(grid_varphi))

  # Lambda space correlation matrix and its inverse calculated based on distance matrix
  H_lambda=corr_exp(W,lambda,jitter=jitter_H)
  H_lambda.inv=solve(H_lambda)


  #Compute the matrices for each varphi value in grid_varphi
  for (k in 1:length(grid_varphi)) {
    print(c("step: ", k))

    #Varphi time correlation matrix and its inverses for the selected varphi value.
    varphi=grid_varphi[k]
    distance_time=as.matrix(dist(sort(unique(grid_train[,1])),diag=TRUE,upper=TRUE))
    H_varphi=corr_exp(distance_time,varphi,jitter=jitter_varphi)
    H_varphi.inv=solve_cpp(H_varphi,FALSE)

    # Initiate a matrix of lists that will store all computed matrices for each time points
    matrix_varphi=matrix(data=list(),nrow =length(unique(grid[,1])) ,ncol =2 )

    for (i in sort(unique(grid[,1]))) {
      time_i=which(unique(grid[,1])==i)
      mat1=kronecker_cpp(t(as.matrix(H_varphi[time_i,-time_i])),H_lambda)
      mat2=kronecker_cpp(solve_cpp(as.matrix(H_varphi[-time_i,-time_i]),FALSE),solve_cpp(H_lambda,FALSE))

      matrix_varphi[[time_i,1]]=mat1%*%mat2
      matrix_varphi[[time_i,2]]=solve_cpp(H_lambda-mat1%*%mat2%*%t(mat1),FALSE)
    }

    list_matrices[[k]]=matrix_varphi
  }

  return(list_matrices)
}

#' STPPMX model
#'
#' Run the STPPMX model MCMC training described in chapter 4 of my paper (citation in package description) and returns the MCMC parameters samples to be used with the
#' STPPMx_predict function to get useful metric from the predictive distribution.
#'
#'@param W distance matrix between stations coordinates.
#'@param y vector of training response values.
#'@param X_sim training covariates matrix used for similarity (rows = observations, columns = covariates).
#'@param X_reg training covariates matrix used for regression (rows = observations, columns = covariates, first column = unit vector).
#'@param grid data.frame or matrix with rows = observations, first column = time points and second columns = stations.
#'@param nb_clusters_start number of starting clusters.
#'@param alpha fixed alpha parameter for the Dirichlet process (Could be used as a starting parameter in the future).
#'@param v shape parameter of the inverse-gamma prior.
#'@param s rate parameter of the inverse-gamma prior.
#'@param m_0 mean parameter for the similarity Gaussian prior.
#'@param Phi variance matrix parameter for the similarity inverse-Wishart prior.
#'@param w degrees of freedom for the similarity inverse-Wishart prior.
#'@param M number of auxiliary parameters to draw in Neal algorithm.
#'@param B number of MCMC iterations.
#'@param R number of burn-in iterations.
#'@param tau_time variance parameter for the inverse-gamma prior of the random time effects.
#'@param tau_station variance parameter for the inverse-gamma prior of the random space effects.
#'@param tau variance parameter for the inverse-gamma prior of the random spatio-temporal effects.
#'@param lambda space correlation parameter.
#'@param varphi time correlation parameter.
#'@param MC_param vector of TRUE/FALSE values indicating whether to use fixed value or MCMC draw for the
#'following parameters (in order) : varphi, lambda, tau_time, tau_station, alpha.
#'@param MH_vec vector containing values of the instrumental distribution step used for each of the above parameters (same order).
#'@param jitter_lambda if TRUE, add a jitter to the diagonal of the space correlation matrix.
#'@param jitter_varphi if TRUE, add a jitter to the diagonal of the time correlation matrix.
#'@param add_time if TRUE, use the random time effects model.
#'@param add_station if TRUE, use the random space effects model.
#'@param add_mixed if TRUE, use the random spatio-temporal effects model.
#'@param precalc if TRUE, use a grid of varphi values instead of a fixed value and pre-compute the associate matrices for speed.
#'@param fast if TRUE, use a faster Neal's algorithm partition sampling implementation that uses more memory to prevent redundant computations.
#'@param print.state if TRUE, print important information on the MCMC draws at each iteration.
#'@return List containing MCMC draws and acceptance probabilities.
STPPMx <- function(W,y,X_sim,X_reg,grid,nb_clusters_start,alpha,v,s,m_0,Phi,w,M,B,R,tau_time,tau_station,tau,lambda,varphi,MC_param,MH_vec,
                   jitter_lambda=TRUE,jitter_varphi = TRUE,add_time=FALSE,add_station=FALSE, add_mixed=FALSE,precalc=TRUE,fast = FALSE,print.state=TRUE){

  ################## Important values #################################################
  # Regression/similarity covariates dimensions, number of time and space points, long format time and space effects vectors.
  p_sim=dim(X_sim)[2]
  p_reg=dim(X_reg)[2]-1
  nb_time=length(unique(grid[,1]))
  nb_station=length(unique(grid[,2]))
  n=dim(X_reg)[1]
  time_long=rep(0,n)
  station_long=rep(0,n)

  ################################ Parameters initialization ###########################################
  # kmeans with nb_clusters_start clusters to initialize the partition c
  Kmeans=kmeans(X_sim,nb_clusters_start)
  c=Kmeans$cluster
  k=max(c)

  # Initializing the observation vector y density parameters beta and sigma in each clusters with linear regression.
  # If number of observation in a cluster is smaller than number of covariates, draw beta randomly and take sigma_2 as sample variance.
  beta=matrix(rep(0,k*(p_reg+1)),nrow=p_reg+1,ncol=k)
  sigma=c(rep(0,k))
  for (j in 1:k) {
    if((p_reg+1)<length(which(c==j))){
      regression=lm(y[which(c==j)]~X_reg[which(c==j),-1])
      beta[,j]=replace(regression$coefficients,which(is.na(regression$coefficients)),0)
      sigma[j]=(summary(regression)$sigma)^2
    }
    else{
      print(c("Careful!, less observation than covariates in cluster h = ",j))
      sigma[j]=var(y[which(c==j)])
      beta[,j]=rmvnorm(1,rep(0,p_reg+1),diag(p_reg+1))
    }
  }


  # Initializing random effects vectors as 0's vectors and compute associated correlation matrices.
  f_station=rep(0,nb_station);f_time=rep(0,nb_time);f_mixed=rep(0,nb_station*nb_time)
  distance_time=as.matrix(pmin(365-dist(sort(unique(grid[,1])),upper = TRUE,diag=TRUE),dist(sort(unique(grid[,1])),upper = TRUE,diag=TRUE)))

  # Spatial effects model
  if(add_station==TRUE){
    H_lambda=corr_exp(W,lambda,jitter=jitter_lambda)
    H_lambda.inv=solve_cpp(H_lambda,FALSE)
  }

  # Temporal effects model
  if(add_time==TRUE){
    H_varphi=corr_exp(distance_time,varphi,jitter=jitter_varphi)
    H_varphi.inv=solve_cpp(H_varphi,FALSE)
  }

  # Spatio-temporal effects model
  if(add_mixed==TRUE){
    H_lambda=corr_exp(W,lambda,jitter=jitter_lambda)
    H_lambda.inv=solve_cpp(H_lambda,FALSE)

    H_varphi=corr_exp(distance_time,varphi,jitter=jitter_varphi)
    H_varphi.inv=solve_cpp(H_varphi,FALSE)

    H_mixed=kronecker_cpp(H_varphi,H_lambda)
    H_mixed.inv=kronecker_cpp(H_varphi.inv,H_lambda.inv)

    f_mixed=rnorm(nb_station*nb_time)

    # If precalc is TRUE, compute fastidious Kronecker product for a grid of fixed varphi time effects parameter and keep them in memory.
    if(precalc==TRUE){
      list_matrices=vector("list",length(varphi))

      # Lambda space correlation matrix and its inverse calculated based on distance matrix
      H_lambda=corr_exp(W,lambda,jitter=jitter_lambda)
      H_lambda.inv=solve(H_lambda)

      #Compute the matrices for each varphi value in grid_varphi
      for (k in 1:length(varphi)) {
        #Varphi time correlation matrix and its inverses for the selected varphi value.
        print(c("varphi step: ", k, "on" , length(varphi)))
        varphi_k=varphi[k]
        distance_time=as.matrix(pmin(365-dist(sort(unique(grid[,1])),diag=TRUE,upper=TRUE ) ,dist(sort(unique(grid[,1])),diag=TRUE,upper=TRUE )))
        H_varphi_k=corr_exp(distance_time,varphi_k,jitter=jitter_varphi)
        H_varphi_k.inv=solve_cpp(H_varphi_k,FALSE)

        # Initiate a matrix of lists that will store all computed matrices for each time points
        matrix_varphi=matrix(data=list(),nrow =length(unique(grid[,1])) ,ncol =2 )

        # Compute mat1 = Sigma_12 and mat2 = Sigma_22 inverse used to compute mu_bar and Sig_bar posterior parameters.
        for (i in sort(unique(grid[,1]))) {
          time_i=which(unique(grid[,1])==i)
          mat1=kronecker_cpp(t(as.matrix(H_varphi_k[time_i,-time_i])),H_lambda)
          mat2=kronecker_cpp(solve_cpp(as.matrix(H_varphi_k[-time_i,-time_i]),FALSE),solve_cpp(H_lambda,FALSE))
          matrix_varphi[[time_i,1]]=mat1%*%mat2
          matrix_varphi[[time_i,2]]=solve_cpp(H_lambda-mat1%*%mat2%*%t(mat1),FALSE)
        }
        list_matrices[[k]]=matrix_varphi
      }
    }
  }

  # If fast = TRUE, set initial quantities for weights w used for the fast STPPMx model.
  if (fast == TRUE){
    sizes_mem=as.numeric(table(c))
    sums_mem=matrix(data=NA,nrow =p_sim ,ncol=k)
    products_mem=array(data=NA,dim=c(p_sim,p_sim,k))
    for (j in 1:k) {
      sums_mem[,j]=colSums(X_sim[which(c==j),])
      products_mem[,,j]=t(X_sim[which(c==j),])%*%X_sim[which(c==j),]
    }
  }



  # Matrix and data.frame that will be used to store MCMC parameters after burn-in.
  Parameters = vector("list",13)
  #time_calcul=rep(0,5)
  #proba_accept=matrix(data = NA,nrow=B-R,ncol = 5)


  ########################################################### MCMC algorithm iterations #######################################################
  for (l in 1:B) {

    # Long format random effects vectors
    for (r in 1:n) {
      if(add_time==TRUE){
        time_r=which(sort(unique(grid[,1]))==grid[r,1])
        time_long[r]=f_time[time_r]
      }
      if(add_station==TRUE){
        station_long[r]=f_station[grid[r,2]]
      }
    }

    ########################## Partition updating step using modified Neal's algorithm  #############################
    Start_Neal=Sys.time()

    # Sample partition using Neal's modified algorithm with the fast or regular method and keep updated partition and associated parameters.
    if(fast == TRUE){
      Tirage=partitions_draw_cpp_fast(as.matrix(y),X_sim,X_reg[,-1],M,alpha,c-1,v,s,sigma,
                                                  beta,as.matrix(m_0),Phi,w,as.matrix(time_long),as.matrix(station_long),as.matrix(f_mixed),sizes_mem,sums_mem,products_mem)
      sizes_mem=Tirage$sizes_mem;sums_mem=Tirage$sums_mem;products_mem=Tirage$products_mem
    } else{
      Tirage=partitions_draw_cpp(as.matrix(y),X_sim,X_reg[,-1],M,alpha,c-1,v,s,sigma,
                                             beta,as.matrix(m_0),Phi,w,as.matrix(time_long),as.matrix(station_long),as.matrix(f_mixed))
      }
    c=as.vector((Tirage$partition)+1)
    beta=Tirage$beta;sigma=as.vector(Tirage$sigma)

    # If we are after burn-in, store the partition, the clusters sizes and if p = 2, print the covariates graph colored by clusters.
    if(l>R){
      Parameters[[1]]=c
      taille<-NULL
      for (a in 1:length(unique(c))) {
        taille[a]=sum(c==a)
      }
      Parameters[[2]]=taille
    }
    End_Neal=Sys.time()

    ###################### Parameters and hyperparameters Gibbs/Metropolis sampling #################################
    Start_post=Sys.time()
    # In each cluster, compute the posterior beta and sigma_2 parameters.
    for (j in sort(unique(c))) {
      n_c=length(which(c==j));y_c=t(t(y[which(c==j)]));w_c=y_c;grid_c=grid[which(c==j),]
      X_c=X_reg[which(c==j),];f_c=matrix(time_long[which(c==j)]+station_long[which(c==j)]+f_mixed[which(c==j)] ,ncol = 1)
      if(n_c==1){
        X_c=matrix(X_c,nrow=1)
      }
      sum_X_c=vector("list",n_c)

      for(j_2 in 1:n_c){
        sum_X_c[[j_2]]=t(t(X_c[j_2,]))%*%X_c[j_2,]
      }
      lambda_c=Reduce('+',sum_X_c)+diag(p_reg+1)
      mu_c=(chol2inv(chol(lambda_c)))%*%(t(X_c)%*%(w_c-f_c))
      sigma[j]=rinvgamma(1,shape=(v+n_c)/2,0.5*(s+t(w_c)%*%w_c-2*(t(f_c)%*%w_c)+t(f_c)%*%f_c-(t(mu_c)%*%(lambda_c)%*%mu_c)))

      # If only 1 cluster, use a vector instead of a matrix to store parameters.
      if(n_c==n){
        beta=t(rmvnorm(1,mean = mu_c, sigma = chol2inv(chol(lambda_c))*sigma[j]))
      } else{
        beta[,j]=t(rmvnorm(1,mean = mu_c, sigma = chol2inv(chol(lambda_c))*sigma[j]))
      }
    }

    ############################################# If add_time is TRUE, compute posterior random time effects #################################
    # random time effects Gibbs sampling
    if(add_time==TRUE){
      for (t in sort(unique(grid[,1]))) {
        # Indexes of observations at time t
        index_time=which(grid[,1]==t)
        time_i=which(sort(unique(grid[,1]))==t)
        c_time=c[which(grid[,1]==t)];y_time=betaX_time=sigma_time=rep(0,length(index_time))

        # Posterior time effects Gibbs sampling
        sigma_12=H_varphi[time_i,-time_i];sigma_21=t(t(sigma_12));sigma_22.inv=solve_cpp(H_varphi[-time_i,-time_i],FALSE)
        u_post=sigma_12%*%sigma_22.inv%*%f_time[-time_i];sigma_post=tau_time*(1-sigma_12%*%sigma_22.inv%*%sigma_21)
        for (j in 1:length(index_time)) {
          sigma_time[j]=sigma[c_time[j]];betaX_time[j]=X_reg[index_time[j],]%*%beta[,c_time[j]];y_time[j]=y[index_time[j]]
        }
        s_v=(1/sigma_post+sum(1/sigma_time))^(-1);mu_v=(sum((y_time-betaX_time)/sigma_time)+u_post/sigma_post)*s_v

        # Sample random time effects
        f_time[time_i]=rnorm(1,mean = mu_v,sd=sqrt(s_v))
      }
    }

    # tau_time variance parameter Gibbs sampling
    if(add_time==TRUE){
      if(MH_vec[3]==TRUE){
        a_time=1+nb_time
        b_time=0.01+t(f_time)%*%H_varphi.inv%*%f_time
        tau_time=rinvgamma(1,shape = a_time/2,rate = b_time/2)
        print(c("tau_time",tau_time))
        cond_station=TRUE
      } else{
        cond_station="fixe"
      }
    }

    # Varphi Metropolis-Hasting sampling
    if(add_time==TRUE){
      if(MH_vec[1]==TRUE ){
        c_varphi=MC_param[1]
        varphi_metropolis=runif(1,varphi-c_varphi,varphi+c_varphi)

        # Reject the varphi candidate if its not in the designed region
        if(varphi_metropolis< 0 | varphi_metropolis>0.95){
          proba_varphi=0
          cond_varphi=FALSE
        } else{
          H_varphi_metropolis=corr_exp(distance_time,varphi_metropolis,jitter = jitter_varphi)
          H_varphi_metropolis.inv=solve_cpp(H_varphi_metropolis,FALSE)
          ratio_varphi=-0.5*(log(det(H_varphi_metropolis)-log(det(H_varphi))))-(1/(2*tau_time))*(f_time%*%H_varphi_metropolis.inv%*%f_time-f_time%*%H_varphi.inv%*%f_time )
          proba_varphi=min(0,ratio_varphi)
          cond_varphi=log(runif(1,0,1))<proba_varphi
        }

        if(cond_varphi==TRUE){
          varphi=varphi_metropolis
        }
        print(c("varphi",varphi))
      } else{
        cond_varphi="fixe"
      }
    }

    ############################################# If add_station is TRUE, compute posterior random space effects ###########################################
    # random space effects vector Gibbs sampling
    if(add_station==TRUE){
      for (j in 1:nb_station) {
        index_station=which(grid[,2]==j);c_station=c[which(grid[,2]==j)]
        sigma_station=betaX_station=y_station=time_vec=rep(0,length(index_station))

        for (i in 1:length(index_station)){
          sigma_station[i]=sigma[c_station[i]];betaX_station[i]=X_reg[index_station[i],]%*%beta[,c_station[i]]
          y_station[i]=y[index_station[i]]
          time_vec[i]=f_time[grid[i,1]]
        }
        sigma_12=H_lambda[-j,j]
        sigma_22.inv=solve(H_lambda[-j,-j])
        sigma_21=t(t(sigma_12))

        s_u=(1/((1-sigma_12%*%sigma_22.inv%*%sigma_21)*tau_station)+sum(1/sigma_station))^(-1)
        mu_u=((sigma_12%*%sigma_22.inv%*%f_station[-j])/((1-sigma_12%*%sigma_22.inv%*%sigma_21)*tau_station)+sum((y_station-betaX_station-time_vec)/sigma_station))*s_u
        f_station[j]=rnorm(1,mean = mu_u,sd=sqrt(s_u))
      }
    }

    # tau_station variance parameter Gibbs sampling
    if(add_station==TRUE){
      if(MH_vec[4]==TRUE){
        a_station=1+nb_station
        b_station=0.001+t(f_station)%*%H_lambda.inv%*%f_station
        tau_station=rinvgamma(1,shape = a_station/2,rate = b_station/2)
        print(c("var",ifelse(a_station/2>2,((b_station/2)^2)/((a_station/2-2)*(a_station/2-1)^2),0) ))
        print(c("mode",(b_station/2)/(a_station/2-1)) )
        cond_time=TRUE

      } else{
        cond_time="fixe"
      }
      print(c("tau_station",tau_station))
    }

    # lambda parameter Metropolis-Hastings sampling
    if(add_station==TRUE){
      if(MH_vec[2]==TRUE){
        c_lambda=MC_param[2]
        lambda_metropolis=runif(1,lambda-c_lambda,lambda+c_lambda)
        if(lambda_metropolis<0 | lambda_metropolis>3){
          proba_lambda=0
          cond_lambda=FALSE
        } else{
          H_lambda_metropolis=corr_exp(W,lambda_metropolis,jitter = TRUE)
          H_lambda_metropolis.inv=solve(H_lambda_metropolis)

          ratio_lambda=-0.5*(log(det(H_lambda_metropolis))-log(det(H_lambda)))
          -(1/(2*tau_time))*(t(f_station)%*%(H_lambda_metropolis.inv)%*%f_station-t(f_station)%*%(H_lambda.inv)%*%f_station)

          proba_lambda=min(0,ratio_lambda)
          cond_lambda=log(runif(1,0,1))<proba_lambda
        }
        if(cond_lambda==TRUE){
          lambda=lambda_metropolis
          H_lambda=H_lambda_metropolis
          H_lambda.inv=H_lambda_metropolis.inv
        }

      } else{
        cond_lambda="fixe"
      }
    }

    #############################################  If add_mixed is TRUE, compute posterior random mixed effects  #########################################
    if(add_mixed==TRUE){
      # For each time i, compute posterior spatial effects vector
      for (i in sort(unique(grid[,1])) ){

        # Observations, covariates and beta/sigma parameters at time i
        X_time=X_reg[which(grid[,1]==i),]
        y_time=y[which(grid[,1]==i)]
        sigma_long=rep(0,length(unique(grid[which(grid[,1]==i),2])))
        beta_long=matrix(data=NA,nrow=p_reg+1,ncol=length(unique(grid[which(grid[,1]==i),2])))
        X_beta_long=rep(0,length(unique(grid[which(grid[,1]==i),2])))
        id_test=1
        for (j in unique(grid[which(grid[,1]==i),2]) ){
          sigma_long[id_test]=sigma[c[which(grid[,1]==i & grid[,2]==j)]]
          beta_long[,id_test]=beta[,c[which(grid[,1]==i & grid[,2]==j)]]
          X_beta_long[id_test]=X_reg[which(grid[,1]==i & grid[,2]==j),]%*%beta_long[,id_test]
          id_test=id_test+1
        }

        # If precalc is TRUE, use pracalculated matrices to omputed posterioir random effect conditional distribution parameters
        if(precalc==TRUE){
          mu_bar=list_matrices[[id_varphi]][i,1][[1]]%*%f_mixed[which(grid[,1]!=i)]
          Sig_bar.inv=(1/tau)*list_matrices[[id_varphi]][i,2][[1]]
          Sig_post=solve_cpp(Sig_bar.inv+diag(1/sigma_long),FALSE)
          mu_post=Sig_post%*%(Sig_bar.inv%*%mu_bar+diag(1/sigma_long)%*%(y_time-X_beta_long))
        }
        # Else, compute it all without using the precalcutated matrices.
        else{
          Sigma_11=H_lambda
          Sigma_12=kronecker_cpp(t(as.matrix(H_varphi[i,-i])),H_lambda);Sigma_21=t(Sigma_12)
          Sigma_22_inv=kronecker_cpp(solve_cpp(as.matrix(H_varphi[-i,-i]),FALSE),solve_cpp(H_lambda,FALSE))
          mu_bar=Sigma_12%*%Sigma_22_inv%*%f_mixed[which(grid[,1]!=i)]
          Sig_bar.inv=(1/tau)*solve_cpp(H_lambda-Sigma_12%*%Sigma_22_inv%*%Sigma_21,FALSE)
        }
        # random spatio-temporal effects Gibbs sampling.
        Sig_post=solve_cpp(Sig_bar.inv+diag(1/sigma_long),FALSE)
        mu_post=Sig_post%*%(Sig_bar.inv%*%mu_bar+diag(1/sigma_long)%*%(y_time-X_beta_long))
        f_mixed[which(grid[,1]==i)]=rmvnorm(1,mu_post,Sig_post,method = "svd")
      }
      print(c("f_mixed",mean(f_mixed),sqrt(var(f_mixed))))

      # tau variance parameter Gibbs sampling
      if(MH_vec[6]==TRUE){
        if(tau>0.005){
          a_mixed=1+nb_time*nb_station
          b_mixed=0.01+as.numeric(t(f_mixed)%*%H_mixed.inv%*%f_mixed)
          tau=rinvgamma(1,shape = a_mixed/2,rate = b_mixed/2)
        } else{
          tau=0.005
        }
        print(c("tau",tau))
        if(tau>10){
          tau=9
        }
      }

      # varphi parameter Metropolis Hastings sampling
      if(MH_vec[1]==TRUE ){
        c_varphi=MC_param[1]
        varphi_metropolis=runif(1,varphi-c_varphi,varphi+c_varphi)
        print(c("varphi_metropolis",varphi_metropolis,varphi))

        if(varphi_metropolis< (-0.95) | varphi_metropolis>0.95){
          proba_varphi=0
          cond_varphi=FALSE
        } else{
          H_varphi_metropolis=corr_exp(distance_time,varphi_metropolis,jitter=jitter_varphi)
          H_varphi_metropolis.inv=solve_cpp(H_varphi_metropolis,FALSE)
          H_mixed_metropolis.inv=kronecker_cpp(H_varphi_metropolis.inv,H_lambda.inv)

          ratio_varphi=-(t(f_mixed)%*%H_mixed_metropolis.inv%*%f_mixed-t(f_mixed)%*%H_mixed.inv%*%f_mixed)/(2*tau)
          ratio_varphi=as.numeric(ratio_varphi)
          proba_varphi=min(0,ratio_varphi)
          cond_varphi=log(runif(1,0,1))<proba_varphi
        }
        if(cond_varphi==TRUE){

          if(approx_varphi==TRUE){
            varphi=grid_varphi[which.min(abs(grid_varphi - varphi_metropolis))]
            id_varphi=which(grid_varphi==varphi)
          } else{
            varphi=varphi_metropolis
          }
          H_mixed.inv=H_mixed_metropolis.inv
        }
      } else{
        cond_varphi="fixe"
      }
      print(c("varphi",varphi))
    }
    End_post=Sys.time()

    ################################## Saving parameters after burn-in #################################
    # Keep parameters and print them if print.state is TRUE
    if(print.state==TRUE){
      if(l%%1==0 | l==1){
        print(c("partition c drawing time ",difftime(End_Neal,Start_Neal,units = "secs")))
        print(c("parameters/hyperparameters drawing time",difftime(End_post,Start_post,units = "secs")))
        print(paste("iteration: ",l,sep = ""))
        print(paste("number of clusters: ",length(unique(c)),sep=""))
      }
    # If only two covariates, print covariates graph colored by clusters
      if(p_sim==2){
        plot(X_sim,col=c,pch=19,cex=1)
      }
    }

    if(l>R){
      Parameters[[3]]=beta
      Parameters[[4]]=sigma
      Parameters[[5]]=f_station
      Parameters[[6]]=f_time
      Parameters[[7]]=varphi
      Parameters[[8]]=lambda
      Parameters[[9]]=alpha
      Parameters[[10]]=tau_time
      Parameters[[11]]=tau_station
      Parameters[[12]]=f_mixed
      Parameters[[13]]=tau

    }
  }
  Results_MCMC=as.data.frame(Parameters)
  colnames(Results_MCMC)=c("partition","clusters_sizes","beta","sigma","space_effects","time_effects","varphi","lambda","alpha","tau_time","tau_station","mixed_effects","tau")


  #Results=vector("list",2)
  #Results[[1]]=Results_MCMC
  #Results[[2]]=proba_accept

  return(Results_MCMC)
}


#' STPPMx model predictions and scores.
#'
#' Get an STPPMx model predictions and scores (MAE/CRPS) based on results given by the STPPMx functions and testing data.
#'
#'@param X_train_sim training covariates matrix used for similarity (rows = observations, columns = covariates).
#'@param X_test_sim testing covariates matrix used for similarity (rows = observations, columns = covariates).
#'@param X_test_reg testing covariates matrix used for regression (rows = observations, columns = covariates, first column = unit vector).
#'@param y_test vector of testing response values.
#'@param grid_test data.frame or matrix with rows = testing observations, first column = time points and second columns = stations.
#'@param MCMC_results results list of MCMC parameters samples returned by the STPPMx function.
#'@param w degrees of liberty for the similarity inverse-Wishart prior.
#'@param Phi variance matrix parameter used in training for the similarity inverse-Wishart prior.
#'@param v shape parameter used in training of the inverse-gamma prior.
#'@param s rate parameter used in training of the inverse-gamma prior.
#'@param m_0 mean parameter used in training for the similarity Gaussian prior.
#'@param nb_draws number of draws to approximate the predictive distribution using the method described in my thesis.
#'@return list containing a data.frame with predicted values, CRPS and AE at each spacetime points, a vector with the mean CRPS and MAE and a matrix containing the
#' samples that were drawn to approximate the CRPS.
STPPMx_predict<- function(X_train_sim,X_test_sim,X_test_reg,y_test,grid_test,MCMC_results,w,Phi,v,s,m_0,nb_draws){

  # Extract all te MCMC parameters in the MCMC_results list
  beta=MCMC_results$beta
  sigma=MCMC_results$sigma
  spatial_effects=MCMC_results$space_effects
  time_effects=MCMC_results$time_effects
  mixed_effects=MCMC_results$mixed_effects
  varphi=MCMC_results$varphi
  lambda=MCMC_results$lambda
  clusters_sizes=MCMC_results$clusters_sizes
  partition=MCMC_results$partition
  alpha=unlist(MCMC_results$alpha)

  # Shift the partition by 1 since indexes starts at 0 in C++.
  partition_shifted=lapply(partition, function(x, subt_amnt) x - subt_amnt, subt_amnt = 1)
  p_sim=ncol(X_test_sim)
  p_reg=ncol(X_test_reg)-1

  # Number of MCMC samples and observations to predict
  nb_density=dim(MCMC_results)[1]
  nb_obs=nrow(X_test_reg)

  # Initialize results vectors and matrix
  CRPS_vec=rep(0,nb_obs)
  predi_vec=rep(0,nb_obs)
  matrix_sample=matrix(data = NA, nrow=nb_obs, ncol = nb_draws)

  # For each observation i, we emit a prediction and sample nb_draws values from the predictive distribution to compute the approximate CRPS.
  for (i in 1:nb_obs) {
    vec_time=sort(unique(grid_test[,1]))
    station=grid_test[i,2]
    time=grid_test[i,1]
    time_i=which(vec_time == time)


    predi_i=predict_cpp_new(i-1,station-1,time_i-1,X_train_sim,as.matrix(X_test_sim[i,]),as.matrix(X_test_reg[i,-1]),y_test[i] ,beta,sigma,spatial_effects,time_effects,mixed_effects,
                                     clusters_sizes,partition_shifted,alpha,nb_density,nb_draws,w,Phi,v,s,m_0)
    CRPS_vec[i]=predi_i$CRPS
    predi_vec[i]=predi_i$predictions
    matrix_sample[i,]=predi_i$sample
  }

  Results_predictions=vector("list",3)

  # Dataframe with observations, predictions, absolute error and CRPS values.
  Results_predictions[[1]]=data.frame(grid_test, observation = y_test, prediction = predi_vec, AE = abs(y_test - predi_vec), CRPS = CRPS_vec)

  # MAE and mean CRPS
  Results_predictions[[2]]=c(MAE = MSE(predi_vec,y_test,TRUE), mean_CRPS = mean(CRPS_vec))

  #Samples used to approximate CRPS
  Results_predictions[[3]]=matrix_sample

  names(Results_predictions) = c("Results", "mean_metrics", "predictive_samples")
  return(Results_predictions)
}



