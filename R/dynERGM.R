#########################################################################################################
## Loading the required packages
require(lda)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(ggplot2)
library(abind)
library(reshape)
library(combinat)
require(network)
library(Matrix)
library(abind)

#########################################################################################################
##################################   Dynamic network code (HMM)   #######################################
#########################################################################################################
## Defining a wrapper function for HMM undirected density case to run for different number of clusters 
wrapper_HMM_dyn_undir_Dens<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=1,sim_indicator_TERGM=0,theta_true=NA,pi_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(HMMEMundir)
  library(Matrix)
  
  #################################################
  ## Defining a function to update variational parameters Tau 
  Tau.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_data){
    Tau.next.temp<-Tau_update_HMM_undir(gamma=gamma.curr, log_pi=log(pi.curr), theta=theta.curr, network=network, N=N, K=K, T_data = T_data)
    ## converting in 4-dimensional array and normalizing Tau(i,k,)
    Tau.next<-array(NA_real_,dim=c(K,K,N,T_data-1))
    for (t in 1:(T_data-1)){
      for (i in 1:N){
        for (k in 1:K){
          Tau.next[k,,i,t]<-Tau.next.temp[[N*(t-1)+i]][k,]/sum(Tau.next.temp[[N*(t-1)+i]][k,])
        }
      }
    }
    
    return(Tau.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  gamma_init.update.wrapper<-function(gamma_init, log_alpha, theta, network_init, N, K){
    gamma_init.next<-gamma_init_update_HMM_undir(gamma_init=gamma_init, log_alpha=log_alpha, theta=theta, network_init=network_init, N=N, K=K)
    ## Normalization of gamma_init
    for (i in 1:N){
      gamma_init.next[i,]<-gamma_init.next[i,]/sum(gamma_init.next[i,])
    }
    return(gamma_init.next)
  }
  
  #################################################
  ## Defining a function to update K*K matrix of pi
  pi.update<-function(gamma.curr,Tau.curr,N,K,T_data){
    pi.next<-array(NA_real_,dim=c(K,K,T_data-1))
    for (t in 2:T_data){
      for (k in 1:K){
        for (l in 1:K){
          pi.next[k,l,t-1]<-as.numeric((gamma.curr[,k,t-1]%*%Tau.curr[k,l,,t-1]))
          #denom<-denom+sum(gamma.curr[,k,t-1])
        }
        pi.next[k,,t-1]<-pi.next[k,,t-1]/sum(pi.next[k,,t-1]) ## Normalization
      }
    }
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_HMM_undir(theta=as.vector(theta.curr), gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    hess<-hess_HMM_undir(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,T_data,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    pi<-array(NA_real_,dim=c(K,K,T_data-1,n_iter))
    Tau<-array(NA_real_,dim=c(K,K,N,T_data-1,n_iter))
    theta<-matrix(NA_real_,K,n_iter)
    
    ## Assigning the starting values
    gamma[,,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    pi[,,,1]<-start[[3]]
    theta[,1]<-start[[4]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<120)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Updating Tau(i,k,k_dash)
      Tau[,,,,iter_index]<-Tau.update.wrapper(gamma.curr=gamma[,,,iter_index-1],pi.curr=pi[,,,iter_index-1], theta.curr=theta[,iter_index-1],network=network,N=N,K=K,T_data = T_data)
      
      ## Updating the intial distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,,iter_index-1],N,K,T_data)
      
      ## Updating gamma[,,1] i.e. for first time point
      gamma[,,1,iter_index]<-gamma_init.update.wrapper(gamma_init=gamma[,,1,iter_index-1], log_alpha=log(alpha[,iter_index]), theta=theta[,iter_index-1], network_init=network[,,1], N=N, K=K)
      
      ## Converting Tau from a 4-dim array to a list
      Tau_cube<-array(NA_real_,dim=c(K,K,N*(T_data-1)))
      k<-1
      for (t in 1:(T_data-1)){
        for (i in 1:N){
          Tau_cube[,,k]<-Tau[,,i,t,iter_index]
          k<-k+1
        }
      }
      
      
      ## Updating gamma for time points >=2
      gamma[,,,iter_index]<-gamma_update_HMM_undir(gamma_init=gamma[,,1,iter_index], Tau=Tau_cube, N=N, K=K, T_data=T_data)
      
      ## Updating the K*K matrix pi
      pi[,,,iter_index]<-pi.update(gamma.curr = gamma[,,,iter_index],Tau.curr = Tau[,,,,iter_index],N = N,K = K,T_data = T_data)
      
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1], pi=pi[,,iter_index],gamma=gamma[,,,iter_index], network=network, N=N, K=K, T_data = T_data)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_undir(gamma=gamma[,,,iter_index], alpha=alpha[,iter_index], pi=pi[,,,iter_index], Tau=Tau_cube, theta = theta[,iter_index], network=network, N=N, K=K, T_data = T_data)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,alpha,pi,Tau,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N,T_data){
    gradient<-grad_HMM_undir_K1(theta=theta.curr, network=network, N=N, T_data=T_data)
    hess<-hess_HMM_undir_K1(theta=theta.curr, N=N, T_data=T_data)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N,T_data=T_data)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_undir_K1(theta = theta[iter_index], network=network, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i,t]
            cluster_id_j<-cluster_ids_est[j,t]
            exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
            t1<-t1+((network[i,j,t]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
          }
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i,1]])
      }
      ## 3rd term
      t3<-0
      for (t in 2:T_data){
        for (i in 1:N){
          t3<-t3+log(pi[cluster_ids_est[i,t-1],cluster_ids_est[i,t],t-1])
        }
      }
      comp_val<-t1+t2+t3
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            comp_val<-comp_val+((network[i,j,t]*(2*theta))-log_exp_val)
          }
        }
      }
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, pi = pi, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-K*(K-1)*(T_data-1)*log(N)
      t3<-K*log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2-t3
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K, T_data = T_data)
      t2<-log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  RASE_pi<-function(pi_est, pi_true){
    #diff_2<-rep(NA_real_,dim(pi_est)[3])
    diff_F<-rep(NA_real_,dim(pi_est)[3])
    for (t in 1:(dim(pi_est)[3])){
      #diff_2[t]<-norm((pi_est[,,t]-pi_true[,,t]),"2")
      diff_F[t]<-norm((pi_est[,,t]-pi_true[,,t]),"F")
    }
    #RASE_2<-mean(diff_2)
    RASE_F<-mean(diff_F)
    return(RASE_F)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  K<-nclust ## Defining the number of clusters
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-array(rep(1/K,N*K*T_data),dim=c(N,K,T_data))
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-array(rep(1/K,K^2*(T_data-1)),dim=c(K,K,T_data-1)) ## pi
    start[[4]]<-theta_init ## theta
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,param[1:n_last],ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,,n_last]
    param_converge[[4]]<-param[[4]][,,,,n_last]
    param_converge[[5]]<-param[[5]][,n_last]
    cluster_ids_est<-matrix(NA_integer_,N,T_data)
    for (t in 1:T_data){
      cluster_ids_est[,t]<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        cluster_id<-which.max(param_converge[[1]][,,t][x,])
        return(cluster_id)
      }))
    }
    cluster_ids_true_matrix<-matrix(rep(cluster_ids_true,T_data),N,T_data)
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] , pi = param_converge[[3]],theta = param_converge[[5]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true_matrix)
        Rand_val_final<-ARI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        RASE_pi_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[5]][K_permute_mat[k,]]
          #pi_est<-array(NA_real_,dim=c(K,K,T_data-1))
          for (t in 1:(T_data-1)){
            pi_true[,,t]<-as(as.integer(K_permute_mat[k,]), "pMatrix")%*%pi_true[,,t]%*%t(as(as.integer(K_permute_mat[k,]), "pMatrix"))
          }
          
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
          RASE_pi_vec[k]<-RASE_pi(pi_est = param_converge[[3]],pi_true = pi_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        permute_true_id_pi<-which.min(RASE_pi_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        RASE_pi_val<-RASE_pi_vec[permute_true_id_pi]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val,RASE_pi_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for HMM undirected stability case to run for different number of clusters 
wrapper_HMM_dyn_undir_Stab<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=1,sim_indicator_TERGM=0,theta_true=NA,pi_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(HMMEMundirStab)
  library(Matrix)
  
  #################################################
  ## Defining a function to update variational parameters Tau 
  Tau.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_data){
    Tau.next.temp<-Tau_update_HMM_undir_Stab(gamma=gamma.curr, log_pi=log(pi.curr), theta=theta.curr, network=network, N=N, K=K, T_data = T_data)
    ## converting in 4-dimensional array and normalizing Tau(i,k,)
    Tau.next<-array(NA_real_,dim=c(K,K,N,T_data-2))
    for (t in 1:(T_data-2)){
      for (i in 1:N){
        for (k in 1:K){
          Tau.next[k,,i,t]<-Tau.next.temp[[N*(t-1)+i]][k,]/sum(Tau.next.temp[[N*(t-1)+i]][k,])
        }
      }
    }
    
    return(Tau.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  gamma_init.update.wrapper<-function(gamma_init, log_alpha, theta, network_12, N, K){
    gamma_init.next<-gamma_init_update_HMM_undir_Stab(gamma_init=gamma_init, log_alpha=log_alpha, theta=theta, network_12=network_12, N=N, K=K)
    ## Normalization of gamma_init
    for (i in 1:N){
      gamma_init.next[i,]<-gamma_init.next[i,]/sum(gamma_init.next[i,])
    }
    return(gamma_init.next)
  }
  
  #################################################
  ## Defining a function to update K*K matrix of pi
  pi.update<-function(gamma.curr,Tau.curr,N,K,T_data){
    pi.next<-array(NA_real_,dim=c(K,K,T_data-2))
    for (t in 1:(T_data-2)){
      for (k in 1:K){
        for (l in 1:K){
          pi.next[k,l,t]<-as.numeric((gamma.curr[,k,t]%*%Tau.curr[k,l,,t]))
          #denom<-denom+sum(gamma.curr[,k,t-1])
        }
        pi.next[k,,t]<-pi.next[k,,t]/sum(pi.next[k,,t]) ## Normalization
      }
    }
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_HMM_undir_Stab(theta=as.vector(theta.curr), gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    hess<-hess_HMM_undir_Stab(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,T_data-1,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    pi<-array(NA_real_,dim=c(K,K,T_data-2,n_iter))
    Tau<-array(NA_real_,dim=c(K,K,N,T_data-2,n_iter))
    theta<-matrix(NA_real_,K,n_iter)
    
    ## Assigning the starting values
    gamma[,,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    pi[,,,1]<-start[[3]]
    theta[,1]<-start[[4]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Updating Tau(i,k,k_dash)
      Tau[,,,,iter_index]<-Tau.update.wrapper(gamma.curr=gamma[,,,iter_index-1],pi.curr=pi[,,,iter_index-1], theta.curr=theta[,iter_index-1],network=network,N=N,K=K,T_data = T_data)
      
      ## Updating the intial distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,,iter_index-1],N,K,T_data)
      
      ## Updating gamma[,,1] i.e. for first time point
      gamma[,,1,iter_index]<-gamma_init.update.wrapper(gamma_init=gamma[,,1,iter_index-1], log_alpha=log(alpha[,iter_index]), theta=theta[,iter_index-1], network_12=network[,,c(1,2)], N=N, K=K)
      
      ## Converting Tau from a 4-dim array to a list
      Tau_cube<-array(NA_real_,dim=c(K,K,N*(T_data-2)))
      k<-1
      for (t in 1:(T_data-2)){
        for (i in 1:N){
          Tau_cube[,,k]<-Tau[,,i,t,iter_index]
          k<-k+1
        }
      }
      
      ## Updating gamma for time points >=2
      gamma[,,,iter_index]<-gamma_update_HMM_undir_Stab(gamma_init=gamma[,,1,iter_index], Tau=Tau_cube, N=N, K=K, T_data=T_data)
      
      ## Updating the K*K matrix pi
      pi[,,,iter_index]<-pi.update(gamma.curr = gamma[,,,iter_index],Tau.curr = Tau[,,,,iter_index],N = N,K = K,T_data = T_data)
      
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1], pi=pi[,,iter_index],gamma=gamma[,,,iter_index], network=network, N=N, K=K, T_data = T_data)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_undir_Stab(gamma=gamma[,,,iter_index], alpha=alpha[,iter_index], pi=pi[,,,iter_index], Tau=Tau_cube, theta = theta[,iter_index], network=network, N=N, K=K, T_data = T_data)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,alpha,pi,Tau,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N,T_data){
    gradient<-grad_HMM_undir_K1_Stab(theta=theta.curr, network=network, N=N, T_data=T_data)
    hess<-hess_HMM_undir_K1_Stab(theta=theta.curr, N=N, T_data=T_data)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N,T_data=T_data)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_undir_K1_Stab(theta = theta[iter_index], network=network, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-0
      for (t in 2:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i,t-1]
            cluster_id_j<-cluster_ids_est[j,t-1]
            Stab_stat<-network[i,j,t]*network[i,j,t-1]+((1-network[i,j,t])*(1-network[i,j,t-1]))
            exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
            t1<-t1+((Stab_stat*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
          }
        }
      }
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i,1]])
      }
      t3<-0
      for (t in 3:T_data){
        for (i in 1:N){
          t3<-t3+log(pi[cluster_ids_est[i,t-2],cluster_ids_est[i,t-1],t-2])
        }
      }
      comp_val<-t1+t2+t3
    }else if(K==1){
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      comp_val<-0
      for (t in 2:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            Stab_stat<-network[i,j,t]*network[i,j,t-1]+((1-network[i,j,t])*(1-network[i,j,t-1]))
            comp_val<-comp_val+((Stab_stat*2*theta)-log_exp_val)
          }
        }
      }
    }
    return(comp_val)
  }
  
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, pi = pi, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-K*(K-1)*(T_data-2)*log(N)
      t3<-K*log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2-t3
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N, K=K, T_data = T_data)
      t2<-log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  RASE_pi<-function(pi_est, pi_true){
    K<-dim(pi_est)[1]
    diff_F<-rep(NA_real_,dim(pi_est)[3])
    for (t in 1:(dim(pi_est)[3])){
      diff_F[t]<-sum((pi_est-pi_true)^2)/(K^2)
    }
    RASE_F<-mean(diff_F)
    return(RASE_F)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  K<-nclust ## Defining the number of clusters
  gamma.start<-array(rep(1/K,N*K*(T_data-1)),dim=c(N,K,T_data-1))
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-array(rep(1/K,K^2*(T_data-2)),dim=c(K,K,T_data-2))
    start[[4]]<-theta_init ## theta
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K=K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,,n_last]
    param_converge[[4]]<-param[[4]][,,,,n_last]
    param_converge[[5]]<-param[[5]][,n_last]
    cluster_ids_est<-matrix(NA_integer_,N,T_data-1)
    for (t in 1:(T_data-1)){
      cluster_ids_est[,t]<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        cluster_id<-which.max(param_converge[[1]][,,t][x,])
        return(cluster_id)
      }))
    }
    cluster_ids_true_matrix<-matrix(rep(cluster_ids_true,T_data-1),N,T_data-1)
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] , pi = param_converge[[3]],theta = param_converge[[5]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true_matrix)
        Rand_val_final<-ARI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true[,2:10])
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[5]][K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        for (t in 1:(T_data-2)){
          param_converge[[3]][,,t]<-as(as.integer(K_permute_mat[permute_true_id_theta,]), "pMatrix")%*%param_converge[[3]][,,t]%*%t(as(as.integer(K_permute_mat[permute_true_id_theta,]), "pMatrix"))
        }
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        RASE_pi_val<-RASE_pi(pi_est = param_converge[[3]],pi_true = pi_true)
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val,RASE_pi_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for HMM directed density case to run for different number of clusters 
wrapper_HMM_dyn_dir_Dens<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=1,sim_indicator_TERGM=0,theta_true=NA,pi_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(HMMEMdir)
  library(Matrix)
  
  #################################################
  ## Defining a function to update variational parameters Tau 
  Tau.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K,T_data){
    Tau.next.temp<-Tau_update_HMM_dir(gamma=gamma.curr, log_pi=log(pi.curr), theta=theta.curr, network=network, N=N, K=K, T_data = T_data)
    ## converting in 4-dimensional array and normalizing Tau(i,k,)
    Tau.next<-array(NA_real_,dim=c(K,K,N,T_data-1))
    for (t in 1:(T_data-1)){
      for (i in 1:N){
        for (k in 1:K){
          Tau.next[k,,i,t]<-Tau.next.temp[[N*(t-1)+i]][k,]/sum(Tau.next.temp[[N*(t-1)+i]][k,])
        }
      }
    }
    
    return(Tau.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  gamma_init.update.wrapper<-function(gamma_init, log_alpha, theta, network_init, N, K){
    gamma_init.next<-gamma_init_update_HMM_dir(gamma_init=gamma_init, log_alpha=log_alpha, theta=theta, network_init=network_init, N=N, K=K)
    ## Normalization of gamma_init
    for (i in 1:N){
      gamma_init.next[i,]<-gamma_init.next[i,]/sum(gamma_init.next[i,])
    }
    return(gamma_init.next)
  }
  
  #################################################
  ## Defining a function to update K*K matrix of pi
  pi.update<-function(gamma.curr,Tau.curr,N,K,T_data){
    pi.next<-array(NA_real_,dim=c(K,K,T_data-1))
    for (t in 2:T_data){
      for (k in 1:K){
        for (l in 1:K){
          pi.next[k,l,t-1]<-as.numeric((gamma.curr[,k,t-1]%*%Tau.curr[k,l,,t-1]))
          #denom<-denom+sum(gamma.curr[,k,t-1])
        }
        pi.next[k,,t-1]<-pi.next[k,,t-1]/sum(pi.next[k,,t-1]) ## Normalization
      }
    }
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K,T_data){
    theta.next<-matrix(NA_real_,K,2)
    gradient_oe<-grad_HMM_dir_oe(theta=theta.curr, gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    gradient_re<-grad_HMM_dir_re(theta=theta.curr, gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    #gradient<-c(gradient_oe,gradient_re)
    hess_oe<-hess_HMM_dir_oe(theta=theta.curr, gamma=gamma, N=N, K=K, T_data=T_data)
    hess_re<-hess_HMM_dir_re(theta=theta.curr, gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next[,1]<-theta.curr[,1]-as.vector(solve(hess_oe+diag((10^(-6)),K))%*%gradient_oe)
    theta.next[,2]<-theta.curr[,2]-as.vector(solve(hess_re+diag((10^(-6)),K))%*%gradient_re)
    # hess_oe_re<-hess_HMM_dir_oe_re(theta=theta.curr, gamma=gamma, N=N, K=K, T_data=T_data)
    # hess<-matrix(NA_real_,2*K,2*K)
    # hess[1:K,1:K]<-hess_oe
    # hess[((K+1):(2*K)),((K+1):(2*K))]<-hess_re
    # hess[(1:K),((K+1):(2*K))]<-hess_oe_re
    # hess[((K+1):(2*K)),(1:K)]<-t(hess_oe_re)
    # theta.next<-matrix(c(theta.curr[,1],theta.curr[,2])-as.vector(solve(hess)%*%gradient),K,2)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,T_data,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    pi<-array(NA_real_,dim=c(K,K,T_data-1,n_iter))
    Tau<-array(NA_real_,dim=c(K,K,N,T_data-1,n_iter))
    theta<-array(NA_real_,dim=c(K,2,n_iter))
    
    ## Assigning the starting values
    gamma[,,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    pi[,,,1]<-start[[3]]
    theta[,,1]<-start[[4]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Updating Tau(i,k,k_dash)
      Tau[,,,,iter_index]<-Tau.update.wrapper(gamma.curr=gamma[,,,iter_index-1],pi.curr=pi[,,,iter_index-1], theta.curr=theta[,,iter_index-1],network=network,N=N,K=K,T_data = T_data)
      
      ## Updating the intial distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,,iter_index-1],N,K,T_data)
      
      ## Updating gamma[,,1] i.e. for first time point
      gamma[,,1,iter_index]<-gamma_init.update.wrapper(gamma_init=gamma[,,1,iter_index-1], log_alpha=log(alpha[,iter_index]), theta=theta[,,iter_index-1], network_init=network[,,1], N=N, K=K)
      
      ## Converting Tau from a 4-dim array to a list
      Tau_cube<-array(NA_real_,dim=c(K,K,N*(T_data-1)))
      k<-1
      for (t in 1:(T_data-1)){
        for (i in 1:N){
          Tau_cube[,,k]<-Tau[,,i,t,iter_index]
          k<-k+1
        }
      }
      
      
      ## Updating gamma for time points >=2
      gamma[,,,iter_index]<-gamma_update_HMM_dir(gamma_init=gamma[,,1,iter_index], Tau=Tau_cube, N=N, K=K, T_data=T_data)
      
      ## Updating the K*K matrix pi
      pi[,,,iter_index]<-pi.update(gamma.curr = gamma[,,,iter_index],Tau.curr = Tau[,,,,iter_index],N = N,K = K,T_data = T_data)
      
      ## Updating the theta vector
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1], pi=pi[,,iter_index],gamma=gamma[,,,iter_index], network=network, N=N, K=K, T_data = T_data)
      #print(theta[,,iter_index])
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_dir(gamma=gamma[,,,iter_index], alpha=alpha[,iter_index], pi=pi[,,,iter_index], Tau=Tau_cube, theta = theta[,,iter_index], network=network, N=N, K=K, T_data = T_data)
      print(ELBO_grid.curr)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,alpha,pi,Tau,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  theta.update_K1<-function(theta.curr,network,N,T_data){
    gradient_oe<-grad_HMM_dir_oe_K1(theta=theta.curr, network=network, N=N, T_data=T_data)
    gradient_re<-grad_HMM_dir_re_K1(theta=theta.curr, network=network, N=N, T_data=T_data)
    hess_oe<-hess_HMM_dir_oe_K1(theta=theta.curr, N=N, T_data=T_data)
    hess_re<-hess_HMM_dir_re_K1(theta=theta.curr, N=N, T_data=T_data)
    hess_oe_re<-hess_HMM_dir_oe_re_K1(theta=theta.curr, N=N, T_data=T_data)
    gradient<-c(gradient_oe,gradient_re)
    hess_mat<-matrix(c(hess_oe,hess_oe_re,hess_oe_re,hess_re),2,2)
    theta.next<-theta.curr-as.vector(solve(hess_mat)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    
    ## initializing the arrays for parameters
    theta<-matrix(NA_real_,2,n_iter)
    theta[,1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], network=network, N=N, T_data=T_data)
      print(theta[,iter_index])
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_dir_K1(theta = theta[,iter_index], network=network, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  #########################################################################################################
  #########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-0
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i,t]
            cluster_id_j<-cluster_ids_est[j,t]
            exp_val_1=exp(theta[cluster_id_i,1])
            exp_val_2=exp(theta[cluster_id_j,1])
            exp_val_3=exp(theta[cluster_id_i,2]+theta[cluster_id_j,2])
            indicator_10=(network[i,j,t]==1)&(network[j,i,t]==0)
            indicator_01=(network[i,j,t]==0)&(network[j,i,t]==1)
            indicator_11=(network[i,j,t]==1)&(network[j,i,t]==1)
            t1<-t1+((indicator_10*theta[cluster_id_i,1])+(indicator_01*theta[cluster_id_j,1])+(indicator_11*(theta[cluster_id_i,2]+theta[cluster_id_j,2]))-log(1+exp_val_1+exp_val_2+exp_val_3))
          }
        }
      }
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i,1]])
      }
      t3<-0
      for (t in 2:T_data){
        for (i in 1:N){
          t3<-t3+log(pi[cluster_ids_est[i,t-1],cluster_ids_est[i,t],t-1])
        }
      }
      comp_val<-t1+t2+t3
    }else if(K==1){
      comp_val<-0
      exp_val_1=exp(theta[1])
      exp_val_2=exp(2*theta[2])
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            indicator_10=(network[i,j,t]==1)&(network[j,i,t]==0)
            indicator_01=(network[i,j,t]==0)&(network[j,i,t]==1)
            indicator_11=(network[i,j,t]==1)&(network[j,i,t]==1)
            comp_val<-comp_val+((indicator_10*theta[1])+(indicator_01*theta[1])+(indicator_11*(2*theta[2]))-log(1+2*exp_val_1+exp_val_2))
          }
        }
      }
    }
    return(comp_val)
  }
  
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, pi = pi, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-K*(K-1)*(T_data-1)*log(N)
      t3<-2*K*log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2-t3
    }else if(K==1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, pi = pi, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-2*log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  RASE_pi<-function(pi_est, pi_true){
    #diff_2<-rep(NA_real_,dim(pi_est)[3])
    diff_F<-rep(NA_real_,dim(pi_est)[3])
    for (t in 1:(dim(pi_est)[3])){
      #diff_2[t]<-norm((pi_est[,,t]-pi_true[,,t]),"2")
      diff_F[t]<-norm((pi_est[,,t]-pi_true[,,t]),"F")
    }
    #RASE_2<-mean(diff_2)
    RASE_F<-mean(diff_F)
    return(RASE_F)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  
  K<-nclust ## Defining the number of clusters
  
  #################################################
  gamma.start<-array(rep(1/K,N*K*T_data),dim=c(N,K,T_data))
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    #debug(iterator_K1)
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-array(rep(1/K,K^2*(T_data-1)),dim=c(K,K,T_data-1)) ## pi
    start[[4]]<-theta_init
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K=K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,,n_last]
    param_converge[[4]]<-param[[4]][,,,,n_last]
    param_converge[[5]]<-param[[5]][,,n_last]
    cluster_ids_est<-matrix(NA_integer_,N,T_data)
    for (t in 1:T_data){
      cluster_ids_est[,t]<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        cluster_id<-which.max(param_converge[[1]][,,t][x,])
        return(cluster_id)
      }))
    }
    cluster_ids_true_matrix<-matrix(rep(cluster_ids_true,T_data),N,T_data)
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] , pi = param_converge[[3]],theta = param_converge[[5]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true_matrix)
        Rand_val_final<-ARI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[5]][K_permute_mat[k,],]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id]
        pi_true_rotated<-array(NA_real_,dim=c(K,K,T_data-1))
        for (t in 1:(T_data-1)){
          pi_true_rotated[,,t]<-as(as.integer(K_permute_mat[permute_true_id,]), "pMatrix")%*%pi_true[,,t]%*%t(as(as.integer(K_permute_mat[permute_true_id,]), "pMatrix"))
        }
        RASE_pi_val<-RASE_pi(pi_est = param_converge[[3]],pi_true = pi_true_rotated)
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val,RASE_pi_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for HMM directed transitivity case to run for different number of clusters
wrapper_HMM_dyn_dir_Trans<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=1,sim_indicator_TERGM=0,theta_true=NA,pi_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(HMMEMdirTrans)
  library(Matrix)
  
  #################################################
  ## Defining a function to update variational parameters Tau 
  Tau.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,Trans_stat=Trans_stat,N,K,T_data){
    Tau.next.temp<-Tau_update_HMM_dir_Trans(gamma=gamma.curr, log_pi=log(pi.curr), theta=theta.curr, network=network, Trans_stat=Trans_stat,N=N, K=K, T_data = T_data)
    ## converting in 4-dimensional array and normalizing Tau(i,k,)
    Tau.next<-array(NA_real_,dim=c(K,K,N,T_data-2))
    for (t in 1:(T_data-2)){
      for (i in 1:N){
        for (k in 1:K){
          Tau.next[k,,i,t]<-Tau.next.temp[[N*(t-1)+i]][k,]/sum(Tau.next.temp[[N*(t-1)+i]][k,])
        }
      }
    }
    
    return(Tau.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  gamma_init.update.wrapper<-function(gamma_init, log_alpha, theta, network_12, Trans_stat_12, N, K){
    gamma_init.next<-gamma_init_update_HMM_dir_Trans(gamma_init=gamma_init, log_alpha=log_alpha, theta=theta, network_12=network_12, Trans_stat_12=Trans_stat_12, N=N, K=K)
    ## Normalization of gamma_init
    for (i in 1:N){
      gamma_init.next[i,]<-gamma_init.next[i,]/sum(gamma_init.next[i,])
    }
    return(gamma_init.next)
  }
  
  #################################################
  ## Defining a function to update K*K matrix of pi
  pi.update<-function(gamma.curr,Tau.curr,N,K,T_data){
    pi.next<-array(NA_real_,dim=c(K,K,T_data-2))
    for (t in 1:(T_data-2)){
      for (k in 1:K){
        for (l in 1:K){
          pi.next[k,l,t]<-as.numeric((gamma.curr[,k,t]%*%Tau.curr[k,l,,t]))
          #denom<-denom+sum(gamma.curr[,k,t-1])
        }
        pi.next[k,,t]<-pi.next[k,,t]/sum(pi.next[k,,t]) ## Normalization
      }
    }
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,Trans_stat=Trans_stat,N,K,T_data){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_HMM_dir_Trans(theta=as.vector(theta.curr), gamma=gamma, network=network, Trans_stat=Trans_stat, N=N, K=K, T_data=T_data)
    hess<-hess_HMM_dir_Trans(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,T_data-1,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    pi<-array(NA_real_,dim=c(K,K,T_data-2,n_iter))
    Tau<-array(NA_real_,dim=c(K,K,N,T_data-2,n_iter))
    theta<-matrix(NA_real_,K,n_iter)
    
    ## Assigning the starting values
    gamma[,,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    pi[,,,1]<-start[[3]]
    theta[,1]<-start[[4]]
    
    Trans_stat<-Trans_stat_cal(network=network,N=N,T_data=T_data)
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Updating Tau(i,k,k_dash)
      Tau[,,,,iter_index]<-Tau.update.wrapper(gamma.curr=gamma[,,,iter_index-1],pi.curr=pi[,,,iter_index-1], theta.curr=theta[,iter_index-1],network=network,Trans_stat=Trans_stat,N=N,K=K,T_data = T_data)
      
      ## Updating the intial distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,,iter_index-1],N,K,T_data)
      
      ## Updating gamma[,,1] i.e. for first time point
      gamma[,,1,iter_index]<-gamma_init.update.wrapper(gamma_init=gamma[,,1,iter_index-1], log_alpha=log(alpha[,iter_index]), theta=theta[,iter_index-1], network_12=network[,,c(1,2)], Trans_stat_12=Trans_stat[,,1], N=N, K=K)
      
      ## Converting Tau from a 4-dim array to a list
      Tau_cube<-array(NA_real_,dim=c(K,K,N*(T_data-2)))
      k<-1
      for (t in 1:(T_data-2)){
        for (i in 1:N){
          Tau_cube[,,k]<-Tau[,,i,t,iter_index]
          k<-k+1
        }
      }
      
      ## Updating gamma for time points >=2
      gamma[,,,iter_index]<-gamma_update_HMM_dir_Trans(gamma_init=gamma[,,1,iter_index], Tau=Tau_cube, N=N, K=K, T_data=T_data)
      
      ## Updating the K*K matrix pi
      pi[,,,iter_index]<-pi.update(gamma.curr = gamma[,,,iter_index],Tau.curr = Tau[,,,,iter_index],N = N,K = K,T_data = T_data)
      
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1], pi=pi[,,iter_index],gamma=gamma[,,,iter_index], network=network, Trans_stat=Trans_stat, N=N, K=K, T_data = T_data)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_dir_Trans(gamma=gamma[,,,iter_index], alpha=alpha[,iter_index], pi=pi[,,,iter_index], Tau=Tau_cube, theta = theta[,iter_index], network=network, Trans_stat=Trans_stat, N=N, K=K, T_data = T_data)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,alpha,pi,Tau,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,Trans_stat,N,T_data){
    gradient<-grad_HMM_dir_K1_Trans(theta=theta.curr, network=network, Trans_stat=Trans_stat, N=N, T_data=T_data)
    hess<-hess_HMM_dir_K1_Trans(theta=theta.curr, N=N, T_data=T_data)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    Trans_stat<-Trans_stat_cal(network=network,N=N,T_data=T_data)
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, Trans_stat=Trans_stat, N=N,T_data=T_data)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_dir_K1_Trans(theta = theta[iter_index], network=network, Trans_stat=Trans_stat, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  #########################################################################################################
  #########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    Trans_stat<-Trans_stat_cal(network = network,N = N,T_data = T_data)
    if(K!=1){
      t1<-0
      for (t in 2:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i,t-1]
            cluster_id_j<-cluster_ids_est[j,t-1]
            exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
            t1<-t1+((Trans_stat[i,j,t-1]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
          }
        }
      }
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i,1]])
      }
      t3<-0
      for (t in 3:T_data){
        for (i in 1:N){
          t3<-t3+log(pi[cluster_ids_est[i,t-2],cluster_ids_est[i,t-1],t-2])
        }
      }
      comp_val<-t1+t2+t3
    }else if(K==1){
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      comp_val<-0
      for (t in 2:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            comp_val<-comp_val+((Trans_stat[i,j,t-1]*2*theta)-log_exp_val)
          }
        }
      }
    }
    return(comp_val)
  }
  
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,pi=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, pi = pi, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-K*(K-1)*(T_data-2)*log(N)
      t3<-K*log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2-t3
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K, T_data = T_data)
      t2<-log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  RASE_pi<-function(pi_est, pi_true){
    #diff_2<-rep(NA_real_,dim(pi_est)[3])
    diff_F<-rep(NA_real_,dim(pi_est)[3])
    for (t in 1:(dim(pi_est)[3])){
      #diff_2[t]<-norm((pi_est[,,t]-pi_true[,,t]),"2")
      diff_F[t]<-norm((pi_est[,,t]-pi_true[,,t]),"F")
    }
    #RASE_2<-mean(diff_2)
    RASE_F<-mean(diff_F)
    return(RASE_F)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  K<-nclust ## Defining the number of clusters
  gamma.start<-array(rep(1/K,N*K*(T_data-1)),dim=c(N,K,T_data-1))
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-array(rep(1/K,K^2*(T_data-2)),dim=c(K,K,T_data-2))
    start[[4]]<-theta_init ## theta
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K=K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,,n_last]
    param_converge[[4]]<-param[[4]][,,,,n_last]
    param_converge[[5]]<-param[[5]][,n_last]
    cluster_ids_est<-matrix(NA_integer_,N,T_data-1)
    for (t in 1:(T_data-1)){
      cluster_ids_est[,t]<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
        cluster_id<-which.max(param_converge[[1]][,,t][x,])
        return(cluster_id)
      }))
    }
    cluster_ids_true_matrix<-matrix(rep(cluster_ids_true,T_data),N,T_data)
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] , pi = param_converge[[3]],theta = param_converge[[5]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true_matrix)
        Rand_val_final<-ARI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        RASE_pi_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[5]][K_permute_mat[k,]]
          #pi_est<-array(NA_real_,dim=c(K,K,T_data-1))
          for (t in 1:(T_data-2)){
            pi_true[,,t]<-as(as.integer(K_permute_mat[k,]), "pMatrix")%*%pi_true[,,t]%*%t(as(as.integer(K_permute_mat[k,]), "pMatrix"))
          }
          
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
          RASE_pi_vec[k]<-RASE_pi(pi_est = param_converge[[3]],pi_true = pi_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        permute_true_id_pi<-which.min(RASE_pi_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        RASE_pi_val<-RASE_pi_vec[permute_true_id_pi]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val,RASE_pi_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
#########################################################################################################
##################################   Dynamic network code (TERGM)   #######################################
#########################################################################################################
#########################################################################################################
## Defining a wrapper function for TERGM undirected density case to run for different number of clusters 
wrapper_TERGM_undir_Dens<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=0,sim_indicator_TERGM=1,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(quadprog)
  library(TERGMundir)
  library(Matrix)
  #sourceCpp("TERGMundir.cpp")
  
  #################################################
  gamma.update.wrapper<-function(gamma.curr,alpha.curr,theta.curr,network,N,K,T_data){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_TERGM_undir(gamma=gamma.curr, alpha=alpha.curr, theta=theta.curr, network=network, N=N, K=K, T_data=T_data)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,gamma,network,N,K,T_data){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_TERGM_undir(theta=as.vector(theta.curr), gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    hess<-hess_TERGM_undir(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    
    ## Assigning the starting values
    gamma[,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Updating the intial distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,iter_index-1],N,K,T_data)
      
      ## Updating gamma
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr = gamma[,,iter_index-1],alpha.curr = alpha[,iter_index],theta.curr = theta[,iter_index-1],network = network,N = N,K = K,T_data = T_data)
      
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data = T_data)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_undir(gamma=gamma[,,iter_index], alpha=alpha[,iter_index], theta = theta[,iter_index], network=network, N=N, K=K, T_data = T_data)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,alpha,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N,T_data){
    gradient<-grad_TERGM_undir_K1(theta=theta.curr, network=network, N=N, T_data=T_data)
    hess<-hess_TERGM_undir_K1(theta=theta.curr, N=N, T_data=T_data)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while(error>thres){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N,T_data=T_data)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_undir_K1(theta = theta[iter_index], network=network, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i]
            cluster_id_j<-cluster_ids_est[j]
            exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
            t1<-t1+((network[i,j,t]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
          }
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            comp_val<-comp_val+((network[i,j,t]*(2*theta))-log_exp_val)
          }
        }
      }
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-K*log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K, T_data = T_data)
      t2<-log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  K<-nclust ## Defining the number of clusters
  
  #################################################
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (mixture proportions)
    start[[3]]<-theta_init ## theta
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,param[1:n_last],ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    cluster_ids_est_matrix<-matrix(rep(cluster_ids_est,T_data),N,T_data)
    
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]], theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-RI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est_matrix,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,]]
          #pi_est<-array(NA_real_,dim=c(K,K,T_data-1))
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for TERGM directed density case to run for different number of clusters 
wrapper_TERGM_dir_Dens<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=0,sim_indicator_TERGM=1,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(quadprog)
  library(TERGMdir)
  library(Matrix)
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K*2 array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,alpha.curr,theta.curr,network,N,K,T_data){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    quad_lin_coeff<-gamma_update_TERGM_dir(gamma=gamma.curr, pi=alpha.curr, theta=theta.curr, network=network, N=N, K=K, T_data=T_data)
    
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,gamma,network,N,K,T_data){
    theta.next<-matrix(NA_real_,K,2)
    gradient_oe<-grad_TERGM_dir_oe(theta=theta.curr, gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    gradient_re<-grad_TERGM_dir_re(theta=theta.curr, gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    #gradient<-c(gradient_oe,gradient_re)
    hess_oe<-hess_TERGM_dir_oe(theta=theta.curr, gamma=gamma, N=N, K=K, T_data=T_data)
    hess_re<-hess_TERGM_dir_re(theta=theta.curr, gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next[,1]<-theta.curr[,1]-as.vector(solve(hess_oe+diag((10^(-6)),K))%*%gradient_oe)
    theta.next[,2]<-theta.curr[,2]-as.vector(solve(hess_re+diag((10^(-6)),K))%*%gradient_re)
    # hess_oe_re<-hess_HMM_dir_oe_re(theta=theta.curr, gamma=gamma, N=N, K=K, T_data=T_data)
    # hess<-matrix(NA_real_,2*K,2*K)
    # hess[1:K,1:K]<-hess_oe
    # hess[((K+1):(2*K)),((K+1):(2*K))]<-hess_re
    # hess[(1:K),((K+1):(2*K))]<-hess_oe_re
    # hess[((K+1):(2*K)),(1:K)]<-t(hess_oe_re)
    # theta.next<-matrix(c(theta.curr[,1],theta.curr[,2])-as.vector(solve(hess)%*%gradient),K,2)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    theta<-array(NA_real_,dim=c(K,2,n_iter))
    
    ## Assigning the starting values
    gamma[,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    theta[,,1]<-start[[3]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the mixture distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,iter_index-1],N = N,K = K,T_data = T_data)
      
      ## Updating gamma
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr = gamma[,,iter_index-1],alpha.curr = alpha[,iter_index],theta.curr = theta[,,iter_index-1],network = network,N = N,K = K,T_data = T_data)
      
      ## Updating the theta matrix (outgoing and reciprocity)
      theta[,,iter_index]<-theta.update(theta.curr=theta[,,iter_index-1],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data = T_data)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_dir(gamma=gamma[,,iter_index], pi=alpha[,iter_index],theta = theta[,,iter_index], network=network, N=N, K=K, T_data = T_data)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    print(theta[,,iter_index-1])
    return(list(gamma,alpha,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  theta.update_K1<-function(theta.curr,network,N,T_data){
    gradient_oe<-grad_TERGM_dir_K1_oe(theta=theta.curr, network=network, N=N, T_data=T_data)
    gradient_re<-grad_TERGM_dir_K1_re(theta=theta.curr, network=network, N=N, T_data=T_data)
    hess_oe<-hess_TERGM_dir_K1_oe(theta=theta.curr, N=N, T_data=T_data)
    hess_re<-hess_TERGM_dir_K1_re(theta=theta.curr, N=N, T_data=T_data)
    hess_oe_re<-hess_TERGM_dir_K1_oe_re(theta=theta.curr, N=N, T_data=T_data)
    gradient<-c(gradient_oe,gradient_re)
    hess_mat<-matrix(c(hess_oe,hess_oe_re,hess_oe_re,hess_re),2,2)
    theta.next<-theta.curr-as.vector(solve(hess_mat)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    
    ## initializing the arrays for parameters
    theta<-matrix(NA_real_,2,n_iter)
    theta[,1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update_K1(theta.curr=theta[,iter_index-1], network=network, N=N, T_data=T_data)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_dir_K1(theta = theta[,iter_index], network=network, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    print(theta[,iter_index-1])
    return(theta)
  }
  
  #########################################################################################################
  #########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-0
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i]
            cluster_id_j<-cluster_ids_est[j]
            exp_val_1=exp(theta[cluster_id_i,1])
            exp_val_2=exp(theta[cluster_id_j,1])
            exp_val_3=exp(theta[cluster_id_i,2]+theta[cluster_id_j,2])
            indicator_10=(network[i,j,t]==1)&(network[j,i,t]==0)
            indicator_01=(network[i,j,t]==0)&(network[j,i,t]==1)
            indicator_11=(network[i,j,t]==1)&(network[j,i,t]==1)
            t1<-t1+((indicator_10*theta[cluster_id_i,1])+(indicator_01*theta[cluster_id_j,1])+(indicator_11*(theta[cluster_id_i,2]+theta[cluster_id_j,2]))-log(1+exp_val_1+exp_val_2+exp_val_3))
          }
        }
      }
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val_1=exp(theta[1])
      exp_val_2=exp(2*theta[2])
      for (t in 1:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            indicator_10=(network[i,j,t]==1)&(network[j,i,t]==0)
            indicator_01=(network[i,j,t]==0)&(network[j,i,t]==1)
            indicator_11=(network[i,j,t]==1)&(network[j,i,t]==1)
            comp_val<-comp_val+((indicator_10*theta[1])+(indicator_01*theta[1])+(indicator_11*(2*theta[2]))-log(1+2*exp_val_1+exp_val_2))
          }
        }
      }
    }
    return(comp_val)
  }
  
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-2*K*log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-2*log((N*(N-1)*(T_data))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  K<-nclust ## Defining the number of clusters
  
  #################################################
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (mixture proportions)
    start[[3]]<-theta_init ## theta
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[1,n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[,n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    cluster_ids_est_matrix<-matrix(rep(cluster_ids_est,T_data),N,T_data)
    
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]], theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-RI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est_matrix,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,],]
          #pi_est<-array(NA_real_,dim=c(K,K,T_data-1))
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for HMM undirected stability case to run for different number of clusters 
wrapper_TERGM_undir_Stab<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=0,sim_indicator_TERGM=1,theta_true=NA,pi_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(quadprog)
  library(TERGMundirStab)
  library(Matrix)
  
  #################################################
  gamma.update.wrapper<-function(gamma.curr,alpha.curr,theta.curr,network,N,K,T_data){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_TERGM_undir_Stab(gamma=gamma.curr, pi=alpha.curr, theta=theta.curr, network=network, N=N, K=K, T_data=T_data)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,gamma,network,N,K,T_data){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_TERGM_undir_Stab(theta=as.vector(theta.curr), gamma=gamma, network=network, N=N, K=K, T_data=T_data)
    hess<-hess_TERGM_undir_Stab(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    
    ## Assigning the starting values
    gamma[,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<120)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Updating the intial distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,iter_index-1],N = N,K = K,T_data = T_data)
      
      ## Updating gamma
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr = gamma[,,iter_index-1],alpha.curr = alpha[,iter_index],theta.curr = theta[,iter_index-1],network = network,N = N,K = K,T_data = T_data)
      
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1],gamma=gamma[,,iter_index], network=network, N=N, K=K, T_data = T_data)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_undir_Stab(gamma=gamma[,,iter_index], alpha=alpha[,iter_index],theta = theta[,iter_index], network=network, N=N, K=K, T_data = T_data)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    print(theta[,iter_index-1])
    return(list(gamma,alpha,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N,T_data){
    gradient<-grad_TERGM_undir_K1_Stab(theta=theta.curr, network=network, N=N, T_data=T_data)
    hess<-hess_TERGM_undir_K1_Stab(theta=theta.curr, N=N, T_data=T_data)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    ## adjacency matrix
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N,T_data=T_data)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_undir_K1_Stab(theta = theta[iter_index], network=network, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for (t in 1:(T_data-1)){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i]
            cluster_id_j<-cluster_ids_est[j]
            exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
            Stab_stat<-(network[i,j,t+1]*network[i,j,t])+((1-network[i,j,t+1])*(1-network[i,j,t]));
            t1<-t1+((Stab_stat*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
          }
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for (t in 1:(T_data-1)){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            Stab_stat<-(network[i,j,t+1]*network[i,j,t])+((1-network[i,j,t+1])*(1-network[i,j,t]));
            comp_val<-comp_val+((Stab_stat*(2*theta))-log_exp_val)
          }
        }
      }
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-K*log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K, T_data = T_data)
      t2<-log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  K<-nclust ## Defining the number of clusters
  
  #################################################
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (mixture proportions)
    start[[3]]<-theta_init ## theta
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,param[1:n_last],ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    cluster_ids_est_matrix<-matrix(rep(cluster_ids_est,T_data-1),N,T_data-1)
    
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]], theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-RI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est_matrix,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for HMM directed transitivity case to run for different number of clusters
wrapper_TERGM_dir_Trans<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,sim_indicator_HMM=0,sim_indicator_TERGM=1,theta_true=NA,pi_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(combinat)
  library(quadprog)
  library(TERGMdirTrans)
  library(Matrix)
  
  #################################################
  gamma.update.wrapper<-function(gamma.curr,alpha.curr,theta.curr,network,Trans_stat,N,K,T_data){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_TERGM_dir_Trans(gamma=gamma.curr, pi=alpha.curr, theta=theta.curr, network=network, Trans_stat=Trans_stat,N=N, K=K, T_data=T_data)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of alpha
  alpha.update<-function(gamma.curr,N,K,T_data){
    alpha.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## Normalization of alpha
    alpha.next<-alpha.next/sum(alpha.next)
    return(alpha.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,gamma,network,Trans_stat=Trans_stat,N,K,T_data){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_TERGM_dir_Trans(theta=as.vector(theta.curr), gamma=gamma, network=network, Trans_stat=Trans_stat, N=N, K=K, T_data=T_data)
    hess<-hess_TERGM_dir_Trans(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K, T_data=T_data)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    alpha<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    
    ## Assigning the starting values
    gamma[,,1]<-start[[1]]
    alpha[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    
    Trans_stat<-Trans_stat_cal(network=network,N=N,T_data=T_data)
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<100)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Updating the intial distribution alpha
      alpha[,iter_index]<-alpha.update(gamma.curr = gamma[,,iter_index-1],N,K,T_data)
      
      ## Updating gamma
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr = gamma[,,iter_index-1],alpha.curr = alpha[,iter_index],theta.curr = theta[,iter_index-1],network = network,Trans_stat=Trans_stat,N = N,K = K,T_data = T_data)
      
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1],gamma=gamma[,,iter_index], network=network, Trans_stat=Trans_stat, N=N, K=K, T_data = T_data)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_dir_Trans(gamma=gamma[,,iter_index], alpha=alpha[,iter_index], theta = theta[,iter_index], network=network, Trans_stat=Trans_stat, N=N, K=K, T_data = T_data)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    print(theta[,iter_index-1])
    return(list(gamma,alpha,theta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  theta.update_K1<-function(theta.curr,network,Trans_stat,N,T_data){
    gradient<-grad_TERGM_dir_K1_Trans(theta=theta.curr, network=network, Trans_stat=Trans_stat, N=N, T_data=T_data)
    hess<-hess_TERGM_dir_K1_Trans(theta=theta.curr, N=N, T_data=T_data)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,n_iter,thres){
    N<-dim(network)[1]
    T_data<-dim(network)[3] ## Total time points in the network for which we have data in the form of
    
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    Trans_stat<-Trans_stat_cal(network=network,N=N,T_data=T_data)
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, Trans_stat=Trans_stat, N=N, T_data=T_data)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_TERGM_dir_K1_Trans(theta = theta[iter_index], network=network, Trans_stat=Trans_stat, N=N, T_data=T_data)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    print(theta[iter_index-1])
    return(theta)
  }
  
  #########################################################################################################
  #########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  ## Defining a function to calculate the estimate of conditional log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    Trans_stat<-Trans_stat_cal(network = network,N = N,T_data = T_data)
    if(K!=1){
      t1<-0
      for (t in 2:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            cluster_id_i<-cluster_ids_est[i]
            cluster_id_j<-cluster_ids_est[j]
            exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
            t1<-t1+((Trans_stat[i,j,t-1]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
          }
        }
      }
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      comp_val<-0
      for (t in 2:T_data){
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            comp_val<-comp_val+((Trans_stat[i,j,t-1]*2*theta)-log_exp_val)
          }
        }
      }
    }
    return(comp_val)
  }
  
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,T_data,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-K*log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, T_data = T_data, cluster_ids_est = cluster_ids_est)
      t2<-log((N*(N-1)*(T_data-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  ARI<-function(cluster_ids_est, cluster_ids_true){
    RI_time<-rep(NA_real_,ncol(cluster_ids_est))
    for (k in 1:ncol(cluster_ids_est)){
      n=length(cluster_ids_est[,k])
      RI_val=0
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          RI_val=RI_val+as.numeric((cluster_ids_est[i,k]==cluster_ids_est[j,k])==(cluster_ids_true[i,k]==cluster_ids_true[j,k]))
        }
      }
      RI_mean=RI_val/(n*(n-1)/2)
      RI_time[k]<-RI_mean
    }
    RI_time_mean<-mean(RI_time)
    return(RI_time_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  T_data<-dim(sim.net)[3]
  K<-nclust ## Defining the number of clusters
  
  #################################################
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net[,,1],alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (mixture proportions)
    start[[3]]<-theta_init ## theta
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating ICL
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K,T_data = T_data)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    cluster_ids_est_matrix<-matrix(rep(cluster_ids_est,T_data),N,T_data)
    
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]], theta = param_converge[[3]],network = sim.net,N = N,K = K,T_data = T_data, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      if(sim_indicator_TERGM==1){
        RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-RI_val
      }else if (sim_indicator_HMM==1){
        ARI_val<-ARI(cluster_ids_est = cluster_ids_est_matrix,cluster_ids_true = cluster_ids_true)
        Rand_val_final<-ARI_val
      }
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,Rand_val_final)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################

