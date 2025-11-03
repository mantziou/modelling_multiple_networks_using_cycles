library(stats)
library(boot)
library(igraph)
#library(sand)
library(coda)
library(nnet)
library(Matrix)
library(parallel)

setwd("YOUR_DIRECTORY")
source("Log_post_SNF.R")
source("IS_est_3.R")
#source("MCMC_CER_draw.R")
source("count_cycles.R")
source("Dist_metric_dir.R")
setwd(paste(dir_path,"Sally_project",sep = ""))
source("mat_to_vec.R")
source("Hamming_dist.R")

MH_SNF_ham_symdif<-function(G_data_list,N_data,n_nodes,start_centr_v,start_graph_cer,start_gamm,iterations,burn_in,prob_vec,alpha_cer,centroid_cer,gamma_0,centr_0,alpha_0,beta_0,lam,iter_is,burn_is,K_is,bound,pert_centr,pert_gamma,cyc_size){
  #create chains
  chain_centr<-matrix(nrow=iterations-burn_in,ncol = length(G_data_list[[1]]))
  chain_cyc_centr<-list()
  chain_gamma<-matrix(nrow = iterations-burn_in,ncol = 1)  
  
  #initialize current values
  centr_current<-start_centr_v
  gamma_current<-start_gamm
  
  cyc_centr_current<-all_unique_cycles(graph_reform(centr_current,n_nodes),cyc_size) #all cycles for current centroid
  cyc_centr_0<-all_unique_cycles(graph_reform(centr_0,n_nodes),cyc_size) #all cycles for hyperparameter centroid (fixed throughout simulation)
  cyc_G_data<-count_cyc_dir(G_data_list,n_nodes,cyc_size)#list containing at each entry the directed cycles for each network data (fixed throughout simulation)
  
  sam<-1 #counter for samples after burn_in
  
  for(i in 1:iterations){
    print(i)
    v<-rmultinom(1,1,c(0.3,0.3,0.3,0.2,0.5,0.4,0.3,0.2)) # rep(1/6, 6), 6 categories for the different proposals for the parameters
    #print(v)
    if (v[1]==1){ # 1st case: Updating vector for adjac
      #proposal for Gm centroid
      prop_centr<-abs(centr_current-rbinom(length(G_data_list[[1]]),1,pert_centr[1])) 
      cyc_centr_prop<-all_unique_cycles(graph_reform(prop_centr,n_nodes),cyc_size)
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,prop_centr,cyc_centr_prop,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam))
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        centr_current<-prop_centr # Updating ONLY the vec of adj, gamma stays the same
        cyc_centr_current<-cyc_centr_prop 
      }
    }
    
    if (v[2]==1){ # 1st case: Updating vector for adjac
      #proposal for Gm centroid
      prop_centr<-abs(centr_current-rbinom(length(G_data_list[[1]]),1,pert_centr[2])) 
      cyc_centr_prop<-all_unique_cycles(graph_reform(prop_centr,n_nodes),cyc_size)
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,prop_centr,cyc_centr_prop,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam))
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        centr_current<-prop_centr # Updating ONLY the vec of adj, gamma stays the same
        cyc_centr_current<-cyc_centr_prop 
      }
    }
    
    if (v[3]==1){ # 1st case: Updating vector for adjac
      #proposal for Gm centroid
      prop_centr<-abs(centr_current-rbinom(length(G_data_list[[1]]),1,pert_centr[3])) 
      cyc_centr_prop<-all_unique_cycles(graph_reform(prop_centr,n_nodes),cyc_size)
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,prop_centr,cyc_centr_prop,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam))
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        centr_current<-prop_centr # Updating ONLY the vec of adj, gamma stays the same
        cyc_centr_current<-cyc_centr_prop 
      }
    }
    
    if(v[4]==1){ # 3rd case: Updating vec of adj by Bernoulli with diff probs for each edge
      prop_centr<-rbinom(length(G_data_list[[1]]),1,prob_vec) # generate edges of proposal with diff probs according to prob_vec
      cyc_centr_prop<-all_unique_cycles(graph_reform(prop_centr,n_nodes),cyc_size)
      
      d1<-sum(dbinom(centr_current,1,prob_vec,log = TRUE)) #proposal distr density for current state
      d2<-sum(dbinom(prop_centr,1,prob_vec,log = TRUE)) #proposal distr density for proposal vector
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,prop_centr,cyc_centr_prop,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)+d1-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-d2)
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        centr_current<-prop_centr # Updating ONLY the vec of adj, gamma stays the same
        cyc_centr_current<-cyc_centr_prop 
      }
    }
    
    if(v[5]==1){ # 5th case: Updating gammas with Uniform proposal 
      y<-gamma_current+runif(1,min = -pert_gamma[1],max = pert_gamma[1]) # candidate proposal y, check the below constraints
      if (y<0){
        y<-(-y)
      }
      if (y>bound){ # !look into bound of gamma
        y<-2*bound-y
      }
      prop_gamm<-y
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,prop_gamm,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam))
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        gamma_current<-prop_gamm # Updating ONLY gamma
      }
    }
    
    if(v[6]==1){ # 5th case: Updating gammas with Uniform proposal 
      y<-gamma_current+runif(1,min = -pert_gamma[2],max = pert_gamma[2]) # candidate proposal y, check the below constraints
      if (y<0){
        y<-(-y)
      }
      if (y>bound){ # !look into bound of gamma
        y<-2*bound-y
      }
      prop_gamm<-y
      
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,prop_gamm,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam))
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        gamma_current<-prop_gamm # Updating ONLY gamma
      }
    }
    
    if(v[7]==1){ # 5th case: Updating gammas with Uniform proposal 
      y<-gamma_current+runif(1,min = -pert_gamma[3],max = pert_gamma[3]) # candidate proposal y, check the below constraints
      if (y<0){
        y<-(-y)
      }
      if (y>bound){ # !look into bound of gamma
        y<-2*bound-y
      }
      prop_gamm<-y
      
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,prop_gamm,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam))
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        gamma_current<-prop_gamm # Updating ONLY gamma
      }
    }
    
    if(v[8]==1){ # 5th case: Updating gammas with Uniform proposal
      y<-gamma_current+runif(1,min = -pert_gamma[4],max = pert_gamma[4]) # candidate proposal y, check the below constraints
      if (y<0){
        y<-(-y)
      }
      if (y>bound){ # !look into bound of gamma
        y<-2*bound-y
      }
      prop_gamm<-y
      
      
      # #draw IS sample from CER with parameters centroid_cer, alpha_cer
      # mh_cer<-MCMC_CER_G_draw(start_graph_cer,alpha_cer,centroid_cer,iter_is,burn_is,K_is,pert_centr)
      # is_sample<-lapply(1:K_is,function(j) mh_cer[j, ]) #make the matrix (chain) a list of vectors
      
      # draw IS sample from CER without MCMC
      is_sample <- replicate(K_is,abs(centroid_cer-rbinom(length(centroid_cer),1,alpha_cer)),simplify = FALSE)
      is_sample_cyc<-count_cyc_dir(is_sample,n_nodes,cyc_size)
      
      mhr=exp(LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,prop_gamm,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)-
                LogPost(N_data,G_data_list,cyc_G_data,is_sample,is_sample_cyc,centr_current,cyc_centr_current,gamma_current,centroid_cer,alpha_cer,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam))
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if(gen==1){
        gamma_current<-prop_gamm # Updating ONLY gamma
      }
    }
    
    #Keep updates in chain if iteration>burn_in
    if(i>burn_in){
      chain_centr[sam, ]<-centr_current
      chain_cyc_centr[[sam]]<-cyc_centr_current
      chain_gamma[sam, ]<-gamma_current
      sam<-sam+1
    }
    
  }
  return(list(chain_centr,chain_cyc_centr,chain_gamma))
}
