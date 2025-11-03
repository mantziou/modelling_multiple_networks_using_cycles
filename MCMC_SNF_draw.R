
#MCMC for sampling G* from prior Gm functional form, setting the Go and γο each time
#depending on the chain iteration 

log_target<-function(G,go,Go,c_G,c_Go,lam){
  return(-go*new_metric(G,Go,c_G,c_Go,lam)) 
}

MCMC_priorGm<-function(start_G,c_start_G,go,Go,c_Go,iterations,burn_in,n_nodes,prob_vec,N_data,lam,pert_centr,cycsize){
  #create chain
  chain_G<-matrix(nrow = iterations-burn_in,ncol = length(Go)) #chain for updating Gs (sample from SNF model)
  chain_cyc_G<-list()
  chain_dis<-matrix(nrow = iterations-burn_in,ncol = 1) #chain for keeping the new_metric dist for each draw of the sample with the G_fixed
  
  #initialize current values
  G_current<-start_G
  cyc_G_current<-c_start_G
  cyc_Go<-c_Go
  dis_current<-new_metric(G_current,Go,cyc_G_current,cyc_Go,lam) #get the distance btween G_fixed and G_current
  
  
  sam<-1 #counter for samples after burn_in
  
  for (i in 1:iterations){
    v<-rmultinom(1,1,c(0.4,0.4,0.5))
    if (v[1]==1){
      prop_Gv<-abs(G_current-rbinom(length(Go),1,pert_centr[1]))
      cyc_G_prop<-all_unique_cycles(graph_reform(prop_Gv,n_nodes),cycsize)
      
      mhr=exp(log_target(prop_Gv,go,Go,cyc_G_prop,cyc_Go,lam)-log_target(G_current,go,Go,cyc_G_current,cyc_Go,lam))
      prob1<-min(1, mhr)
      gen<-rbinom(n=1,size =1, prob = prob1)
      if (gen==1){
        G_current<-prop_Gv
        cyc_G_current<-cyc_G_prop
        dis_current<-new_metric(prop_Gv,Go,cyc_G_prop,cyc_Go,lam)
      }
    }
    
    if (v[2]==1){
      prop_Gv<-abs(G_current-rbinom(length(Go),1,pert_centr[2]))
      cyc_G_prop<-all_unique_cycles(graph_reform(prop_Gv,n_nodes),cycsize)
      
      mhr=exp(log_target(prop_Gv,go,Go,cyc_G_prop,cyc_Go,lam)-log_target(G_current,go,Go,cyc_G_current,cyc_Go,lam))
      prob1<-min(1, mhr)
      gen<-rbinom(n=1,size =1, prob = prob1)
      if (gen==1){
        G_current<-prop_Gv
        cyc_G_current<-cyc_G_prop
        dis_current<-new_metric(prop_Gv,Go,cyc_G_prop,cyc_Go,lam)
      }
    }
    
    if (v[3]==1){
      prop_Gv<-abs(G_current-rbinom(length(Go),1,pert_centr[3]))
      cyc_G_prop<-all_unique_cycles(graph_reform(prop_Gv,n_nodes),cycsize)
      
      mhr=exp(log_target(prop_Gv,go,Go,cyc_G_prop,cyc_Go,lam)-log_target(G_current,go,Go,cyc_G_current,cyc_Go,lam))
      prob1<-min(1, mhr)
      gen<-rbinom(n=1,size =1, prob = prob1)
      if (gen==1){
        G_current<-prop_Gv
        cyc_G_current<-cyc_G_prop
        dis_current<-new_metric(prop_Gv,Go,cyc_G_prop,cyc_Go,lam)
      }
    }
    
    # if (v[3]==1){
    #   prop_Gv<-rbinom(length(Go),1,prob_vec)
    #   cyc_G_prop<-all_unique_cycles(graph_reform_undir(prop_Gv,n_nodes))
    # 
    #   d1<-sum(dbinom(G_current,1,prob_vec,log = TRUE)) #proposal distr density for current state
    #   d2<-sum(dbinom(prop_Gv,1,prob_vec,log = TRUE)) #proposal distr density for proposal vector
    # 
    #   mhr=exp(log_target(prop_Gv,go,Go,cyc_G_prop,cyc_Go,lam)+d1-log_target(G_current,go,Go,cyc_G_current,cyc_Go,lam)-d2)
    #   prob1<-min(1, mhr)
    #   gen<-rbinom(n=1,size =1, prob = prob1)
    #   if (gen==1){
    #     G_current<-prop_Gv
    #     cyc_G_current<-cyc_G_prop
    #   }
    # }
    #Keep updates in chain if iteration>burn_in
    if(i>burn_in){
      chain_G[sam, ]<-G_current
      chain_cyc_G[[sam]]<-cyc_G_current
      chain_dis[sam, ]<-dis_current
      sam<-sam+1
    }
  }
  ind_return<-as.integer(seq(iterations-burn_in,1,-((iterations-burn_in)/N_data)))
  
  #return burnt and lagged chain
  #return(chain_dis)
  
  return(list(chain_G[ind_return, ],chain_cyc_G[ind_return],chain_dis))
}
