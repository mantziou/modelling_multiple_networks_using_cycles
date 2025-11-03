log_IS_est<-function(is_sample,cyc_is_sample,centroid_cer,alpha,centroid_snf,centroid_snf_cyc,gamma,lambda){
  num_cores<-detectCores()
  ham_vec<-sapply(is_sample,function(x) Hamm_dist(x,centroid_cer)) #calculate vector of Hamm dist of IS sample from centroid (did not make it unlist(mclapply()) as already quick and mclapply() makes it slower)
  new_m_vec<-unlist(mclapply(1:length(is_sample),function(i) new_metric(is_sample[[i]],centroid_snf,cyc_is_sample[[i]],centroid_snf_cyc,lambda),mc.cores = num_cores)) #calculate new metr of IS sample from centroid (updated or not)
  denom_vec<-1/((alpha^ham_vec)*(1-alpha)^(length(centroid_cer)-ham_vec)) #denominator of Z()est
  exp_power<-(-gamma*new_m_vec)
  max_exp_power<-max(exp_power)# to avoid overflow in the case of min()
  nomin_vec<-exp(exp_power-max_exp_power) #nominator of Z()est
  #print(c(nomin_vec,"and",denom_vec))
  tot<-nomin_vec%*%denom_vec #calculate the sum as the inner product of the two vectors
  return(max_exp_power+log(tot)-log(length(is_sample)))
}


