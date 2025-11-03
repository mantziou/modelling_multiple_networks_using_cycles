
# Log Likelihood of data
LogLikel<-function(N_data,G_data_list,cyc_G_data_list,is_sample,cyc_is_sample,centr_snf,cyc_centr_snf,gamm,centr_cer,alpha,lam){
  h<-0
  for(i in 1:N_data){
    h<-h+new_metric(G_data_list[[i]],centr_snf,cyc_G_data_list[[i]],cyc_centr_snf,lam)
  }
  log_z_est<-log_IS_est(is_sample,cyc_is_sample,centr_cer,alpha,centr_snf,cyc_centr_snf,gamm,lam)
  Likel<-(-(gamm*h)-(N_data*log_z_est))
  return(Likel)
}


#log-priors for Gm and γ
LogPriors<-function(centr_snf,cyc_centr_snf,gamm,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam){ #centr should be a vector
  centr_prior<-(-gamma_0*new_metric(centr_snf,centr_0,cyc_centr_snf,cyc_centr_0,lam)) 
  gamm_prior<- ((alpha_0-1)*log(gamm))-(beta_0*gamm) # Gamma prior
  return(sum(centr_prior,gamm_prior))
}


#Log posterior for Gm, γ, G*
LogPost<-function(N_data,G_data_list,cyc_G_data_list,is_sample,cyc_is_sample,centr_snf,cyc_centr_snf,gamm,centr_cer,alpha,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam){
  return(sum(LogLikel(N_data,G_data_list,cyc_G_data_list,is_sample,cyc_is_sample,centr_snf,cyc_centr_snf,gamm,centr_cer,alpha,lam),LogPriors(centr_snf,cyc_centr_snf,gamm,gamma_0,centr_0,cyc_centr_0,alpha_0,beta_0,lam)))
}