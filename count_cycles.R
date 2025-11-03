
count_cyc<-function(samp,nodes){
  num_cores<-detectCores()
  fin_list<-mclapply(samp,function(x) all_unique_cycles(graph_reform(x,nodes)),mc.cores=num_cores)
  fin_list<-lapply(fin_list, function(x) lapply(x, function(x) ifelse(is.null(x), list(), x))) # substitute NULL in list of lists (obtained for vectors with only 0 entries) with an empty list
  return(fin_list)
}

count_cyc_dir<-function(samp,nodes,cycsize){
  num_cores<-detectCores()
  fin_list<-mclapply(samp,function(x) all_unique_cycles(graph_reform(x,nodes),cycsize),mc.cores=num_cores)
  fin_list<-lapply(fin_list, function(x) lapply(x, function(x) ifelse(is.null(x), list(), x))) # substitute NULL in list of lists (obtained for vectors with only 0 entries) with an empty list
  return(fin_list)
}
