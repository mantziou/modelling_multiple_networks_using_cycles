####################################
###### DIRECTED NETWORKS ######
####################################

#function that vectorizes adjacency matrices, excluding diagonal for directed networks
mat_vect<-function(mat){
  return(c(mat[which(upper.tri(mat),arr.ind = TRUE)],mat[which(lower.tri(mat),arr.ind = TRUE)])) 
}

mat_reform<-function(vec,n_nodes){ 
  vec_len<-length(vec)
  t<-matrix(rep(0,n_nodes*n_nodes),nrow = n_nodes,ncol = n_nodes)
  t[upper.tri(t)]<-vec[1:(vec_len/2)]
  t[lower.tri(t)]<-vec[((vec_len/2)+1):vec_len]
  return(t)
}

graph_reform<-function(vec,n_nodes){ 
  vec_len<-length(vec)
  t<-matrix(rep(0,n_nodes*n_nodes),nrow = n_nodes,ncol = n_nodes)
  t[upper.tri(t)]<-vec[1:(vec_len/2)]
  t[lower.tri(t)]<-vec[((vec_len/2)+1):vec_len]
  g<-graph_from_adjacency_matrix(t,mode = "directed",weighted = NULL)
  return(g)
}

####################################
###### UNDIRECTED NETWORKS ######
####################################

graph_reform_undir<-function(vec,n_nodes){ 
  t<-matrix(rep(0,n_nodes*n_nodes),nrow = n_nodes,ncol = n_nodes)
  t[upper.tri(t)]<-vec
  t<-t(t)+t
  g<-graph_from_adjacency_matrix(t,mode = "undirected")
  return(g)
}

#function that vectorizes adjacency matrices, excluding diagonal for UNdirected networks
mat_vect_undir<-function(mat){
  return(mat[which(upper.tri(mat),arr.ind = TRUE)])
}

#get adjacency from vector (Undirected case)
adjac_reform_undir<-function(vec,num_nodes){ 
  t<-matrix(rep(0,num_nodes*num_nodes),nrow = num_nodes,ncol = num_nodes)
  t[upper.tri(t)]<-vec
  t<-t(t)+t
  return(t)
}