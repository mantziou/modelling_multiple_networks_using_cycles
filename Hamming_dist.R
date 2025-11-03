#function for hamming distance btwn two adjacencies
Hamm_dist<-function(v1,v2){ 
  return(sum(abs(v1-v2))) 
}
