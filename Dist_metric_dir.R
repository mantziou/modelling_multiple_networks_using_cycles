# version for cycles of fixed length and directed graphs
all_unique_cycles<-function(g,cyc_size){
  cycles = NULL
  for(v1 in V(g)){
    for(v2 in neighbors(g, v1, mode = "out")){# NOTE: for UNDIRECTED graphs mode="out" ignored
      cycles = c(cycles,lapply(all_simple_paths(g, v2,v1,mode = "out",cutoff = cyc_size), function(p) c(v1,p))) # NOTE: for UNDIRECTED graphs mode="out" ignored
    }
  }
  LongCycles = cycles[which(sapply(cycles, length) > 3)]
  return(LongCycles[sapply(LongCycles, min) == sapply(LongCycles, `[`, 1)]) #exclude duplicate cycles by considering only the smallest number involved in a cycle, and then keep only the cycles that start with this smallest number to avoid duplicates
}

#symmetric difference between directed cycles of two graphs (directed cycles not in common)
symm_diff_cycles<-function(c1,c2){
  comm_cylc<-intersect(c1,c2) #common cycles in two lists
  return(length(c1)+length(c2)-2*length(comm_cylc))# OR length(union())-length(intersection())
}

#Hamming distance and symmetric difference
new_metric<-function(v1,v2,c1,c2,lambda){ #v1,v2 vectorized adjacencies, c1,c2 directed cycles of the two graphs respectively 
  h<-sum(abs(v1-v2)) #hamming distance
  if(lambda!=0){
    symm<-symm_diff_cycles(c1,c2)
    return(h+lambda*symm)
  }else{
    return(h)
  }
}


