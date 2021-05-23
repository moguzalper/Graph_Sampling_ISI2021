# -----------------------------------------------------------------
# Subsidiary functions
# -----------------------------------------------------------------
checkSum <- function(x,vrb)
{
  sum(x %in% vrb)
}

# Identify # of cliques observed via any of its incident nodes in s0
cliqueFun <- function(graph,v,ref.nodes)
{
  a <- cliques(graph,min=v,max=v)
  if(length(a)!=0) {temp <- as.data.frame(t(apply(sapply(a, as_ids),2,sort)))
  cnt <- sum(apply(temp,1,checkSum,vrb=as_ids(ref.nodes))>0)
  }
  if(length(a)==0) {cnt <- 0}
  return(cnt)
}


# Identify # of motifs observed via any of its incident nodes in s0
motifFun <- function(pattern,graph,ref.nodes)
{
  a <- subgraph_isomorphisms(pattern,graph,induced=T)
  
  if(length(a)!=0) {temp <- as.data.frame(t(apply(sapply(a, as_ids),2,sort)))
  temp.red <- unique(temp)
  cnt <- sum(apply(temp.red,1,checkSum,vrb=as_ids(ref.nodes))>0)
  }
  if(length(a)==0) {cnt <- 0}
  return(cnt)
}

# Variance of the HTE of a clique under SRS with induced observation procedure
varSRSCliqueFun <- function(graph,v,temp.n)
{
  temp.N <- length(V(graph))
  a <- cliques(graph,min=v,max=v)
  temp <- as.data.frame(t(apply(sapply(a, as_ids),2,sort)))
  prk <- choose(temp.N-v,temp.n-v)/choose(temp.N,temp.n)
  var.clique <- 0
  for(i in 1:nrow(temp)){
    for(j in 1:nrow(temp)){
      temp.m <- length(unique(unlist(c(temp[i,],temp[j,]))))
      prkl <- choose(temp.N-temp.m,temp.n-temp.m)/choose(temp.N,temp.n)
      temp.sum <- (prkl-prk^2)/prk/prk
      var.clique <- var.clique + temp.sum    
    }
  }
  return(var.clique)
}

# Variance of the HTE of a motif under SRS with induced observation procedure
varSRSMotifFun <- function(pattern,graph,temp.n)
{
  temp.N <- length(V(graph))
  v <- length(V(pattern))
  a <- subgraph_isomorphisms(pattern,graph,induced=T)
  temp.all <- as.data.frame(t(apply(sapply(a, as_ids),2,sort)))
  temp <- unique(temp.all)
  prk <- choose(temp.N-v,temp.n-v)/choose(temp.N,temp.n)
  
  var.motif <- 0
  for(i in 1:nrow(temp)){
    for(j in 1:nrow(temp)){
      temp.m <- length(unique(unlist(c(temp[i,],temp[j,]))))
      prkl <- choose(temp.N-temp.m,temp.n-temp.m)/choose(temp.N,temp.n)
      temp.sum <- (prkl-prk^2)/prk/prk
      var.motif <- var.motif + temp.sum    
    }
  }
  return(var.motif)
}


# Identify # of cliques observed via any of its ancestors selected in st
# Calculate zi-values based on equal-share weighting for all i in s0
cliquezimkFun <- function(graph,v,ref.nodes,ref.node.s0,order){ 
  n.s0 <- length(ref.node.s0)
  zi <- rep(0,n.s0); mk <- 0
  a <- cliques(graph,min=v,max=v)
  if(length(a)!=0){ 
    select <- lapply(lapply(a,as_ids),checkSum,vrb=as_ids(ref.nodes))>0
    if(sum(select)>0){
      a.filtered <- a[select]
      ego.k <- lapply(a.filtered,ego,graph=graph,order=order)
      ego.k.ids <- lapply(ego.k,function(x) lapply(x,as_ids))
      ego.k.ids.unique <- lapply(lapply(ego.k.ids,unlist),unique)
      
      if(v==1){ff <- lapply(ego.k.ids.unique,function(x){x[[1]]})
      ego.k.ids.unique <- ego.k.ids.unique[order(as.numeric(ff))]
      }
      
      mk <- unlist(lapply(ego.k.ids.unique,length))
      
      zi <- NULL
      for(i in 1:n.s0){
        temp.ref.node <- ref.node.s0[i]
        temp <- lapply(ego.k.ids.unique,checkSum,vrb=as_ids(temp.ref.node))
        temp.zi <- sum(unlist(temp)/mk)
        zi <- c(zi,temp.zi)
      }
    }
  }
  return(list(zi,mk))
}


# Identify # of motifs observed via any of its ancestors selected in st
# Calculate zi-values based on equal-share weighting for all i in s0
motifzimkFun <- function(pattern,graph,ref.nodes,ref.node.s0,order){ 
  n.s0 <- length(ref.node.s0)
  zi <- rep(0,n.s0); mk <- 0
  a <- subgraph_isomorphisms(pattern,graph,induced=T)
  if(length(a)!=0){  a.sorted <- lapply(a,sort)
  select <- lapply(lapply(a.sorted,as_ids),checkSum,vrb=as_ids(ref.nodes))>0
  if(sum(select)>0){
    a.filtered <- a.sorted[select]
    temp <- as.data.frame(t(sapply(a.filtered,as_ids)))
    dup <- duplicated(temp)
    a.filtered <- a.filtered[!dup]
    ego.k <- lapply(a.filtered,ego,graph=graph,order=order)
    ego.k.ids <- lapply(ego.k,function(x) lapply(x,as_ids))
    ego.k.ids.unique <- lapply(lapply(ego.k.ids,unlist),unique)
    mk <- unlist(lapply(ego.k.ids.unique,length))
    
    zi <- NULL
    for(i in 1:n.s0){
      temp.ref.node <- ref.node.s0[i]
      temp <- lapply(ego.k.ids.unique,checkSum,vrb=as_ids(temp.ref.node))
      temp.zi <- sum(unlist(temp)/mk)
      zi <- c(zi,temp.zi)
    } 
  }
  }
  return(list(zi,mk))
}



# Identify # of cliques observed via any of its ancestors selected in st
# Calculate zi-values based on PIDA weights for all i in s0
cliqueziUnSFun <- function(graph,ref.graph,v,ref.nodes,ref.node.s0,order){ 
  n.s0 <- length(ref.node.s0)
  zi <- rep(0,n.s0)
  a <- cliques(ref.graph,min=v,max=v)
  b <- cliques(graph,min=v,max=v)
  if(length(a)!=0){ 
    select <- lapply(lapply(a,as_ids),checkSum,vrb=as_ids(ref.nodes))>0
    if(sum(select)>0){
      a.filtered <- a[select]
      a.filtered
      
      ego.k <- lapply(a.filtered,ego,graph=ref.graph,order=order)
      ego.k.ids <- lapply(ego.k,function(x) lapply(x,as_ids))
      ego.k.ids.unique <- lapply(lapply(ego.k.ids,unlist),unique)
      
      if(v==1){ff <- lapply(ego.k.ids.unique,function(x){x[[1]]})
      ego.k.ids.unique <- ego.k.ids.unique[order(as.numeric(ff))]
      }
      
      alphai.k <- rep(list(),length(ego.k.ids.unique))
      for(i in 1:length(ego.k.ids.unique)){ 
        alphai <- NULL
        temp.ego.k.ids.unique <- ego.k.ids.unique[[i]]
        
        for(j in 1:length(temp.ego.k.ids.unique)){ 
          
          temp.node <- temp.ego.k.ids.unique[j]
          
          temp.vrb <- unlist(lapply(ego(graph,order=order,nodes=V(graph)[as_ids(V(graph))==temp.node]),as_ids))
          
          temp.alphai <- sum(lapply(lapply(b,as_ids),checkSum,vrb=temp.vrb)>0)
          alphai <- c(alphai,temp.alphai)
        }
        alphai.k[[i]] <- alphai
      }
      
      Pik <- lapply(alphai.k,function(x){(1/x)/(sum(1/x))})
      
      zi <- NULL
      for(i in 1:n.s0){
        temp.ref.node <- ref.node.s0[i]
        temp <- lapply(ego.k.ids.unique,function(x){as.numeric(x %in% as_ids(temp.ref.node))})
        temp.zi <- sum(unlist(Pik)*unlist(temp))
        zi <- c(zi,temp.zi)
      } 
    }
  }
  return.result <- zi
  if(v==1) {return.result <- list(zi,Pik)}
  return(return.result)
}

# Identify # of motifs observed via any of its ancestors selected in st
# Calculate zi-values based on PIDA weights for all i in s0
motifziUnSFun <- function(pattern,graph,ref.graph,ref.nodes,ref.node.s0,order){ 
  n.s0 <- length(ref.node.s0)
  zi <- rep(0,n.s0)
  a <- subgraph_isomorphisms(pattern,ref.graph,induced=T)
  b <- subgraph_isomorphisms(pattern,graph,induced=T)
  b.sorted <- lapply(b,sort)
  b.df <-  as.data.frame(t(sapply(b.sorted,as_ids)))
  b.dup <- duplicated(b.df)
  b.unique <- b.sorted[!b.dup]
  
  if(length(a)!=0){a.sorted <- lapply(a,sort)
  select <- lapply(lapply(a.sorted,as_ids),checkSum,vrb=as_ids(ref.nodes))>0 
  if(sum(select)>0){
    a.filtered <- a.sorted[select]
    temp <- as.data.frame(t(sapply(a.filtered,as_ids)))
    dup <- duplicated(temp)
    a.filtered <- a.filtered[!dup]
    
    ego.k <- lapply(a.filtered,ego,graph=ref.graph,order=order)
    ego.k.ids <- lapply(ego.k,function(x) lapply(x,as_ids))
    ego.k.ids.unique <- lapply(lapply(ego.k.ids,unlist),unique)
    
    alphai.k <- rep(list(),length(ego.k.ids.unique))
    for(i in 1:length(ego.k.ids.unique)){ 
      alphai <- NULL
      temp.ego.k.ids.unique <- ego.k.ids.unique[[i]]
      
      for(j in 1:length(temp.ego.k.ids.unique)){ 
        temp.node <- temp.ego.k.ids.unique[j]
        temp.vrb <- unlist(lapply(ego(graph,order=order,nodes=V(graph)[as_ids(V(graph))==temp.node]),as_ids))
        temp.alphai <- sum(lapply(lapply(b,as_ids),checkSum,vrb=temp.vrb)>0)
        alphai <- c(alphai,temp.alphai)
      }
      alphai.k[[i]] <- alphai
    }
    
    Pik <- lapply(alphai.k,function(x){(1/x)/(sum(1/x))})
    
    
    zi <- NULL
    for(i in 1:n.s0){
      temp.ref.node <- ref.node.s0[i]
      temp <- lapply(ego.k.ids.unique,function(x){as.numeric(x %in% as_ids(temp.ref.node))})
      temp.zi <- sum(unlist(Pik)*unlist(temp))
      zi <- c(zi,temp.zi)
    } 
  }
  }
  return(zi)
}

