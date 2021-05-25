# Chapter 2.4 - examples
library(igraph)
##

# Same vertices as the graph in Chapter 2.4.1 in Lecture Notes
# Edges created randomly may differ from those in L.Notes
skthBIG <- function(showplot=FALSE){
tmp.card_alphai <- c(3,4,1,5)
tmp.card_betak <- c(1,2,3,3,2,1,1)



# Generate random bipartite graph.
# g <- sample_bipartite(6, 9, type="Gnm",m=15,directed=TRUE) # this may generate a graph with isolated nodes
g <- sample_degseq(out.deg=c(tmp.card_alphai,rep(0,length(tmp.card_betak))), in.deg = c(rep(0,length(tmp.card_alphai)),tmp.card_betak))


# Simplify graph by removing loops and multiple edges
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
              edge.attr.comb = igraph_opt("edge.attr.comb"))

# Apply bipartite layout
LO_bipart <- layout_as_bipartite(g,types=bipartite_mapping(g)$type)
LO_bipart[bipartite_mapping(g)$type==FALSE,2] <- 0
LO_bipart[bipartite_mapping(g)$type==TRUE,2] <- 1



nodecolor <- rep("yellow",length(V(g)))
nodecolor[bipartite_mapping(g)$type==TRUE] <- "orange"

# Plot BIG
if(showplot){
plot(g, vertex.label=V(g), vertex.size=10,vertex.label.dist=0,vertex.label.cex=1.25,
     vertex.color=nodecolor, layout=LO_bipart[,2:1])
}
return(list(G=g))
}


# To make total number of in and out-degrees equal
degcorrection <- function(controldeg,targetdeg){
  check=as.numeric(sum(targetdeg)!=sum(controldeg))
  while(check==1){
    if(sum(targetdeg)>sum(controldeg)){diff <- sum(targetdeg)-sum(controldeg)
    selones <- targetdeg==1
    targetdeg[!selones] <- trunc(targetdeg[!selones]*(sum(controldeg)-sum(selones))/(sum(targetdeg)-sum(selones)))
    }
    
    if(sum(targetdeg)<sum(controldeg)){diff <- sum(controldeg)-sum(targetdeg)
    sel <- sample(1:length(targetdeg),diff,replace=TRUE,prob=targetdeg)
    targetdeg[unique(sel[order(sel)])] <- targetdeg[unique(sel[order(sel)])] + as.vector(table(sel))
    }
    check=as.numeric(sum(targetdeg)!=sum(controldeg))
  }
  return(targetdeg)
}

# Generate two random graphs with deg. dist.: uniform and skewed
mainrndBIG <- function(sizeF=50,sizeOmega=10,meanoutdeg=10,showplot=FALSE){
  ## Exponential degree distribution
  set.seed(260521)
  out_degs_exp <- sample(1:sizeF, sizeF, replace=TRUE, prob=exp(-(1/meanoutdeg)*(1:sizeF)))
  
  
  meanindeg <- trunc(sum(out_degs_exp)/sizeOmega+0.5)
  set.seed(260521)
  in_degs <- sample(1:sizeOmega, sizeOmega, replace=TRUE, prob=exp(-(1/meanindeg)*(1:sizeOmega)))
  in_degs <- degcorrection(out_degs_exp,in_degs)   
  
  
  set.seed(260521)
  grnd <- sample_degseq(out.deg=c(out_degs_exp,rep(0,length(in_degs))),in.deg = c(rep(0,length(out_degs_exp)),in_degs))
  g_exp <- make_empty_graph(n=length(V(grnd)))
  g_exp <- add_edges(g_exp,edges=c(t(unique(as_edgelist(grnd)))))
  
  ## Uniform degree distribution: with same F, Omega and the total number of edges as the exponential one
  set.seed(260521)
  out_degs_uni <- trunc(runif(sizeF,1,trunc((2*sum(out_degs_exp)/sizeF-1)+0.5)))
  out_degs_uni <- degcorrection(out_degs_exp,out_degs_uni) 
  
  set.seed(260521)
  in_degs <- trunc(runif(sizeOmega,1,trunc((2*sum(out_degs_uni)/sizeOmega-1)+0.5)))
  in_degs <- degcorrection(out_degs_uni,in_degs)   
  
  set.seed(260521)
  grnd <- sample_degseq(out.deg=c(out_degs_uni,rep(0,length(in_degs))),in.deg = c(rep(0,length(out_degs_uni)),in_degs))
  g_uni <- make_empty_graph(n=length(V(grnd)))
  g_uni <- add_edges(g_uni,edges=c(t(unique(as_edgelist(grnd)))))
  
  if(showplot){par(mfrow=c(1,2))
    hist(out_degs_uni,xlab='Degree',main='Uniform',breaks=max(out_degs_uni)*5)
    hist(out_degs_exp,xlab='Degree',main='Exponential',breaks=max(out_degs_exp)*5)}
  
  return(list(Guniform=g_uni,Gskewed=g_exp))
  
}


# zi-values
zFun <- function(popgraph,coefgamma=0,n=2,multiplicity=FALSE){
  edgeik <- data.frame(as_edgelist(popgraph))
  colnames(edgeik) <- c('i','k')
  idx_F <- unique(edgeik$i)
  idx_omega <- unique(edgeik$k)
  card_alphai <- NULL
  for(i in idx_F){
    card_alphai <- c(card_alphai,sum(edgeik$i %in% i))  
  }
  
  card_betak <- NULL
  for(k in idx_omega){
    card_betak <- c(card_betak,sum(edgeik$k %in% k))  
  }
  yk <- rep(1,length(idx_omega))
  probi <- rep(n/length(idx_F),length(idx_F))
  zi <-  NULL
  for(i in idx_F){
    if(i %in% edgeik$i){
      tmp.k <- edgeik$k[edgeik$i %in% i]
      if(multiplicity){tmp.zi <- sum(yk[idx_omega %in% tmp.k]/card_betak[idx_omega %in% tmp.k])}
      if(!multiplicity){
        tmp.zi <- 0
        for(k in tmp.k){
          betak <- edgeik$i[edgeik$k %in% k]
          wik <- (probi[idx_F==i]*(1/(card_alphai[unique(edgeik$i)==i])^coefgamma)/(sum(probi[idx_F %in% betak]*(1/card_alphai[unique(edgeik$i) %in% betak])^coefgamma)))
          tmp.zi <- tmp.zi + yk[idx_omega==k]*wik
        }
      }
    }
    if(!(i %in% edgeik$i)){
      tmp.zi <- 0}
    zi <- c(zi,tmp.zi)
  }
  return(zi)
}



# HTE, HH-type estimators and priority-rule estimators
mainBIGSIWE <- function(popgraph,coefgamma=0,n=2){
  edgeik <- data.frame(as_edgelist(popgraph))
  colnames(edgeik) <- c('i','k')
  idx_F <- unique(edgeik$i)
  N <- length(idx_F)
  idx_omega <- unique(edgeik$k)
  card_alphai <- NULL
  for(i in idx_F){
    card_alphai <- c(card_alphai,sum(edgeik$i %in% i))  
  }
  
  card_betak <- NULL
  for(k in idx_omega){
    card_betak <- c(card_betak,sum(edgeik$k %in% k))  
  }
  
  yk <- rep(1,length(idx_omega))
  probi <- rep(n/N,N)
  probk <- 1-choose(N-card_betak,n)/choose(N,n)
  
  
  # All possible subsets
  all.subsets <- combn(N,n)
  B <- dim(all.subsets)[2]
  
  # Estimates over random samples
  YhatHH_alpha <- NULL
  YhatHT <- NULL
  Yhatp <- matrix(0,nrow=B,ncol=3)
  idx_F_rnd <- idx_F[order(runif(N,0,1))]
  for(b in 1:B){
    s0 <- idx_F[all.subsets[,b]]
    s1 <- unique(edgeik$k[edgeik$i %in% s0])
    zi_alpha <- zFun(popgraph,coefgamma)
    probi <- rep(n/N,N)
    YhatHH_alpha <- c(YhatHH_alpha,sum(zi_alpha[idx_F %in% s0]/probi[idx_F %in% s0])) 
    YhatHT <- c(YhatHT,sum(yk[idx_omega %in% s1]/probk[idx_omega %in% s1])) 
    # For priority-rule estimators 
    col.id <- 1
    for(orderF in c('random','ascending','descending')){
      idx_F_ordered <- idx_F_rnd
      if(orderF=="ascending"){idx_F_ordered <- idx_F[order(card_alphai)]}
      if(orderF=="descending"){idx_F_ordered <- rev(idx_F[order(card_alphai)])}
      tmp.Yhatp <- 0
      for(k in s1){
        betak <- edgeik$i[edgeik$k %in% k]
        prioriloc <- min(which(idx_F_ordered %in% betak[betak %in% s0]))
        priori <- idx_F_ordered[prioriloc]
        dik <- sum(which(idx_F_ordered %in% betak)<prioriloc)
        pik <- choose(N-1-dik,n-1)/choose(N-1,n-1)
        tmp.Yhatp <- tmp.Yhatp + yk[idx_omega %in% k]/pik/card_betak[idx_omega %in% k]/probi[idx_F %in% priori]
      }
      Yhatp[b,col.id] <- tmp.Yhatp
      col.id <- col.id + 1
    }
  }
  
  cat("gamma =",coefgamma,"\n")
  cat("idxF: ",idx_F, '\t',"alphai: ", card_alphai,"\n")
  
  resultIWE <- t(array(c(mean(YhatHT),mean(YhatHH_alpha),mean(Yhatp[,1]),mean(Yhatp[,2]),mean(Yhatp[,3]),c(var(YhatHT),var(YhatHH_alpha),var(Yhatp[,1]),var(Yhatp[,2]),var(Yhatp[,3]))*(B-1)/B),c(5,2)))
  colnames(resultIWE) <- c('YhatHTE','ZhatHH','Zhatp_rnd','Zhatp_asc','Zhatp_desc')
  rownames(resultIWE) <- c('ExpectedValue','Variance')
  print(resultIWE)
}

