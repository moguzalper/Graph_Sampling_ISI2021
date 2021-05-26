## Load R-package igraph
library(igraph) # If not installed previously, use install.packages("igraph")
##


# Examples similar to those in Chapter 2.4
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
  colnames(resultIWE) <- c('YhatHTE','ZhatHH','Yhatp_rnd','Yhatp_asc','Yhatp_desc')
  rownames(resultIWE) <- c('ExpectedValue','Variance')
  print(resultIWE)
}


# Generate a population graph
set.seed(270521)
popg <- skthBIG(showplot=TRUE)$G
V(popg)
E(popg)
as_edgelist(popg)
cat('|Beta_kappa| =',table(as_edgelist(popg)[,2]),'\n')
cat('|alpha_i| =', table(as_edgelist(popg)[,1]),'\n')

# zFun(popgraph,coefgamma=0,n=2,multiplicity=FALSE)
coefg <- 1
cat('g:',coefg,'\t','zi_alpha:',zFun(popg,coefg),'; var(zi_alpha):',var(zFun(popg,coefg)),"\n",'\t','zi_beta:',zFun(popg,multiplicity = TRUE),'; var(zi_beta):',var(zFun(popg,multiplicity = TRUE)),"\n")


# Variance of zi_alpha for different choices of gamma
max.gamma <- 25
range.gamma <- seq(0,max.gamma,by=0.1)
varzi_alpha <- NULL
for(tmp.gamma in range.gamma){tmp <- var(zFun(popg,tmp.gamma))
varzi_alpha <- c(varzi_alpha,tmp)}
# Gamma value which gives minimum variance
gamma.minvar <- range.gamma[which(varzi_alpha==min(varzi_alpha))]
print(gamma.minvar)
par(mfrow=c(1,1),xpd=TRUE)
plot(range.gamma,varzi_alpha,xlab='gamma',ylab='V(zi_alpha)')
abline(h=var(zFun(popg,0,multiplicity=TRUE)))
text(max.gamma*0.95,var(zFun(popg,0,multiplicity=TRUE))*0.98,label='V(zi_beta)',pos=1)
abline(v=gamma.minvar,lty=2)
text(gamma.minvar,max(varzi_alpha)*1.05,label=paste('gamma.minvar=',gamma.minvar),adj=1)


# Minimum variance
var(zFun(popg,gamma.minvar))
# Variance with multiplicity weights
var(zFun(popg,0))

mainBIGSIWE(popg)

mainBIGSIWE(popg,coefgamma=gamma.minvar)

mainBIGSIWE(popg,coefgamma=0.5)


