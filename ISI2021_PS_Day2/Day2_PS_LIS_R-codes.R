##
library(igraph)

# Skecth of the wolverine tracks, Becker (1991)
skthLISBecker <- function()
{
  plot(0,xaxt="n",yaxt="n",type="l",ylab="",xlab="Baseline",xlim=c(0,120),ylim=c(0,60),bty="n")
  lines(c(0,0),c(0,60)); lines(c(0,120),c(0,0)); lines(c(0,120),c(60,60)); lines(c(120,120),c(0,60))
  lines(c(0,10,12,19),c(40,40,36,42),lty=2); text(19,45,labels="k1") 
  lines(c(0,5,15,25),c(15,18,14,18),lty=2); text(20,12,labels="k2") 
  lines(c(31,34,38),c(25,30,20),lty=2); text(29,30,labels="k3") 
  lines(c(75,68,85,90),c(10,15,20,15),lty=2); text(86,23,labels="k4") 
  abline(v=c(2,42,82)); text(2,4,label="A1"); text(42,4,label="A2"); text(82,4,label="A3")
  abline(v=c(10,50,90)); text(10,1,label="B1"); text(50,1,label="B2"); text(90,1,label="B3")
  abline(v=c(35,75,115)); text(35,3,label="C1"); text(75,3,label="C2"); text(115,3,label="C3")
  abline(v=c(38,78,118)); text(38,-1,label="D1"); text(78,-1,label="D2"); text(118,-1,label="D3")

}


# Skecth of BIG representation of LIS, Becker (1991)
skthLISBeckerBIG <- function(showplot=FALSE)
{
  # Projection of wolverine tracks on the baseline; 
  # Areas without tracks assigned length of 1
  idx_F <- paste('i',1:7,sep='')
  idx_omega <- paste('k',1:4,sep='')
  edgeik <- data.frame(i=c('i1','i1','i2','i4','i6'),k=c('k1','k2','k2','k3','k4'))
  g <- graph_from_data_frame(edgeik, directed = TRUE)
  g <- add_vertices(g,3,attr=list(name=c('i3','i5','i7')))

  # Apply bipartite layout
  LO_bipart <- layout_as_bipartite(g,types=bipartite_mapping(g)$type)
  LO_bipart[bipartite_mapping(g)$type==FALSE,2] <- 0
  LO_bipart[bipartite_mapping(g)$type==TRUE,2] <- 1
  
  nodecolor <- rep("yellow",length(V(g)))
  nodecolor[bipartite_mapping(g)$type==TRUE] <- "orange"
  
  # Plot BIG
  if(showplot){
    plot(g, vertex.label=V(g)$name, vertex.size=10,vertex.label.dist=0,vertex.label.cex=1.25,
         vertex.color=nodecolor, layout=LO_bipart[,2:1]) }

return(list(G=g))  
}




# zi-values
zLISBecker <- function(graphstar,coefgamma=0,probi,multiplicity=FALSE){
  edgeik <- data.frame(as_edgelist(graphstar))
  colnames(edgeik) <- c('i','k')
  idx_F <- as_ids(V(graphstar)[bipartite.mapping(graphstar)$type==FALSE])
  idx_F <- idx_F[order(idx_F)]
  idx_omega <- as_ids(V(graphstar)[bipartite.mapping(graphstar)$type==TRUE])
  card_alphai <- NULL
  for(i in idx_F){
    card_alphai <- c(card_alphai,sum(edgeik$i %in% i))  
  }
  
  card_betak <- NULL
  for(k in idx_omega){
    card_betak <- c(card_betak,sum(edgeik$k %in% k))  
  }
  yk <- c(1,2,2,1)
  zi <-  NULL
  for(i in idx_F){
    if(i %in% edgeik$i){
      tmp.k <- edgeik$k[edgeik$i %in% i]
      if(multiplicity){tmp.zi <- sum(yk[idx_omega %in% tmp.k]/card_betak[idx_omega %in% tmp.k])}
      if(!multiplicity){
        tmp.zi <- 0
        for(k in tmp.k){
          betak <- edgeik$i[edgeik$k %in% k]
          wik <- (probi[idx_F==i]*(1/(card_alphai[idx_F==i])^coefgamma)/(sum(probi[idx_F %in% betak]*(1/card_alphai[idx_F %in% betak])^coefgamma)))
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






# HH-estimators over all draws
mainLISBecker <- function(graphstar,coefgamma=0,probi,multiplicity=FALSE,showcat=TRUE){
  edgeik <- data.frame(as_edgelist(graphstar))
  colnames(edgeik) <- c('i','k')
  idx_F <- as_ids(V(graphstar)[bipartite.mapping(graphstar)$type==FALSE])
  idx_F <- idx_F[order(idx_F)]
  idx_omega <- as_ids(V(graphstar)[bipartite.mapping(graphstar)$type==TRUE])

all.subsets <- list(c(1,5,6),c(1,5,6),c(4,6,7),c(4,6,7))
B <- length(all.subsets)

# Estimates over random samples
  YhatHH_alpha <- NULL
  for(b in 1:B){
    s0 <- idx_F[all.subsets[[b]]]
    s1 <- unique(edgeik$k[edgeik$i %in% s0])
    zi_alpha <- zLISBecker(graphstar,coefgamma,probi,multiplicity = multiplicity)
    YhatHH_alpha <- c(YhatHH_alpha,sum(zi_alpha[idx_F %in% s0]/probi[idx_F %in% s0])) 
  }

if(showcat){
  cat('g:',coefgamma,'\t','zi:',zi_alpha,"\n")
  
  cat("Estimate per draw: ",YhatHH_alpha, "\n")
  cat("(Estimate over all draws,VarEst): ",mean(YhatHH_alpha), '\t',var(YhatHH_alpha)/B, "\n")
}
  return(list(varest=var(YhatHH_alpha)/B))
}  

