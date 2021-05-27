#
## Load R-packages igraph and latex2exp
library(igraph) # If not installed previously, use install.packages("igraph") and install.packages("latex2exp"), respectively
library(latex2exp)
##



# -----------------------------------------------------------------
# Create a random population graph
# -----------------------------------------------------------------
skthG <- function(sizeF=40,p=0.1,showplot=FALSE){
  set.seed(051119)
  g <- erdos.renyi.game(n=sizeF,p,type='gnp',directed=F)
  vertex_attr(g) <- list(name = seq(1:length(V(g))))
  if(showplot){plot(g)}
  cat('Diameter: ',diameter(g,directed=FALSE),'\n')
  return(list(G=g))
}


# Generate a subgraph with a chosen order and type 
makeMotif <- function(orderM=2,motif=c('clique','cycle','star','path')[1]){
  if(motif=='clique'){subg <- make_full_graph(orderM,directed=FALSE)}
  if(motif=='cycle'){subg <- make_ring(orderM)}
  if(motif=='star'){subg <- make_star(orderM,mode='undirected')}
  if(motif=='path'){E <- NULL; for(a in 1:(orderM-1)){E <- c(E,a,(a+1))}
  subg <- graph(edges=E,n=orderM,directed=FALSE)}
  return(subg)
}



# Plot motifs up to order 4
skthMotifs <- function(){
  # Make subgraphs
  dyad.g <- makeMotif(2,'clique')
  star2.g <- makeMotif(3,'star')
  triangle.g <- makeMotif(3,'clique')
  cycle4.g <- makeMotif(4,'cycle')
  star3.g <- makeMotif(4,'star')
  clique4.g <- makeMotif(4,'clique')
  path3.g  <-  makeMotif(4,'path')
  
  # Plot subgraphs
  par(mfrow=c(3,4))
  par(cex.lab=1.5)
  plot(dyad.g,vertex.label=NA,xlab=TeX('2-clique (dyad,$K_2$)'))
  plot(star2.g,vertex.label=NA,xlab=TeX('2-star ($S_2$)'))
  plot(triangle.g,vertex.label=NA,xlab=TeX('3-clique (triangle,$K_3$)'))
  plot(clique4.g,vertex.label=NA,xlab=TeX('4-clique ($K_4$)'))
  plot(cycle4.g,vertex.label=NA,xlab=TeX('4-cycle ($C_4$)'))
  plot(star3.g,vertex.label=NA,xlab=TeX('3-star ($S_3$)'))
  plot(path3.g,vertex.label=NA,xlab=TeX('3-path ($P_3$)'))
  
  par(mfrow=c(1,1))
  par(cex.lab=1)
}



# Count motifs in the population graph
countMotif <- function(popgraph,orderM=2,motif=c('clique','cycle','star','path')[1]){
  if(orderM==2 | (orderM==3 & motif=='cycle') | motif=='clique'){cnt <- cliqueFun(popgraph,orderM,V(popgraph))} else {
    tmp.g <- makeMotif(orderM,motif)
    cnt <- motifFun(tmp.g,popgraph,V(popgraph))
  }
  return(cnt)
}


# Variance of the HTE of a motif under SRS with induced observation procedure
varSRSinduced <- function(popgraph,orderM=2,motif=c('clique','cycle','star','path')[1],sizes0=2){
  if(orderM==2 | (orderM==3 & motif=='cycle') | motif=='clique'){tmp.var <- varSRSCliqueFun(popgraph,orderM,sizes0)} else {
    tmp.g <- makeMotif(orderM,motif)
    tmp.var <-  varSRSMotifFun(tmp.g,popgraph,sizes0)
  }
  return(tmp.var)
}


# Diagnostics of a motif: its diameter and observation diameter
diagMotif <- function(orderM=2,motif=c('clique','cycle','star','path')[1]){
  tmp.g <- makeMotif(orderM,motif=motif)
  diam.g <- diameter(tmp.g,directed=FALSE)
  dist.g <- distances(tmp.g)
  sel <- lower.tri(dist.g)
  argmax.dist <- sum(dist.g[sel]==max(dist.g[sel])) 
  if(argmax.dist>1){obsdiam.g <- diam.g +1 } else {
    obsdiam.g <- diam.g
  }
  return(list(diamkappa=diam.g,obsdiamkappa=obsdiam.g))    
} 


# BIGS-IWE for T-wave snowball sampling
mainTSBS <- function(popgraph,n=2,Twave=1,orderM=2,motif=c('clique','cycle','star','path')[1],B=50,includePIDA=FALSE){
  # Diagnostics of the motif: its diameter and observation diameter
  tmp.subgraph <- makeMotif(orderM,motif)
  diagmotif <- diagMotif(orderM,motif)
  
  
  # Error message if Twave<obsdiamkappa
  if(diagmotif$obsdiamkappa>Twave){stop('Observation diameter of the motif is greater than the number of waves','\n',
                                        'Wave must be equal or greater than ',diagmotif$obsdiamkappa,'\n')}
  
  goPIDA <- (diagmotif$obsdiamkappa+diagmotif$diamkappa)<=Twave # For ancestor network = M(kappa)
  gomaxgeo <- trunc((Twave-diagmotif$diamkappa)/2)
  gomaxgeoPIDA <- trunc((Twave-diagmotif$diamkappa-diagmotif$obsdiamkappa)/3)
  
  N <- length(V(popgraph))
  probi <- n/N
  all.subsets <- combn(N,n)
  nr.all.samples <- dim(all.subsets)[2]
  
  if(nr.all.samples<=1000){B <- nr.all.samples}
  
  
  # -----------------
  # SIMULATION
  # -----------------
  thetaHT <- array(0,dim=c(B,1+gomaxgeo))
  thetaHHmult <- array(0,dim=c(B,1+gomaxgeo))
  if(goPIDA & includePIDA){thetaHHPIDA  <- array(0,dim=c(B,1+gomaxgeoPIDA))}
  expn <- rep(0,gomaxgeo)
  
  for(b in 1:B){
    if(B<=1000){s0 <- all.subsets[,b]} else {set.seed(b+270521)
      s0 <- sample(N,n,replace=F)
    }
    
    # Sample graphs
    s0.g <- induced_subgraph(popgraph,V(popgraph)[s0],impl='copy_and_delete')
    st.g <- rep(list(NULL),Twave)
    tmp.samplenodes <- as_ids(V(s0.g)) 
    for(t in 1:Twave){
      tmp.g <- subgraph.edges(popgraph,eids=E(popgraph)[from(V(popgraph)[V(popgraph) %in% tmp.samplenodes])])
      st.g[[t]] <- tmp.g
      tmp.samplenodes <- as_ids(V(tmp.g))
    }
    
    
    
    # Ancestor network = M(kappa)
    sref.g <- s0.g
    s.g <- st.g[[diagmotif$obsdiamkappa]]
    
    # For HTE
    if(orderM==2 | (orderM==3 & motif=='cycle') | motif=='clique'){cnt.motif <- cliqueFun(s.g,orderM,V(sref.g))} else {
      cnt.motif <- motifFun(tmp.subgraph,s.g,V(sref.g))
    }
    
    thetaHT[b,1] <- sum(cnt.motif/(1-choose(N-orderM,n)/choose(N,n)))
    
    
    # zi-values for multiplicity estimator
    if(orderM==2 | (orderM==3 & motif=='cycle') | motif=='clique'){motif.zimk <- cliquezimkFun(s.g,orderM,V(sref.g),V(s0.g),0)} else {
      motif.zimk <- motifzimkFun(tmp.subgraph,s.g,V(sref.g),V(s0.g),0)
    }
    
    thetaHHmult[b,1] <- sum(motif.zimk[[1]]/probi)
    
    
    # zi-values based on PIDA weights
    if(goPIDA & includePIDA){
      sPIDA.g <- st.g[[diagmotif$obsdiamkappa+diagmotif$diamkappa]]
      if(orderM==2 | (orderM==3 & motif=='cycle') | motif=='clique'){motif.zimk.PIDA <- cliqueziUnSFun(sPIDA.g,s.g,orderM,V(sref.g),V(s0.g),0)} else {
        motif.zimk.PIDA <- motifziUnSFun(tmp.subgraph,sPIDA.g,s.g,V(sref.g),V(s0.g),0)
      }
      thetaHHPIDA[b,1] <- sum(motif.zimk.PIDA/probi)
    }
    
    
    # Extended ancestor network = M(kappa) U beta^t_kappa(M); HTE & multiplicity estimator
    if(gomaxgeo>0){for(geodist in 1:gomaxgeo){
      sref.g <- st.g[[geodist]]
      s.g <- st.g[[diagmotif$diamkappa+2*geodist]]
      if(orderM==2 | (orderM==3 & motif=='cycle') | motif=='clique'){motif.zimk <- cliquezimkFun(s.g,orderM,V(sref.g),V(s0.g),geodist)} else {
        motif.zimk <- motifzimkFun(tmp.subgraph,s.g,V(sref.g),V(s0.g),geodist)
      }
      motifmult <- motif.zimk[[2]]
      probk <- 1-choose(N-motifmult,n)/choose(N,n)
      thetaHT[b,1+geodist] <- sum(1/probk[probk>0])
      thetaHHmult[b,1+geodist] <- sum(motif.zimk[[1]]/probi)
    }
    }
    
    
    # Extended ancestor network = M(kappa) U beta^t_kappa(M); IWE with PIDA weights
    if(gomaxgeoPIDA>0 & includePIDA){for(geodist in 1:gomaxgeoPIDA){
      sref.g <- st.g[[geodist]]
      s.g <- st.g[[diagmotif$diamkappa+2*geodist]]
      sPIDA.g <- st.g[[diagmotif$obsdiamkappa+diagmotif$diamkappa+3*geodist]]
      if(orderM==2 | (orderM==3 & motif=='cycle') | motif=='clique'){motif.zimk.PIDA <- cliqueziUnSFun(sPIDA.g,s.g,orderM,V(sref.g),V(s0.g),geodist)} else {
        motif.zimk.PIDA <- motifziUnSFun(tmp.subgraph,sPIDA.g,s.g,V(sref.g),V(s0.g),geodist)
      }
      thetaHHPIDA[b,1+geodist] <- sum(motif.zimk.PIDA/probi)
    }
    }
    
    # sample sizes in incident samples
    if(gomaxgeo>0){tmp.n <- NULL
    for(geodist in 1:gomaxgeo){
      tmp.n <- c(tmp.n,length(V(st.g[[geodist]])))
    }
    expn <- expn + tmp.n/B
    }
  }
  
  thetapop <- countMotif(popgraph,orderM,motif)
  meanthetaHT <- apply(thetaHT,2,mean)
  meanthetamult <- apply(thetaHHmult,2,mean)
  
  
  varthetaHT <- apply(thetaHT,2,var)*(B-1)/B
  varthetamult <- apply(thetaHHmult,2,var)*(B-1)/B
  
  
  msethetaHT <- (meanthetaHT-thetapop)^2 + varthetaHT
  msethetamult <- (meanthetamult-thetapop)^2 + varthetamult
  
  
  
  if(goPIDA & includePIDA){meanthetaPIDA <- apply(thetaHHPIDA,2,mean)
  varthetaPIDA <- apply(thetaHHPIDA,2,var)*(B-1)/B
  msethetaPIDA <- (meanthetaPIDA-thetapop)^2 + varthetaPIDA
  }
  
  
  
  expn <- trunc(expn)
  cat('(Order,Motif) =', orderM,motif,'\n') 
  cat('Theta =', thetapop,'\t', '(Diameter, Obs.diam) =', diagmotif$diamkappa,diagmotif$obsdiamkappa,'\n') 
  if(gomaxgeo>0){cat('Expected n for t-wave SBS =',expn,'\n')}
  cat('estHT =', meanthetaHT,';','\t','MCvar =',varthetaHT,'\n') 
  cat('estmult =',meanthetamult,';','\t','MCvar =',varthetamult,'\n') 
  if(goPIDA & includePIDA){cat('estPIDA =', meanthetaPIDA,';','\t','MCvar =',varthetaPIDA,'\n')}
  cat('MCMSE(estHT) =', msethetaHT,'\n') 
  cat('MCMSE(estmult) =', msethetamult,'\n') 
  if(goPIDA & includePIDA){cat('MCMSE(estPIDA) =', msethetaPIDA,'\n')}
}