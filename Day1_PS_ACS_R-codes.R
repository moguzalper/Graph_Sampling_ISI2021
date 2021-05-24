# Extended data set from Thompson (1990)
# -------------------------
# Function mainACS
# --------------------
mainACS <- function(first.sample=c(0,2),choose20=T){

yi <- c(1,0,20,2,10,1000)
  
N <- length(yi)
n <- 2
  
pri <- n/N
  
idx_F <- yi
idx_omega <- yi
  
# Strategy with B
edgeik <- data.frame(i=c(1,0,20,20,2,10,10,10,1000,1000,1000),k=c(1,0,20,2,2,2,10,1000,2,10,1000))

# Strategy with B*
edgeik_star <- data.frame(i=c(1,0,20,2,10,10,1000,1000),k=c(1,0,20,2,10,1000,10,1000))


# Strategy with Bdagger, with ancestor ntw. {20}
edgeik_dagger_opta <- data.frame(i=c(1,0,20,20,10,10,1000,1000),k=c(1,0,20,2,10,1000,10,1000))

# Strategy with Bdagger, with ancestor ntw. {10,1000}
edgeik_dagger_optb <- data.frame(i=c(1,0,20,10,10,10,1000,1000,1000),k=c(1,0,20,2,10,1000,2,10,1000))

  
# Multiplicities
card_betak_star <- NULL
card_betak_dagger_opta <- NULL
card_betak_dagger_optb <- NULL
for(k in idx_omega){
  card_betak_star <- c(card_betak_star,sum(edgeik_star$k %in% k))  
  card_betak_dagger_opta <- c(card_betak_dagger_opta,sum(edgeik_dagger_opta$k %in% k))  
  card_betak_dagger_optb <- c(card_betak_dagger_optb,sum(edgeik_dagger_optb$k %in% k))  
}
  

  
all.subsets <- combn(N,n)
nr.all.subsets <- dim(all.subsets)[2]

  
tmp.subset <- c(1:nr.all.subsets)[apply(all.subsets,2,FUN = function(x) sum(x %in% which(idx_F %in% first.sample))==2)]
B <- c(1:nr.all.subsets)
muHT_s0 <- NULL
muHT <- NULL
muHT_SD <- NULL

s1_all <- rep(list(0),nr.all.subsets)


for(b in 1:nr.all.subsets){
if(b==1){s <- tmp.subset}


if(b>1 & b<nr.all.subsets){B <- c(1:nr.all.subsets)[!(c(1:nr.all.subsets) %in% tmp.subset)]
s <- sample(B,1)
tmp.subset <- c(tmp.subset,s)
}

if(b==nr.all.subsets){B <- c(1:nr.all.subsets)[!(c(1:nr.all.subsets) %in% tmp.subset)]
  s <- B}

s0 <- idx_F[all.subsets[,s]]

# HTE based on s0
muHT_s0 <- c(muHT_s0,mean(s0))

# Strategy I, B
  s1 <- unique(edgeik$k[edgeik$i %in% s0])
  
# Strategy II, Bstar
  s1_star <- unique(edgeik_star$k[edgeik_star$i %in% s0])

#  Strategy I, Modified HTE & Strategy II, HTE
  prk <- 1-choose(N-card_betak_star,n)/choose(N,n)
  tmp.muHT <- sum(yi[idx_F %in% s1_star]/prk[idx_F %in% s1_star])/N

# Sample-dependent strategy
  s1_SD <- s1_star
  prk_SD <- prk
  sample.str <- 'Bstar'

if(sum(intersect(c(20, 10, 1000), first.sample))>0){
if(sum(intersect(c(20,10,1000),first.sample) %in% 20)==1){
  s1_SD <- unique(edgeik_dagger_opta$k[edgeik_dagger_opta$i %in% s0])
  prk_SD <- 1-choose(N-card_betak_dagger_opta,n)/choose(N,n)
  sample.str <- 'Bdagger'
} }

if(sum(intersect(c(10, 1000), first.sample))>0){  
if(sum((intersect(c(10,1000),first.sample) %in% c(10,1000)))>0){
    s1_SD <- unique(edgeik_dagger_optb$k[edgeik_dagger_optb$i %in% s0])
    prk_SD <- 1-choose(N-card_betak_dagger_optb,n)/choose(N,n)
    sample.str <- 'Bdagger'
} }

if(sum(c(20,10,1000) %in% first.sample)>0){
# if(20 %in% (intersect(c(20,10,1000),first.sample[first.sample %in% c(20,10,1000)]))){
if(sum(first.sample %in% c(20,10))==2 | sum(first.sample %in% c(20,1000))==2){
  if(choose20==T){
    s1_SD <- unique(edgeik_dagger_opta$k[edgeik_dagger_opta$i %in% s0])
    prk_SD <- 1-choose(N-card_betak_dagger_opta,n)/choose(N,n)
    sample.str <- 'Bdagger'
  }
  
  if(choose20==F){
    s1_SD <- unique(edgeik_dagger_optb$k[edgeik_dagger_optb$i %in% s0])
    prk_SD <- 1-choose(N-card_betak_dagger_optb,n)/choose(N,n)
    sample.str <- 'Bdagger'
  }
} }

tmp.muHT_SD <- sum(yi[idx_F %in% s1_SD]/prk_SD[idx_F %in% s1_SD])/N  

s1_all[[b]] <- s1

muHT <- c(muHT,tmp.muHT)
muHT_SD <- c(muHT_SD, tmp.muHT_SD)  

}


# Rao-Blackwellisation for strategy I
muHT_modified_RB <- NULL
for(j in 1:nr.all.subsets){
  # j <- 14
  sel <- sapply(s1_all,FUN = function(x) sum(s1_all[[j]] %in% x)==length(s1_all[[j]]) & length(s1_all[[j]])==length(x))
  muHT_modified_RB <- c(muHT_modified_RB,mean(muHT[sel]))
}

cat("first.sample =",first.sample,"\n")
cat("sample_dependent_strategy =",sample.str,"\n")
cat("pop.mean =",mean(yi),"\n")


result_ACS <- t(array(c(mean(muHT_s0),mean(muHT),mean(muHT_modified_RB),mean(muHT),mean(muHT_SD),c(var(muHT_s0),var(muHT),var(muHT_modified_RB),var(muHT),var(muHT_SD))*(nr.all.subsets-1)/nr.all.subsets),c(5,2)))
colnames(result_ACS) <- c('HT_s0','ACS_HTmodif','ACS_HTmodif_RB','ACS_HT_Bstar','ACS_HT_SD')
rownames(result_ACS) <- c('ExpectedValue','Variance')
 
print(result_ACS)

return(list(muHT_SD=muHT_SD))

}




