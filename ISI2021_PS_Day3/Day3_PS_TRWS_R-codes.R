#
## Load R-packages igraph
library(igraph) # If not installed previously, use install.packages("igraph")
##

# Generate graph
genG <- function(N=100,tot=20,xi=c(0.7,0.2,0.1))
{
  amat = array(0,c(N,N))
  N11 = tot*(tot-1)/2; R11 = trunc(N11*xi[1]+0.5)
  A11off = rep(0,N11); A11off[sample.int(N11,R11,replace=F)] = 1 # Assign edges randomly to pairs of case-nodes
  A11 <- array(0,c(tot,tot))
  A11[lower.tri(A11)] <- A11off
  A11 = A11+t(A11); A11 = 1*(A11>0)
  amat[1:tot,1:tot] = A11 
  N10 = tot*(N-tot); R10 = trunc(N10*xi[2]+0.5)
  A10 = rep(0,N10); A10[sample(N10,R10,replace=F)] = 1
  amat[1:tot,c(tot+1):N] = A10
  N00 = (N-tot)*(N-tot-1)/2; R00 = trunc(N00*xi[3]+0.5)
  A00off = rep(0,N00); A00off[sample.int(N00,R00,replace=F)] = 1
  A00 <- array(0,c((N-tot),(N-tot)))
  A00[lower.tri(A00)] <- A00off
  A00 = A00+t(A00); A00 = 1*(A00>0)
  amat[c(tot+1):N,c(tot+1):N] = A00
  amat = amat+t(amat)
  amat = 1*(amat>0)
  cat("no. edges =",sum(amat)/2,"\n")
  return(amat)
}  




# Count triangles
triangle <- function(amat)
{
  N = nrow(amat); d = rowSums(amat); ridx = c(1:N)[d>1]
  theta = 0
  for (i in 1:length(ridx)) {
    adj = c(1:N)[amat[ridx[i],]==1]
    for (j in 1:length(adj)) {
      theta = theta + sum(amat[cbind(adj[j],adj[-j])])
  }}
  theta/6
}




# K-step TRW
trw <- function(amat,init=0,K=100,r=0.1)
{
  N = nrow(amat); d = rowSums(amat)
  if (init==0) { 
    p = rowSums(Amat)+r; p = p/sum(p) # Probabilities prop. to degrees
    X0 = sample.int(N,1,replace=F,prob=p) }
  if (init<0) { X0 = sample.int(N,1,replace=F) }
  if (init>0) { X0 = init }
  visited = rep(F,N); visited[X0] = T
  X = rep(X0,K+1)
  for (i in 1:K) {
    j = X[i]
    p = rep(r/(d[j]+r),N)/N  # Applying eq. (6.6)
    if (d[j]>0) { p[amat[j,]==1] = (1+r/N)/(d[j]+r) }  # Applying eq. (6.6)
    j = sample.int(N,1,replace=F,prob=p)
    visited[j] = T
    X[i+1] = j
  }
  list(X=X,visited=visited)
} 
 

# Convergence in terms of E(Y_t)
cnv.y <- function(amat,y,init=0,K=2,r=1,B=100)
{
  p = rowSums(amat)+r; p = p/sum(p)
  cat("stationary E(Y_t):",sum(y*p),"\n")
  x = rep(0,B)
  for (i in 1:B) { x[i] = trw(amat,init,K,r)$X[1+K] }
  cat("(mean, MC_SD) =", c(mean(y[x]),sqrt(var(y[x])/B)), "\n")
} 


# Generalised ratio estimator for 1st order parameter
est.y <- function(amat,y,init=0,K=100,r=0.1,B=100)
{
  N = nrow(amat); d = rowSums(amat)
  p = d+r
  mu = n = sd.w = rep(0,B)
  for (i in 1:B) {
    e = trw(amat,init,K,r); X = e$X
    n[i] = sum(e$visited)
    mu[i] = mean(y[X]/p[X])/mean(1/p[X])
    sd.w[i] = sqrt(var(y[X])/(K+1))
  }
  print(summary(cbind(n,mu,sd.w)))
  cat("(mean, SD) =", c(mean(mu),sqrt(var(mu))), "\n")
  cat("SD of IID sample:",sqrt(mean(y)*(1-mean(y))/(K+1)),"\n")
} 
 

# Subsidiary function called by the function "muFun"
checkSum <- function(x,vrb)
{
  sum(x %in% vrb)
}
# The ratio of case-triangles to the triangles with at least one non-case nodes in the population graph
muFun <- function(amat,y){
  popg <- graph_from_adjacency_matrix(amat,mode='undirected')
  tri.g <- cliques(popg,min=3,max=3)
  tmp <- as.data.frame(unlist(t(apply(sapply(tri.g, as_ids),2,sort))))
  cnt.tri <- sum(apply(tmp,1,checkSum,vrb=as_ids(V(popg)))>0)
  casetri <- sum(apply(tmp,1,FUN=function(x) sum(y[unlist(x)])==3))
  mu <- casetri/(cnt.tri-casetri)
  return(mu)
}



# Ratio estimator: case triangles vs case-noncase triangles
est.tri <- function(amat,y,init=0,K=100,r=0.1,B=100)
{
  N = nrow(amat); d = rowSums(amat)
  p = d+r
  Pmat = array(r/(d+r),c(N,N))/N # Applying eq. (6.6)
  for (i in 1:N) { if (d[i]>0) { Pmat[i,amat[i,]==1] = (1+r/N)/(d[i]+r) }} # Applying eq. (6.6)
  mu = n = ids = rep(0,B)
  for (i in 1:B) {
    e = trw(amat,init,K,r); X = e$X
    n[i] = sum(e$visited)
    idt = z1 = z0 = p.s = rep(0,K-1)
    for (k in 1:(K-1)) {
      a = X[k]; b = X[k+1]; p.s[k] = p[a]*Pmat[a,b] # Applying pi_i*pij in Section 6.4.3
      if (amat[a,b]==1) {
        adj = c(1:N)[amat[a,]*amat[b,]==1 & e$visited] # If a & b adj. find nodes adj. to both as well as selected to TRWS
        idt[k] = (length(adj) > 0) # Tells if triangle observed at step k under TRWS
        if (idt[k]) {
          for (j in 1:length(adj)) {
            h = adj[j]
            p.k = p[a]*(Pmat[a,b]+Pmat[a,h])+p[b]*(Pmat[b,a]+Pmat[b,h])+p[h]*(Pmat[h,a]+Pmat[h,b]) # Applying eqs. (6.9) and (6.14)
            z1[k] = z1[k] + y[a]*y[b]*y[h]*p.s[k]/p.k # Count triangles if all nodes are cases
            z0[k] = z0[k] + (1-y[a]*y[b]*y[h])*p.s[k]/p.k # Count triangles if any one of nodes is non-case
          }}}}
    mu[i] = sum(z1[idt==1]/p.s[idt==1])/sum(z0[idt==1]/p.s[idt==1]) # Ratio of case triangles to triangles with at least one non-case node
    ids[i] = (sum(idt)==0) | (sum(z0[idt==1])==0) # Identify replicates with invalid walks
  }
  n = n[ids!=1]; mu = mu[ids!=1] # Number of nodes visited and estimate of the ratio for replicates with valid walks
  print(summary(cbind(n,mu)))
  cat("valid walks:",length(mu),"\n")
  cat("(mean, SD) =",c(mean(mu),sqrt(var(mu))),"\nMC_SD(mean):",sqrt(var(mu)/B), "\n")
}


theta <- 20
popsize <- 100
y = c(rep(1,theta),rep(0,popsize-theta))
Amat <- genG(N=popsize,tot=theta,xi=c(0.9,0.1,0.05))
cat('Average degree of cases =',mean(colSums(Amat[1:popsize,1:theta])),'\n')
cat('Average degree of non-cases =',mean(colSums(Amat[1:popsize,(theta+1):popsize])),'\n')


# Convergence to equilibrium
# X0=1
cnv.y(amat=Amat,y,init=1,K=1,r=1,B=10000)
cnv.y(amat=Amat,y,init=1,K=1,r=0.1,B=10000)


# Convergence to equilibrium
# X0=1
cnv.y(amat=Amat,y,init=1,K=4,r=1,B=10000)
cnv.y(amat=Amat,y,init=1,K=4,r=0.1,B=10000)

# Convergence to equilibrium
# X0=1
cnv.y(amat=Amat,y,init=1,K=8,r=1,B=10000)
cnv.y(amat=Amat,y,init=1,K=8,r=0.1,B=10000)


# Convergence to equilibrium
# X0=1
cnv.y(amat=Amat,y,init=1,K=16,r=1,B=10000)
cnv.y(amat=Amat,y,init=1,K=16,r=0.1,B=10000)


# Convergence to equilibrium
# Pr(X0=i)=1/N
cnv.y(amat=Amat,y,init=-1,K=1,r=1,B=10000)
cnv.y(amat=Amat,y,init=-1,K=1,r=0.1,B=10000)


# Convergence to equilibrium
# Pr(X0=i)=1/N
cnv.y(amat=Amat,y,init=-1,K=4,r=1,B=10000)
cnv.y(amat=Amat,y,init=-1,K=4,r=0.1,B=10000)


# Convergence to equilibrium
# Pr(X0=i)=1/N
cnv.y(amat=Amat,y,init=-1,K=8,r=1,B=10000)
cnv.y(amat=Amat,y,init=-1,K=8,r=0.1,B=10000)


# Convergence to equilibrium
# Pr(X0=i)=pi_i
cnv.y(amat=Amat,y,init=0,K=1,r=1,B=10000)
cnv.y(amat=Amat,y,init=0,K=1,r=0.1,B=10000)


# Convergence to equilibrium
# Pr(X0=i)=pi_i
cnv.y(amat=Amat,y,init=0,K=4,r=1,B=10000)
cnv.y(amat=Amat,y,init=0,K=4,r=0.1,B=10000)


# est.y(amat,y,init=-1,K=100,r=0.1,B=1000)
# Initial state selected with probabilities 1/N
est.y(amat=Amat,y,init=-1,K=50,r=0.1,B=1000)


# Initial state selected with probabilities 1/N
est.y(amat=Amat,y,init=-1,K=100,r=0.1,B=1000)


# Initial state selected with probabilities 1/N
est.y(amat=Amat,y,init=-1,K=50,r=1,B=1000)


# Initial state selected with probabilities 1/N
est.y(amat=Amat,y,init=-1,K=100,r=1,B=1000)


# Number of triangles in the population graph
triangle(Amat)


# The population ratio between the total number of case-triangles and the total number of triangles with at least one non-case node
muFun(Amat,y)

# Estimation of mu
# Initial state X_0 selected with probabilities pi_i to avoid burn-in states
# est.tri(amat,y,init=0,K=100,r=0.1,B=100)
est.tri(amat=Amat,y,init=0,K=100,r=0.1,B=1000)


# Estimation of mu
# Initial state X_0 selected with probabilities pi_i to avoid burn-in states
# est.tri(amat,y,init=0,K=100,r=0.1,B=100)
est.tri(amat=Amat,y,init=0,K=500,r=0.1,B=1000)


# Estimation of mu
# Initial state X_0 selected with probabilities pi_i to avoid burn-in states
# est.tri(amat,y,init=0,K=100,r=0.1,B=100)
est.tri(amat=Amat,y,init=0,K=1000,r=0.1,B=1000)


# Average degree in the population graph
avedegree <- round(mean(colSums(Amat)),0)
avedegree


# Estimation of mu
# Initial state X_0 selected with probabilities pi_i to avoid burn-in states
# est.tri(amat,y,init=0,K=100,r=0.1,B=100)
# Let's set r = average(degree) in the graph to make a random jump on average at least as probable as an adjacent move at each time step
est.tri(amat=Amat,y,init=0,K=100,r=avedegree,B=1000)


# Estimation of mu
# Initial state X_0 selected with probabilities pi_i to avoid burn-in states
# est.tri(amat,y,init=0,K=100,r=0.1,B=100)
est.tri(amat=Amat,y,init=0,K=500,r=avedegree,B=1000)


