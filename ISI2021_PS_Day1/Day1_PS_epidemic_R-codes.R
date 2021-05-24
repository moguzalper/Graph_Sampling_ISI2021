# prevalence = theta, sampling fraction = f
# Y = no. cases, M = no. case-networks, K = equal-size of case-network
# kidx = index of case-networks, the first Y units in the pop
# nidx: 1. no. cases in network, 2. inclusion probability
# for size-biased sampling: lift = odds of case selection
mainEpi <- function(N=10^5,theta=0.01,f=0.01,M=10,lift=1,B=100)
{
  K = trunc(N*theta/M + 0.5); Y = M*K
  kidx = rep(0,N); kidx[1:Y] = c(t(array(1:M,c(M,K))))
  y = 1*(kidx>0)
  n0 = trunc(N*f+0.5)
  cat("(M, K, Y, n0) =",c(length(unique(kidx))-1,K,sum(y),n0),"\n")
  
  p = rep(n0/N,N); p[1:Y] = lift*p[1:Y]; p = n0*p/sum(p)
  pr.k = p.k = 1-exp(c(1:K)*log(1-max(p))); cat('Summary of pr.k under pois:','\n');print(summary(pr.k))
  if (lift==1) {
    for (m in 1:K) { p.k[m] = 1-exp(sum(log(N-m+1-c(1:n0)) - log(N+1-c(1:n0)))) }
    cat('Summary of p.k if SRS:','\n');
    print(summary(p.k))}
  nidx = array(0,c(Y,2))
  for (k in 1:M) { idk = kidx[1:Y]==k; nidx[idk,1] = sum(idk); nidx[idk,2] = pr.k[sum(idk)] }
  
  v.srs = (1-n0/N)*theta*(1-theta)/n0
  v.pois = sum((1/p-1)*y)/N^2
  cat("SE(initial) by SRS, Pois =",sqrt(c(v.srs,v.pois)),"\n")
  tmp = 1 + exp(2*K*log(1-max(p))) - 2*exp(K*log(1-max(p)))
  v.acs = c(M*(1/pr.k[K]-1), M*(M-1)*(tmp/pr.k[K]^2 -1))*K^2/N^2
  cat("SE(ACS) =",sqrt(sum(v.acs)),"\t leading term:",sqrt(v.acs[1]),"\n")
  cat("RE(analytic) =",sum(v.acs)/v.pois,"\n")
  
  mat = array(0,c(B,5))
  for (i in 1:B) {
    # s0 by SPS (sequential Poisson sampling), s0y = effective sample of cases 
    u = runif(N,0,1)/p
    s0 = sort((c(1:N)[order(u)])[1:n0])
    s0y = s0[s0<=Y]
    mat[i,3] = length(s0y)
    mat[i,5] = n0-mat[i,3]
    # initial sample estimator      
    mat[i,1] = sum(1/p[s0y])/N
    # HT estimator, need only the first unit in s0 for each case-network
    idk = s0y[!duplicated(kidx[s0y])]
    mat[i,2] = sum(nidx[idk,1]/nidx[idk,2])/N
    mat[i,4] = length(idk)
    mat[i,5] = mat[i,5] + sum(nidx[idk,1])
  }
  
  emat = rbind(colMeans(mat), sqrt(diag(var(mat))))
  colnames(emat) = c("Est.s0","Est.ACS","No cases s0",'No sample case-ntw','No sample units')
  rownames(emat) = c("MC-Mean","MC-SD")
  print(emat) 
  cat("RE(simulation) =",var(mat[,2])/c(var(mat[,1]),v.pois),"\n")
}

