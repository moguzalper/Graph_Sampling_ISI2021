# new cases: Y12, finished cases: Y21
# networks: M12 (<M1), emering if emerge=T, or evolving otherwise 
# NB. if Y12 > M1, use Y12 divisible by M1 to ensure Y1 = Y2 
# kidx = index of case-networks, 0 for non-cases
gen.pop <- function(N=10^5,theta=c(0.01,0.01),M1=50,M12=10,Y12=100,emerge=T)
{
# at t
  K = trunc(N*theta[1]/M1 + 0.5); Y1 = M1*K
  kidx = array(0,c(N,2))
  kidx[1:Y1,1] = c(t(array(1:M1,c(M1,K))))
# at t+1  
#  Kny = trunc(min(K, Y12/M12) + 0.5); Y12 = M12*Kny
  Kny = trunc(Y12/M12+0.5); Y12 = M12*Kny
  if (emerge) { kidx[Y1+c(1:Y12),2] = c(t(array(M1+c(1:M12),c(M12,Kny)))) }
  else { kidx[Y1+c(1:Y12),2] = c(t(array(sample.int(M1,M12,replace=F),c(M12,Kny)))) }
  Y21 = N*theta[1] - (N*theta[2] - Y12)
  kidx[1:Y1,2] = kidx[1:Y1,1]
  K21 = trunc(Y21/M1)
  if (K21>0) { s = c(t(array(K*c(1:M1),c(M1,K21))) - c(1:K21) +1) }
  else { s = sample.int(Y1,Y21,replace=F) }
  kidx[s,2] = 0
  list(kidx=kidx)
}

# prevalence = theta, sampling fraction = f
# Y's = no. cases, M's = no. case-networks, K = max size of case-network
# kidx = index of case-networks, the first Y units in the pop at t=1
# n1idx, n2idx: 1. no. cases in network, 2. inclusion probability
# for size-biased sampling: lift = odds of case selection
# L1-L3: (M1, M12, Y12, emerge) = (10, 2, 200, F),(10, 2, 200, T), (10, 5, 100, T)
# M1-M3: (M1, M12, Y12, emerge) = (100, 10, 400, F),(100, 10, 400, T), (100, 10, 100, T)
# S1-S3: (M1, M12, Y12, emerge) = (500, 10, 400, F),(500, 10, 400, T), (500, 50, 100, T)
mainEpiPanel <- function(N=10^5,theta=c(0.01,0.01),M1=50,M12=10,Y12=100,emerge=T,f=0.01,lift=1,B=50,go=F)
{
  kidx = gen.pop(N=N,theta=theta,M1=M1,M12=M12,Y12=Y12,emerge=emerge)$kidx
  Y = colSums(kidx>0)
  y = 1*(kidx>0) 
  K = c(max(table(kidx[y[,1]>0,1])), max(table(kidx[y[,2]>0,2])))
  cat("(M, K) =",c(length(unique(c(kidx)))-1, K),"\t Y =",Y,"\n")

  n0 = trunc(N*f+0.5)
  p = rep(n0/N,N); p[1:Y[1]] = lift*p[1:Y[1]]; p = n0*p/sum(p); p.1 = max(p)
  pr.k = 1-exp(c(1:max(K))*log(1-p.1)); print(summary(pr.k))
  n1idx = n2idx = t(array(c(min(p)),c(2,N)))
  for (k in 1:M1) { idk = kidx[,1]==k; n1idx[idk,1] = sum(idk); n1idx[idk,2] = pr.k[sum(idk)] }
  for (k in unique(kidx[,2])) { if (k>0) {
    idk = kidx[,2]==k; n2idx[idk,1] = sum(idk); n2idx[idk,2] = pr.k[sum(idk)] }}
  
  v.srs = (1-n0/N)*var(y[,2]-y[,1])/n0
  v.pois = sum((1/p-1)*(y[,2]-y[,1])^2)/N^2
  cat("SE(iniital) by SRS, Pois =",sqrt(c(v.srs,v.pois)),"\n")
  tmp = 1 + exp(2*K[1]*log(1-p.1)) - 2*exp(K[1]*log(1-p.1))
  v1.acs = c(M1*(1/pr.k[K[1]]-1), M1*(M1-1)*(tmp/pr.k[K[1]]^2 -1))*K[1]^2/N^2
  cat("SE(ACS_1) =",sqrt(sum(v1.acs)),"\t leading term:",sqrt(v1.acs[1]),"\n")
  i2k = unique(kidx[kidx[,2]>0,2]); v2.acs = 0
  for (k in i2k) { m = sum(kidx[,2]==k); v2.acs = v2.acs + (1/pr.k[m]-1)*m^2/N^2 }
  cat("SE(ACS_2) =",sqrt(v2.acs),"\n")
  if (go) {  
    i1k = unique(kidx[kidx[,1]>0,1]); cov.acs = 0
    for (k1 in i1k) { for (k2 in i2k) {
      m1 = sum(kidx[,1]==k1); m2 = sum(kidx[,2]==k2); m12 = sum(kidx[,1]==k1 | kidx[,2]==k2)
      tmp = 1 - exp(m1*log(1-p.1)) - exp(m2*log(1-p.1)) + exp(m12*log(1-p.1))
      cov.acs = cov.acs + (tmp/(pr.k[m1]*pr.k[m2]) -1)*m1*m2/N^2 }}
    v.pacs = sum(v1.acs) + v2.acs - 2*cov.acs
    cat("SE(ACS_panel) =",sqrt(v.pacs),"\t RE(analytic) =",sum(v.pacs)/v.pois,"\n")
  }
  
  mat = array(0,c(B,3)); t1 = rep(0,B)
  p1 = p; p1[kidx[,1]>0] = n1idx[kidx[,1]>0,2]
  for (i in 1:B) {
# s0 by SPS (sequential Poisson sampling), s0y = effective sample of cases 
    u = runif(N,0,1)/p
    s0 = sort((c(1:N)[order(u)])[1:n0])
# HT estimator under panel, based on initial sample s0      
    mat[i,1] = sum((y[s0,2]-y[s0,1])/p[s0])/N
# HT estimator under panel ACS, based on samples s(0) and A(t+1)
    if (sum(y[s0,1])>0) { s0y = s0[y[s0,1]>0]
      idk = s0y[!duplicated(kidx[s0y,1])]
      t1[i] = sum(n1idx[idk,1]/n1idx[idk,2])/N }
    if (sum(y[s0,2])>0) { s0y = s0[y[s0,2]>0]
      idk = s0y[!duplicated(kidx[s0y,2])]
      mat[i,3] = sum(n2idx[idk,1]/n2idx[idk,2])/N }
    mat[i,2] = mat[i,3] - t1[i]
# iterated ACS, based on samples s(t) and A(t+1)
    s.t = s0[y[s0,1]==0]
    if (sum(y[s0,1])>0) { s0y = s0[y[s0,1]>0]
      idk = unique(kidx[s0y,1])
      for (k in idk) { s.t = c(s.t, c(1:N)[kidx[,1]==k]) }} 
    mat[i,3] = sum(y[s.t,2]/p1[s.t])/N - t1[i]
  }
  
  emat = rbind(colMeans(mat), sqrt(diag(var(mat))))
  colnames(emat) = c("Panel","pACS","iACS")
  rownames(emat) = c("Mean","SE")
  print(emat)
  cat("RE(pACS, iACS) =",emat[2,2:3]^2/emat[2,1]^2,"\n")
}

