# ---------------------
# THOMPSONS DATA (1991)
# ---------------------
mainPSUACS <- function(n=1){
N <- 20 # Nr of PSUs (strips) in sampling frame
Mbar <- 20 # Nr of SSUs (grids) in PSU i
psuidx <- rep(c(1:N),rep(Mbar,N))
ssuidx <- c(1:(Mbar*N))

yij <- c(rep(0,7),2,2,rep(0,17),1,11,22,3,rep(0,16),1,26,19,5,rep(0,6),1,rep(0,9),1,5,8,
         rep(0,6),17,10,5,rep(0,16),1,26,26,2,rep(0,16),1,9,6,1,rep(0,146),2,rep(0,17),14,22,2,rep(0,17),25,38,6,rep(0,17),2,3,1,rep(0,75))
grid_cases <- array(yij,c(Mbar,N))
colnames(grid_cases) <- LETTERS[1:N]

yi <- colSums(grid_cases)
ntwidx <- c(1,2,3)
yk <- c(106,105,115)

card_alphai_nepsu <- c(1,1,2,2,1,1,1,1,1,1)
nempty_psu <- colnames(grid_cases)[colSums(grid_cases)>0]
edgeik_nepsu <- data.frame(psuidx=rep(nempty_psu,card_alphai_nepsu),ntwidx=c(1,1,1,2,1,2,2,2,3,3,3,3))

card_betak <- as.vector(table(edgeik_nepsu$ntwidx))


# -------------------------------------
# Analytic variances of mean estimators
# -------------------------------------
varSRS.s0 <- (1-n/N)*var(yi)/n/Mbar/Mbar


varACS <- 0
  for(k in ntwidx){
   probk <- 1-choose(N-card_betak[k],n)/choose(N,n)
    for(l in ntwidx){
      probl <- 1-choose(N-card_betak[l],n)/choose(N,n)
      card_betakl <- length(unique(edgeik_nepsu$psuidx[edgeik_nepsu$ntwidx %in% c(k,l)]))
      probkl <- probk + probl - (1-choose(N-card_betakl,n)/choose(N,n))
    varACS <- varACS + (probkl-probk*probl)*yk[k]*yk[l]/probk/probl
    }}

varACS <- varACS/N/N/Mbar/Mbar
cat("(varSRS_s0, varHTACS) =",c(varSRS.s0,varACS),"\n")

}

