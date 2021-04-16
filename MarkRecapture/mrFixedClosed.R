mrFixedClosed<- function(npop=5013, L=100, nS=30, quiet=FALSE)
{
  #
  P<- rep(0, npop)
  H<- matrix(0, nrow=nS, ncol=6)
  colnames(H)<- c("M.k", "U.k", "max.allocated", 
                  "n.unmarked", "N.star", "LPest")
  Sk<- rep(L, nS)
  maxAlloc<- 0
  marked<- vector("list", nS)
  Mk<- rep(0,nS)
  termsStar<- rep(0,nS)
  termsStar[1]<- L
  LPv<- rep(NA,nS)
  for (k in (1:nS))
  {
    S<- P[1:L]
    Sk[k]<- L
    labelled<- which(0<S)
    inThisStage<- S[labelled]
    marked[[k]]<- inThisStage
    unlabelled<- which(0 == S)
    M<- length(labelled)
    Mk[k]<- M
    if (!quiet)
    {
      cat(sprintf("\n%.0f marked subjects caught in stage %.0f\n", 
                  M, k))
      cat(sprintf("%.0f subjects caught unmarked.\n", 
                  length(unlabelled)))
    }
    unlabelledM<- length(unlabelled)
    if (0 < unlabelledM)
    {
      P[unlabelled]<- seq((1+maxAlloc),(maxAlloc+unlabelledM),1)
      maxAlloc<- unlabelledM + maxAlloc
    }
    if (1 < k)
    {
      totalMarkedTilNow<- sum(L-Mk[1:(k-1)])
      LPv[k]<- totalMarkedTilNow*L/M
      termsStar[k]<- (L+1)*(1+totalMarkedTilNow)/(1+M) - 1
    }
    H[k,]<- c(M, length(unlabelled), maxAlloc, 
              length(which(0==P)), sum(termsStar)/k, mean(LPv[2:k]))
    # mix
    P<- sample(P, npop, replace=FALSE)
  }
  NStar<- sum(termsStar)/nS
  markedCounts<- sapply(X=marked, FUN=length)
  stopifnot( sum(L-Mk[1:nS]) == length(which(P>0)) )
  return(list(H=H, npop=npop, L=L, nS=nS, n.star=NStar, marked.at.each=markedCounts, N.allocated=maxAlloc))
}
