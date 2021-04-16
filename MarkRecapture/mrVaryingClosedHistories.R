mrVaryingClosedHistories<- function(npop0=5011, L=100, nS=50, quiet=FALSE)
{
  #
  npop<- npop0
  P<- rep(0, npop)
  P.weights<- rep(1/npop, npop)
  H<- matrix(0, nrow=nS, ncol=6)
  colnames(H)<- c("M.k", "U.k", "max.allocated", 
                  "n.unmarked", "N.star", "LPest")
  Histories<- matrix(0, nrow=0, ncol=nS)
  maxAlloc<- 0
  marked<- vector("list", nS)
  Mk<- rep(0,nS)
  termsStar<- rep(0,nS)
  termsStar[1]<- L
  LPv<- rep(NA,nS)
  npopHist<- rep(0, nS)
  for (k in (1:nS))
  {
    choices<- sample.int(npop, L, replace=FALSE, prob=P.weights)
    nonchoices<- (1:npop)[-choices]
    S<- P[choices]
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
      newAllocations<- seq((1+maxAlloc),(maxAlloc+unlabelledM),1)
      P[choices[unlabelled]]<- newAllocations
      maxAlloc<- unlabelledM + maxAlloc
      Histories<- rbind(Histories, 
                        matrix(0, nrow=unlabelledM, ncol=nS))
      stopifnot( maxAlloc == nrow(Histories) )
    }
    if (1 < k)
    {
      totalMarkedTilNow<- sum(L-Mk[1:(k-1)])
      LPv[k]<- totalMarkedTilNow*L/M
      termsStar[k]<- (L+1)*(1+totalMarkedTilNow)/(1+M) - 1
      Histories[inThisStage,k]<- rep(1,M)
    }
    H[k,]<- c(M, length(unlabelled), maxAlloc, 
              length(which(0==P)), sum(termsStar)/k, mean(LPv[2:k]))
    # (no mixing needed when sampling from population with probabilities)
#   P<- sample(P, npop, replace=FALSE)
    #
    npopHist[k]<- npop
    if (11 == k)
    {
      stopifnot( 0 == (npop%%2) )
      half<- sample.int(npop, npop/2, replace=FALSE)
      otherHalf<- (1:npop)[-half]
      P.weights[half]<- 4/(3*npop)
      P.weights[otherHalf]<- 2/(3*npop)
    } else if (31 == k)
    {
      P.weights<- rep(1/npop, npop)
    }
    #
  }
  NStar<- sum(termsStar)/nS
  markedCounts<- sapply(X=marked, FUN=length)
  colnames(Histories)<- 
        sapply(X=1:ncol(Histories), 
               FUN=function(k) sprintf("stage%0.0f", k))
  # Consistency checks:
  stopifnot( sum(colSums(Histories[,1:(nS-1)])) == sum( Mk[1:(nS-1)] ) )
  # Only non-zeroes:
  nonzero<- which(0<rowSums(Histories))
  Histories<- Histories[nonzero,]
  rownames(Histories)<- sapply(X=nonzero, FUN=function(k) sprintf("R%-.0f", k))
  #
  marksMarked<- matrix(c(H[(2:nS),1], rep(L, (nS-1)), H[(1:(nS-1)),3]), nrow=(nS-1), 3, byrow=FALSE)
  colnames(marksMarked)<- c("markedInSample", "sampleSize", "allocatedDistinctMarks")
  #
  return(list(H=H, npop=npop, L=L, nS=nS, n.star=NStar, 
              marks.marked=marksMarked, npop.history=npopHist,
              marked.at.each=markedCounts, N.allocated=maxAlloc,
              histories=Histories, nonzeros=nonzero))
}
