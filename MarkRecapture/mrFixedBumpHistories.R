mrFixedBumpHistories<- function(npop0=20013, L=100, nS=100, quiet=FALSE)
{
  #
  npop<- npop0
  P<- rep(0, npop)
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
    S<- P[1:L]
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
      P[unlabelled]<- newAllocations
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
    # mix
    P<- sample(P, npop, replace=FALSE)
    #
    npopHist[k]<- npop
    if ((20 <= k) && (k <= 35))
    {
      npopPrior<- npop
      stopifnot( npopPrior == length(P) )
      eta<- (k - 20)/(35 - 20)
      npop<- round(npop0+(50000-npop0)*eta)
      extra<- npop - npopPrior
      if (0 < extra)
      {
        P<- c(P, rep(0,extra))
      }
      stopifnot( npop == length(P) )
      P<- sample(P, npop, replace=FALSE)
      cat(sprintf("k %.0f npop %.0f\n", k, npop))
      if (k == 35)
      {
        stopifnot( npop == 50000 )
        stopifnot( length(P) == 50000 )
      }
    } else if ((55 <= k) && (k <= 75))
    {
      if (k == 55)
      {
        stopifnot( npop == 50000 )
        stopifnot( length(P) == 50000 )
      }
      npopPrior<- npop
      stopifnot( npopPrior == length(P) )
      eta<- (k - 55)/(75 - 55)
      npop<- round(50000-eta*(50000-15000))
      extra<- npopPrior - npop
      if (0 < extra)
      {
        cat(sprintf("P %.0f npop %.0f extra %.0f npopPrior %.0f\n",
                    length(P), npop, extra, npopPrior))
        P<- sample(P, npop, replace=FALSE)
      }
#     stopifnot( npop == length(P) )
      if ( npop != length(P) )
      {
        cat("break pop fault:\n")
        browser()
      }
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
