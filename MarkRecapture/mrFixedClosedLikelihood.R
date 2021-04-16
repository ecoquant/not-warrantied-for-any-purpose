
logLikelihood<- function(N.hat, L, Mv, k)
{
  P<- function(i, L, Mv) if (2>i) { 0 } else {i*L - sum(Mv[2:i])}
  lL<- 0
  for (i in (1:k))
  {
    P.i<- P(i, L, Mv)
    # x := number of marked drawn without replacement
    # m := number of marked in population of size N.hat
    # n := number of unmarked in population of size N.hat
    # k := number drawn
    # (Note the built-in stats "dhyper" doesn't tolerate double
    #  datatypes well on its parameters. It expects whole numbers.)
    lL<- lL + dhyper(x=Mv[i], m=P.i, n=round(N.hat-P.i), k=L, log=TRUE)
  }
  return(lL)
}

mrFixedClosedLikelihoodFrom<- function(N, L, nS, Chapman.estimate, 
                                       markedSequence, quiet=FALSE)
{
  #
  if (!quiet)
  {
    cat(sprintf(
    "mrFixedClosedLikelihoodFrom(N=%.0f, L=%.0f, nS=%.0f)\n", 
                N, L, nS))
    cat("counts of marked in sequence:\n")
    print(markedSequence)
    cat("\n")
  }
  stopifnot( nS == length(markedSequence) )
  #
  
  yields<- optimize(f=logLikelihood, 
                    lower=L,  upper=2*Chapman.estimate,
                    maximum=TRUE, tol = .Machine$double.eps^0.25,
                    L=L, Mv=markedSequence, k=nS)
  #
  N.hat<- round(yields$maximum)
  p.hat<- 1/N.hat
  #
  return(list(N=N, L=L, nS=nS, N.hat=N.hat, p.hat=p.hat))
}
