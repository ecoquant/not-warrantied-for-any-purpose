# wlscbsblm.R : code supporting ``When linear systems can't be solved by linear means''
# Jan Galkowski, 11th-14th June 2018.
# Last changed 22nd June 2020.

library(nloptr)
library(random)

randomizeSeed<- function(external=FALSE)
{
  #set.seed(31415)
  # Futz with the random seed
  if (!external)
  {
    rf<-  (as.numeric(Sys.time())*100000)%%99989
    s<- round(rf)
    set.seed(s)
    cat(sprintf("Seed chosen is: %.0f, internal.\n", s))
  } else
  {
    set.seed(randomNumbers(n=1, min=1, max=10000, col=1, base=10, check=TRUE))
  }
  return( sample.int(2000000, size=sample.int(2000, size=1), replace=TRUE)[1] )
}

wonkyRandom<- randomizeSeed(external=TRUE)

A<- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3, byrow=TRUE)

b<- matrix(c(12,4,16), 3, 1)

print(solve(A,b))

print(lm(b ~ A + 0))


B<- matrix(c(b, matrix(c(20, 11, 99), 3, 1)), 3, 2)

B<- cbind(B, matrix(c(101, -1, 10, 200, 3, 9), 3, 2))
 
print(solve(A,B)) 

b1 <- rbind(b, 23)

A1<- rbind(A, c(17, -2, 11))

print(lm(b1 ~ A1 + 0))

L2norm<- function(x)
{
  sqrt( sum(x*x) )
}

L2normForMatrix<- function(X, scaling=1)
{
  stopifnot( is.matrix(X) )
  X.s<- X/scaling
  summing<- sum( colSums(X.s*X.s) )
  stopifnot (0 <= summing)
  yields<- sqrt( summing )
  return( yields )
}

r<- function(sigma, alpha, beta)
{
  stopifnot( (0 <= sigma) && (sigma <= 1) )
  stopifnot( alpha < beta )
  alpha*(1 - sigma) + beta*sigma
}

# Recall original was:
#
# A<- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3, byrow=TRUE)


P1.func<- function(x, alpha.beta)
{
  stopifnot( is.vector(x) )
  stopifnot( 3 == length(x) )
  #
  sigma11<- x[1]
  sigma23<- x[2]
  sigma31<- x[3]
  alpha11<- alpha.beta[1]
  beta11<- alpha.beta[2]
  alpha23<- alpha.beta[3]
  beta23<- alpha.beta[4]
  alpha31<- alpha.beta[5]
  beta31<- alpha.beta[6]
  #
  P1<- matrix( c( r(sigma11,alpha11,beta11),  2,                          3,
                  2,                          1,                          r(sigma23,alpha23,beta23),
                  r(sigma31,alpha31,beta31),  r(sigma23,alpha23,beta23),  1
                ),
              nrow=3, ncol=3, byrow=TRUE )
  return(P1)
}


objective1<- function(x, alpha.beta)
{
  stopifnot( is.vector(x) )
  stopifnot( 3 == length(x) )
  b<- matrix(c(12,4,16),3,1)
  x.right<- matrix(c(-3,6,1),3,1)
  P1<- P1.func(x, alpha.beta)
  d<- b - P1 %*% x.right
  # L2 norm
  return( L2norm(d) )
}

nloptr.options1<- list("algorithm"="NLOPT_GN_ISRES", "xtol_rel"=1.0e-4, "print_level"=0, "maxeval"=100000, "population"=1000)

alpha.beta<- c(-2, 2, -1, 8, -6, 6)

Y1<- nloptr(x0=rep(0.5,3), 
            eval_f=objective1,
            lb=rep(0,3), ub=rep(1,3),
            opts=nloptr.options1,
            alpha.beta=alpha.beta
          )
          
print(Y1)
          
cat(sprintf("Y1 resulting estimates for a_{11}, a_{23}, and a_{31} are: %.2f, %.2f, %2.f\n", 
            r(Y1$solution[1], alpha.beta[1], alpha.beta[2]), r(Y1$solution[2], alpha.beta[3], alpha.beta[4]), 
            r(Y1$solution[3], alpha.beta[5], alpha.beta[6])))
            
A2<- matrix(c(0.8, 0.05, 0.15, 
              0.2, 0.6,  0.2, 
              0.2, 0.3,  0.5
             ),
             nrow=3, ncol=3, byrow=TRUE)
             

            
P2.func<- function(x)
{
  # Sunny, Rainy, Foggy
  stopifnot( is.vector(x) )
  stopifnot( 2 == length(x) )
  #
  a.31<- x[1]
  a.32<- x[2]
  #
  P2<- matrix( c( 0.8,  0.05, 0.15,
                  0.2,  0.6,  0.2,
                  a.31, a.32, 0.5 
                ),
              nrow=3, ncol=3, byrow=TRUE )
  return(P2)
}

objective2<- function(x)
{
  stopifnot( is.vector(x) )
  stopifnot( 2 == length(x) )
  x.right<- matrix(c(1020, 301, 155), 3, 1)
  b<- matrix(c(854, 416, 372),3,1)
  P2<- P2.func(x)
  d<- b - P2 %*% x.right
  # L2 norm
  return( L2norm(d) )
}

constraint2<- function(x)
{
  stopifnot( 2 == length(x) )
  return( (x[1] + x[2] - 0.5 ))
}

nloptr.options2<- list("algorithm"="NLOPT_GN_ISRES", "xtol_rel"=1.0e-4, "print_level"=0, "maxeval"=100000, "population"=1000)

Y2<- nloptr(x0=rep(0.5,2), 
            eval_f=objective2,
            eval_g_eq=constraint2,
            lb=rep(0,2), ub=rep(1,2),
            opts=nloptr.options2
          )
          
print(Y2)
          
cat(sprintf("Y2 resulting estimates for a_{31}, a_{32} are: %.2f, %.2f\n", 
            Y2$solution[1], Y2$solution[2]))
            

P3.func<- function(x)
{
  # Sunny, Rainy, Foggy
  stopifnot( is.vector(x) )
  stopifnot( 3 == length(x) )
  #
  a.31<- x[1]
  a.32<- x[2]
  # There's an x[3] but it isn't used in the P3.func. See
  # the objective3.
  #
  P3<- matrix( c( 0.8,  0.05, 0.15,
                  0.2,  0.6,  0.2,
                  a.31, a.32, 0.5 
                ),
              nrow=3, ncol=3, byrow=TRUE )
  return(P3)
}

objective3<- function(x)
{
  stopifnot( is.vector(x) )
  stopifnot( 3 == length(x) )
  x.right<- matrix(c(1020, r(x[3], 155, 1020), 155), 3, 1)
  b<- matrix(c(854, 416, 372),3,1)
  P3<- P3.func(x)
  d<- b - P3 %*% x.right
  # L2 norm 
  return( L2norm(d) )
}

constraint3<- function(x)
{
  stopifnot( 3 == length(x) )
  return( (x[1] + x[2] - 0.5 ))
}

nloptr.options3<- list("algorithm"="NLOPT_GN_ISRES", "xtol_rel"=1.0e-4, "print_level"=0, "maxeval"=100000, "population"=1000)

Y3<- nloptr(x0=rep(0.5,3), 
            eval_f=objective3,
            eval_g_eq=constraint3,
            lb=rep(0,3), ub=rep(1,3),
            opts=nloptr.options3
          )
          
print(Y3)
          
cat(sprintf("Y3 resulting estimates for a_{31}, a_{32}, and eta are: %.2f, %.2f, %.2f, with that eta corresponding to %.0f\n", 
            Y3$solution[1], Y3$solution[2], Y3$solution[3], r(Y3$solution[3], 155, 1020)))
            
A4<- matrix(c(0.663, 0.138, 0.934,
              0.928, 0.631, 0.710,
              0.334, 0.478, 0.791
             ), 3, 3, byrow=TRUE)


objective4<- function(x)
{
  stopifnot( is.vector(x) )
  stopifnot( 9 == length(x) )
  B<- matrix(c(5215, 13693,  7265, 4217, 9367, 10588, 14372, 12043,
               7528, 17825, 11024, 4989, 14860, 9447, 16162, 13087,
               6161, 12798,  7702, 3023,  9551, 8908, 11429,  8734
              ), 3, 8, byrow=TRUE)
  X.right<- matrix(c(1356, 7505, 4299, 3419, 7132, 1965, 8365, 8031,
                     5689, 8065, 7001,  638, 8977, 1088, 3229, 1016,
                     3777, 8135, 3689, 1993, 3635, 9776, 8967, 7039
                     ), 3, 8, byrow=TRUE)
  P4<- matrix(x, nrow=3, ncol=3, byrow=TRUE)
  d<- B - P4 %*% X.right
  # L2 norm for matrix
  return( L2normForMatrix(d, scaling=1000) )
}


nloptr.options4<- list("algorithm"="NLOPT_GN_ISRES", "xtol_rel"=1.0e-6, "print_level"=0, "maxeval"=300000, "population"=1000)

Y4<- nloptr(x0=rep(0.5,9), 
            eval_f=objective4,
            lb=rep(0,9), ub=rep(1,9),
            opts=nloptr.options4
          )
          
print(Y4)
          
cat("Y4 resulting estimates for mathbf{A}:\n")
print(matrix(Y4$solution, 3, 3, byrow=TRUE))
