# Matrix exponentiation using Pade' approximations and scaling and squaring from:
# "New Scaling and Squaring Algorithm for the Matrix Exponential"
# Awad H. Al-Mohy and Nicholas J. Higham, August 2009


matpow_by_squaring <- base.matpow_by_squaring
matexp_pade <- base.matexp_pade


p_matpow_by_squaring <- function(A, b=1)
{
  b <- as.integer(b)
  
  desca <- base.descinit(dim=A@dim, bldim=A@bldim, ldim=A@ldim, ICTXT=A@ICTXT)
  
  out <- base.p_matpow_by_squaring_wrap(A=A@Data, desca=desca, b=b)
  
  ret <- new("ddmatrix", Data=out, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  
  return( ret )
}


#matpow <- function(A, n)
#{
#  m <- nrow(A)
#  
##  if (n==0)
##    return(diag(1, m))
##  else if (n==1)
##    return(A)
##  
##  if (n >= 2^8)
##  {
##    E <- eigen(A)
##    B <- E$vectors %*% diag(E$values^n) %*% solve(E$vectors)
##    
##    return( B )
##  }
#  
#  B <- matpow_by_squaring(A, n)
#  
#  return(B)
#}



p_matexp_pade <- function(A)
{
  desca <- base.descinit(dim=A@dim, bldim=A@bldim, ldim=A@ldim, ICTXT=A@ICTXT)
  
  out <- base.p_matexp_pade_wrap(A=A@Data, desca=desca)
  
  N <- new("ddmatrix", Data=out$N, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  D <- new("ddmatrix", Data=out$D, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  
  R <- solve(D) %*% N
  
  return( R )
}



matexp_scale_factor <- function(x)
{
#  theta <- c(3.7e-8, 5.3e-4, 1.5e-2, 8.5e-2, 2.5e-1, 5.4e-1, 9.5e-1, 1.5e0, 2.1e0, 2.8e0, 3.6e0, 4.5e0, 5.4e0, 6.3e0, 7.3e0, 8.4e0, 9.4e0, 1.1e1, 1.2e1, 1.3e1)
  theta <- c(1.5e-2, 2.5e-1, 9.5e-1, 2.1e0, 5.4e0)
  
  
  # 1-norm
  x_1 <- max(colSums(abs(x))) 
  
  for (th in theta)
  {
    if (x_1 <= th)
      return( 0 )
  }
  
  j <- ceiling(log2(x_1/theta[5]))
  n <- 2^j
  
  return( n )
}



setMethod("expm", signature(x="matrix"), 
  function(x)
  {
    n <- matexp_scale_factor(x)
    
    if (n == 0)
      return( matexp_pade(x) )
    
    x <- x/n
    
    S <- matexp_pade(x)
    S <- matpow_by_squaring(S, n)
    
    return( S )
  }
)



setMethod("expm", signature(x="ddmatrix"), 
  function(x)
  {
    n <- matexp_scale_factor(x)
    
    if (n == 0)
      return( p_matexp_pade(x) )
    
    x <- x/n
    
    S <- p_matexp_pade(x)
    S <- p_matpow_by_squaring(S, n)
    
    return( S )
  }
)


