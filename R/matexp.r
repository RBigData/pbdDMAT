matexp_pade <- function(A)
{
  n <- as.integer(nrow(A))
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  out <- .Call("R_matexp_pade", n, A, PACKAGE="pbdDMAT")
  
  N <- out$N
  D <- out$D
  
  # FIXME 
  print(N)
  print(D)
  
  R <- solve(D) %*% N
  
  return( R )
}

matpow_by_squaring <- function(A, b=1)
{
  n <- as.integer(nrow(A))
  b <- as.integer(b)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call("R_matpow_by_squaring", n, A, b, PACKAGE="pbdDMAT")
  
  return( ret )
}


# n%%2 == 0
matpow <- function(A, n)
{
  m <- nrow(A)
  
#  if (n==0)
#    return(diag(1, m))
#  else if (n==1)
#    return(A)
#  
#  if (n >= 2^8)
#  {
#    E <- eigen(A)
#    B <- E$vectors %*% diag(E$values^n) %*% solve(E$vectors)
#    
#    return( B )
#  }
  
  B <- matexp_by_squaring(A, n)
  
  return(B)
}


matexp_dumb <- function(A, n=100)
{
  X <- A
  tmp <- A
  
  for (j in 2:n)
  {
    tmp <- tmp %*% (A/j)
    X <- X + tmp
  }
  
  return( X ) 
}


# using pade' approximations
#  tmp1 <- factorial(p+q)
#  tmp2 <- factorial(p)
#  tmp3 <- factorial(q)
#  
#  TMP <- diag(1, m)
#  
#  # j == 0
#  N <- tmp1 / factorial(p+q) * TMP
#  D <- tmp1 / factorial(p+q) * TMP
#  
#  k <- min(p, q)
#  
#  for (j in 1:k)
#  {
#    TMP <- TMP %*% A
#    tmp <- factorial(p+q-j)
#    N <- N + tmp*tmp2 / (tmp1*factorial(j)*factorial(p-j)) * TMP
#    D <- D + tmp*tmp3 / (tmp1*factorial(j)*factorial(q-j)) * (-1)^j * TMP
#  }
#  
#  if (p>q)
#  {
#    for (j in (k+1):p)
#    {
#      TMP <- TMP %*% A
#      tmp <- factorial(p+q-j)
#      N <- N + tmp*tmp2 / (tmp1*factorial(j)*factorial(p-j)) * TMP
#    }
#  }
#  else if (q>p)
#  {
#    for (j in (k+1):q)
#    {
#      TMP <- TMP %*% A
#      tmp <- factorial(p+q-j)
#      D <- D + tmp*tmp3 / (tmp1*factorial(j)*factorial(q-j)) * (-1)^j * TMP
#    }
#  }
  


# Matrix exponentiation using Pade' approximations and scaling and squaring from:
# "New Scaling and Squaring Algorithm for the Matrix Exponential"
# Awad H. Al-Mohy and Nicholas J. Higham, August 2009
matexp <- function(A, t=2.1)
{
  theta <- c(3.7e-8, 5.3e-4, 1.5e-2, 8.5e-2, 2.5e-1, 5.4e-1, 9.5e-1, 1.5e0, 2.1e0, 2.8e0, 3.6e0, 4.5e0, 5.4e0, 6.3e0, 7.3e0, 8.4e0, 9.4e0, 1.1e1, 1.2e1, 1.3e1)
  
  A1 <- max(colSums(abs(A))) # 1-norm
  
  for (m in c(3, 5, 7, 9, 13))
  {
    if (A1 <= theta[m])
      return(matexp_pade(A))
  }
  
  
  j <- ceiling(log2(A1/theta[13]))
  n <- 2^j
  
  A <- t*A/n
  
  S <- matexp_pade(A)
  S <- matpow_by_squaring(S, n)
  
  return( S )
}

