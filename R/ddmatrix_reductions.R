# -------------------
# Reductions
# -------------------

dmat.rcsum <- function(x, na.rm=FALSE, SCOPE, MEAN=FALSE)
{
  if (SCOPE == 'Row'){
    if (MEAN)
      Data <- matrix(rowSums(x@Data / as.double(x@dim[2L]), na.rm=na.rm), ncol=1L)
    else
      Data <- matrix(rowSums(x@Data, na.rm=na.rm), ncol=1L)
    
    if (x@dim[2L]==1)
      return( Data )
    
    n <- nrow(Data)
  }
  else {
    if (MEAN)
      Data <- matrix(colSums(x@Data / as.double(x@dim[1L]), na.rm=na.rm), nrow=1L)
    else
      Data <- matrix(colSums(x@Data, na.rm=na.rm), nrow=1L)
    
    if (x@dim[1L]==1)
      return( Data )
    
    n <- ncol(Data)
  }
    
#    if (SCOPE=='Row' && MEAN) comm.print(Data, all.rank=T)
    
  
  nprows <- base.blacs(ICTXT=x@ICTXT)$NPROW
  
  if (!is.double(Data))
    storage.mode(Data) <- "double"
  
  out <- .Call("R_dgsum2d1", as.integer(x@ICTXT), as.character(SCOPE), as.integer(1L),
                as.integer(n), Data, as.integer(1), PACKAGE="pbdBASE")
  
  return( out )
}

# rowSums
setMethod("rowSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    z <- new("ddmatrix", dim=c(x@dim[1L], 1L), ldim=c(length(x@Data), 1L), bldim=x@bldim) 
    z@Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Row', MEAN=FALSE)
    return( z )
  }
)

# colSums
setMethod("colSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    z <- new("ddmatrix", dim=c(1L, x@dim[2L]), ldim=c(1L,length(x@Data)), bldim=x@bldim) 
    z@Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Col', MEAN=FALSE)
    return( z )
  }
)

# rowMeans
setMethod("rowMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    z <- new("ddmatrix", dim=c(x@dim[1], 1), ldim=c(length(x@Data), 1), bldim=x@bldim) 
    z@Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Row', MEAN=TRUE)
    return( z )
  }
)

# colMeans
setMethod("colMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    z <- new("ddmatrix", dim=c(1, x@dim[2]), ldim=c(1,length(x@Data)), bldim=x@bldim) 
    z@Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Col', MEAN=TRUE)
    return( z )
  }
)



# diag
setMethod("diag", signature(x="ddmatrix"),
  function(x)
  {
    base.ddiagtk(x)
  }
)

# sum
setMethod("sum", signature(x="ddmatrix"),
  function(x, ..., na.rm=FALSE)
  {
    # no need to correct for local storage issues
    other <- list(...)
    if (length(other) > 0)
      other <- sum(
        sapply(other, 
          function(i) {
            if (is.ddmatrix(i)) 
              sum(i@Data, na.rm=na.rm) 
            else {
              if (comm.rank()==0)
                sum(i, na.rm=na.rm)
              else
                0
            }
          }
        ), 
      na.rm=na.rm)
    else
      other <- 0
    local <- sum(x@Data, na.rm=na.rm) + other
    pbdMPI::allreduce(local, op="sum")
  }
)

# mean, with large chunks taken from base:::mean.default
setMethod("mean", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (na.rm) 
        x@Data <- matrix(x@Data[!is.na(x@Data)])
#    if (!is.numeric(trim) || length(trim) != 1L) {
#      comm.print("'trim' must be numeric of length one")
#      stop("")
#    }
    if (!base.ownany(x@dim, x@bldim, x@ICTXT))
      n <- 0
    else
    n <- length(x@Data)
    n <- pbdMPI::allreduce(n, op='sum')
#    if (trim > 0 && n) {
#        if (is.complex(x)) {
#            comm.print("trimmed means are not defined for complex data")
#            stop("")
#          }
#        if (any(is.na(x))) 
#            return(NA_real_)
#        if (trim >= 0.5) 
#            return(median(x, na.rm = FALSE))
##        lo <- floor(n * trim) + 1
##        hi <- n + 1 - lo
##        x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
#    }
    
    sum(x, na.rm=na.rm) / n
  }
)

# prod
setMethod("prod", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      prod <- prod(x@Data, na.rm=na.rm)
    else
      prod <- 1
    pbdMPI::allreduce(prod, op="prod")
  }
)

# min/max
setMethod("min", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      min <- min(x@Data, na.rm=na.rm)
    else
      min <- Inf
    pbdMPI::allreduce(min, op="min")
  }
)

setMethod("max", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      max <- max(x@Data, na.rm=na.rm)
    else
      max <- -Inf
    pbdMPI::allreduce(max(x@Data, na.rm=na.rm), op="max")
  }
)

setMethod("median", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (!na.rm){
      test <- any(is.na(x@Data))
      test <- pbdMPI::allreduce(test, op='max')
      if (test>0)
        return(NA)
    } else
      x@Data <- matrix(x@Data[!is.na(x@Data)])
    lenloc <- length(x@Data)
    if (!base.ownany(x@dim, x@bldim, x@ICTXT))
      lenloc <- 0
    n <- pbdMPI::allreduce(lenloc, op='sum')
    if (n%%2==1)
      ret <- dmat.rank_k(vec=x@Data, k=ceiling(n/2), shouldsort=T)
    else {
      ret1 <- dmat.rank_k(vec=x@Data, k=ceiling(n/2), shouldsort=T)
      ret2 <- dmat.rank_k(vec=x@Data, k=ceiling((n+1)/2), shouldsort=T)
      ret <- mean(c(ret1, ret2))
    }

    return( ret )
  }
)

