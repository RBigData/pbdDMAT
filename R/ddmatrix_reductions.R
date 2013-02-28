# -------------------
# MPI-like BLACS reductions
# -------------------

# Sums
dmat.igsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call("R_igsum2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

dmat.dgsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgsum2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}


# Max value
dmat.igamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call("R_igamx2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

dmat.dgamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgamx2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}


# Min value
dmat.igamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call("R_igamn2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

dmat.dgamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgamn2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}


# point to point communication
dmat.dgesd2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgesd2d1", as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

dmat.dgerv2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgerv2d1", as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}




# Higher level reduction interface
dmat.blacsreduction <- function(x, SCOPE, op, ICTXT, proc.dest=-1)
{
  if (!is.character(SCOPE)){
    if (SCOPE==1)
      SCOPE <- 'Row'
    else if (SCOPE==2)
      SCOPE <- 'Col'
    else 
      comm.stop("ERROR : invalid argument 'scope'")
  }
  else if (SCOPE != 'Row' && SCOPE != 'Col')
    comm.stop("ERROR : invalid argument 'scope'")
  
  if (is.matrix(x)){
    m <- dim(x)[1L]
    n <- dim(x)[2L]
  }
  else if (is.vector(x)){
    m <- length(x)[1L]
    n <- 1L
  }
  else 
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (length(proc.dest==1)){
    if (proc.dest==-1){
      rdest <- cdest <- -1
    } else {
      proc.dest <- base.pcoord(ICTXT=ICTXT, PNUM=proc.dest)
      rdest <- proc.dest[[1L]]
      cdest <- proc.dest[[2L]]
    }
  }
  else {
    rdest <- proc.dest[1L]
    cdest <- proc.dest[2L]
  }
  
  if (op == 'sum'){
    if (is.integer(x))
      out <- dmat.igsum2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
    else
      out <- dmat.dgsum2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
  }
  else if (op == 'max'){
    if (is.integer(x))
      out <- dmat.igamx2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
    else
      out <- dmat.dgamx2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
  }
  else if (op == 'min'){
    out <- dmat.dgamn2d(ICTXT=ICTXT, SCOPE=SCOPE, m=m, n=n, x=x, lda=m, RDEST=rdest, CDEST=cdest)
  }
  else 
    comm.stop("ERROR : invalid argument 'op'")
  
  return( out )
}



dmat.allcolreduce <- function(x, op='sum', ICTXT=0)
{
  dmat.blacsreduction(x=x, SCOPE='Col', op=op, ICTXT=ICTXT, proc.dest=-1)
}

allcolreduce <- dmat.allcolreduce


dmat.allrowreduce <- function(x, op='sum', ICTXT=0)
{
  dmat.blacsreduction(x=x, SCOPE='Row', op=op, ICTXT=ICTXT, proc.dest=-1)
}

allrowreduce <- dmat.allrowreduce


dmat.colreduce <- function(x, op='sum', proc.dest=0, ICTXT=0)
{
  dmat.blacsreduction(x=x, SCOPE='Col', op=op, ICTXT=ICTXT, proc.dest=proc.dest)
}

colreduce <- dmat.colreduce


dmat.rowreduce <- function(x, op='sum', proc.dest=0, ICTXT=0)
{
  dmat.blacsreduction(x=x, SCOPE='Row', op=op, ICTXT=ICTXT, proc.dest=proc.dest)
}

rowreduce <- dmat.rowreduce



dmat.rcsum <- function(x, na.rm=FALSE, SCOPE, MEAN=FALSE)
{
  if (SCOPE == 'Row'){
    if (MEAN)
      Data <- rowSums(x@Data / as.double(x@dim[2L]), na.rm=na.rm)
#      Data <- matrix(rowSums(x@Data / as.double(x@dim[2L]), na.rm=na.rm), ncol=1L)
    else
      Data <- rowSums(x@Data, na.rm=na.rm)
#      Data <- matrix(rowSums(x@Data, na.rm=na.rm), ncol=1L)
    
    dim(Data) <- c(base::length(Data), 1L)
    
    if (x@dim[2L]==1)
      return( Data )
    else
      n <- nrow(Data)
  }
  else {
    if (MEAN)
      Data <- colSums(x@Data / as.double(x@dim[1L]), na.rm=na.rm)
#      Data <- matrix(colSums(x@Data / as.double(x@dim[1L]), na.rm=na.rm), nrow=1L)
    else
      Data <- colSums(x@Data, na.rm=na.rm)
#      Data <- matrix(colSums(x@Data, na.rm=na.rm), nrow=1L)
    
    dim(Data) <- c(1L, base::length(Data))
    
    if (x@dim[1L]==1)
      return( Data )
    else
      n <- ncol(Data)
  }
  
  
  nprows <- base.blacs(ICTXT=x@ICTXT)$NPROW
  
  if (!is.double(Data))
    storage.mode(Data) <- "double"
  
  out <- .Call("R_dgsum2d1", as.integer(x@ICTXT), as.character(SCOPE), as.integer(1L),
                as.integer(n), Data, as.integer(1), as.integer(-1), as.integer(-1), 
                PACKAGE="pbdBASE")
  
  return( out )
}

# -------------------
# Reductions
# -------------------

# rowSums
setMethod("rowSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Row', MEAN=FALSE)
#    dim(Data) <- c(base::length(Data), 1L)
    
    z <- new("ddmatrix", Data=Data, dim=c(x@dim[1L], 1L), ldim=c(length(x@Data), 1L), bldim=x@bldim) 
    
    return( z )
  }
)

# colSums
setMethod("colSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Col', MEAN=FALSE)
    
    z <- new("ddmatrix", Data=Data, dim=c(1L, x@dim[2L]), ldim=c(1L,length(x@Data)), bldim=x@bldim) 
    
    return( z )
  }
)

# rowMeans
setMethod("rowMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Row', MEAN=TRUE)
#    dim(Data) <- c(base::length(Data), 1L)
    
    z <- new("ddmatrix", Data=Data, dim=c(x@dim[1L], 1L), ldim=c(length(x@Data), 1L), bldim=x@bldim) 
    
    return( z )
  }
)

# colMeans
setMethod("colMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    Data <- dmat.rcsum(x, na.rm=na.rm, SCOPE='Col', MEAN=TRUE)
    
    z <- new("ddmatrix", Data=Data, dim=c(1, x@dim[2]), ldim=c(1,length(x@Data)), bldim=x@bldim) 
    
    return( z )
  }
)



# diag
setMethod("diag", signature(x="ddmatrix"),
  function(x)
  {
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    ret <- base.ddiagtk(x=x@Data, descx=descx)
    
    return( ret )
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
#    if (!is.numeric(trim) || length(trim) != 1L) 
#      comm.stop("'trim' must be numeric of length one")
    if (!base.ownany(x@dim, x@bldim, x@ICTXT))
      n <- 0
    else
    n <- length(x@Data)
    n <- pbdMPI::allreduce(n, op='sum')
#    if (trim > 0 && n) {
#        if (is.complex(x)) 
#            comm.stop("trimmed means are not defined for complex data")
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

