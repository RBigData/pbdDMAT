# -------------------
# Reductions
# -------------------

# rowSums
setMethod("rowSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    z <- new("ddmatrix", dim=c(x@dim[1], 1), ldim=c(x@bldim[1],1), bldim=x@bldim) 
    z@Data <- matrix(base.blacs.sum('Row', x@Data, x@dim, na.rm=na.rm))
    return( as.vector(base.as.matrix(z)) )
  }
)

# colSums
setMethod("colSums", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    ldim <- base.numroc(x@dim, x@bldim)
    z <- new("ddmatrix", dim=c(1, x@dim[2]), ldim=ldim, bldim=x@bldim) 
    z@Data <- matrix(base.blacs.sum('Column', x@Data, x@dim, na.rm=na.rm), nrow=1)
    return( as.vector(base.as.matrix(z)) )
  }
)

# rowMeans
setMethod("rowMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    z <- new("ddmatrix", dim=c(x@dim[1], 1), ldim=c(x@bldim[1],1), bldim=x@bldim) 
    z@Data <- matrix(base.blacs.sum('Row', x@Data, x@dim, na.rm=na.rm, means=TRUE, num=x@dim[2]))
    return(as.vector(base.as.matrix(z)))
  }
)

# colMeans
setMethod("colMeans", signature(x="ddmatrix"), 
  function(x, na.rm=FALSE){
    z <- new("ddmatrix", dim=c(1, x@dim[2]), ldim=c(1,x@dim[2]), bldim=x@bldim) 
    z@Data <- matrix(base.blacs.sum('Column', x@Data, x@dim, na.rm=na.rm, means=TRUE, num=x@dim[1]), nrow=1)
    return(as.vector(base.as.matrix(z)))
  }
)

# diag
setMethod("diag", signature(x="ddmatrix"),
  function(x)
  {
    blacs_ <- base.blacs(x@CTXT)
    myP <- c(blacs_$MYROW, blacs_$MYCOL)
    PROCS <- c(blacs_$NPROW, blacs_$NPCOL)
    RSRC <- CSRC <- 0 # processes with first row/col of global A
    ISRCPROC <- 0
    
    dim <- x@dim
    ldim <- x@ldim
    bldim <- x@bldim
    
    out <- .Call("diag_grab", 
           x@Data,
           dim=as.integer(dim),
           bldim=as.integer(bldim),
           gP=as.integer(PROCS),
           myP=as.integer(myP),
           SRC=as.integer(c(RSRC, CSRC)),
           PACKAGE="pbdDMAT"
         )
    
    pbdMPI::allreduce(out, op='sum')
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
    if (!base.ownany(x@dim, x@bldim, x@CTXT))
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
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
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
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
      min <- min(x@Data, na.rm=na.rm)
    else
      min <- Inf
    pbdMPI::allreduce(min, op="min")
  }
)

setMethod("max", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
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
    if (!base.ownany(x@dim, x@bldim, x@CTXT))
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




