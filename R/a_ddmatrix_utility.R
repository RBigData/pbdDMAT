qkmed <- function(vec)
{
  n <- as.integer(length(vec))
  
  ret <- .Fortran("QKMED", n, as.double(vec), MED=as.double(0), dup=F)$MED
  return(ret)
}



# Implementation of http://www.umiacs.umd.edu/research/EXPAR/papers/3494/node18.html#SECTION00051000000000000000
# finds k'th ordered element of the distributed matrix
# desperately needs to be rewritten in C
dmat.rank_k <- function(vec, k, shouldsort=FALSE)
{
  if (shouldsort)
    vec <- sort(vec, na.last=NA)
  
  # FIXME change to numroc call
  if (length(vec)==1)
    if (vec==0)
      vec <- NA
  
    mdmd <- median(unlist(pbdMPI::allgather(median(vec, na.rm=T))), na.rm=T)

    below <- vec[which(vec <= mdmd)]
    lbelow <- length(below)
    test <- pbdMPI::allreduce(lbelow, op='sum')

    if (test < k){
      vec <- vec[which(vec > mdmd)]
      k <- k - test
      mdmd <- dmat.rank_k(vec=vec, k=k)
    }
    else if (test > k){
      vec <- below
      mdmd <- dmat.rank_k(vec=vec, k=k)
    } else {
    
      if (lbelow==0)
        below <- -Inf
      else if (is.na(below)[1])
        below <- -Inf
      mxbl <- max(below)
      closest <- diff(c(mxbl, mdmd))
      allclosest <- pbdMPI::allreduce(closest, op='min')
      if (allclosest==closest)
        mdmd <- mxbl
      else
        mdmd <- 0
      mdmd <- pbdMPI::allreduce(mdmd, op='sum')
    }
    
  return(mdmd)
}

rank_k <- dmat.rank_k



dmat.sweep <- function(x, MARGIN, STATS, FUN="-", check.margin=FALSE)
{
  blacs_ <- blacs(x@CTXT)
  myP <- c(blacs_$MYROW, blacs_$MYCOL)
  PROCS <- c(blacs_$NPROW, blacs_$NPCOL)
  RSRC <- CSRC <- 0 # processes with first row/col of global A
  
  if(!( MARGIN==1 || MARGIN==2 )){
    comm.print("A distributed matrix is a 2 dimensional array")
    stop("")
  }
  
  if (FUN=="+")
    FUN <- 0
  else if (FUN=="-")
    FUN <- 1
  else if (FUN=="*")
    FUN <- 2
  else if (FUN=="/")
    FUN <- 3
  else {
    comm.print("Only simple arithmetic (+, -, *, /) is implemented at this time")
    stop("")
  }

  tx <- x@Data + 0 # demand a copy

  x@Data <- .Call("ddmatrix_sweep", 
         tx, STATS, 
         as.integer(x@dim), as.integer(x@bldim), 
         as.integer(PROCS), as.integer(myP), as.integer(c(RSRC, CSRC)),
         as.integer(MARGIN), as.integer(FUN),
         PACKAGE="pbdDMAT"
       )

  return(x)
}


# create a diagonal ddmatrix
dmat.ddiag <- function(x=1, nrow, ncol, bldim, ICTXT=0)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  
  dim <- c(nrow, ncol)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
  
  blacs_ <- base.blacs(ICTXT=ICTXT)
  myP <- c(blacs_$MYROW, blacs_$MYCOL)
  PROCS <- c(blacs_$NPROW, blacs_$NPCOL)
  RSRC <- CSRC <- 0L # processes with first row/col of global A
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("diag_dmat", 
               x, as.integer(dim), as.integer(ldim), as.integer(bldim),
               as.integer(PROCS), as.integer(myP), as.integer(c(RSRC, CSRC)),
               PACKAGE = "pbdDMAT")
  
  ret <- new("ddmatrix", 
             Data=out, dim=dim, ldim=ldim, bldim=bldim, CTXT=ICTXT)
  
  return( ret )
}

# Diag
setMethod("diag", signature(x="vector"), 
  function(x, nrow, ncol, type="matrix", ..., bldim=4, ICTXT=0){
    type <- match.arg(type, c("matrix", "ddmatrix"))
    
    if (type=="ddmatrix")
      ret <- dmat.ddiag(x=x, nrow=nrow, ncol=ncol, bldim=bldim, ICTXT=ICTXT)
    else
      ret <- base::diag(x=x, nrow=nrow, ncol=ncol)
    
    return( ret )
  }
)

setMethod("diag", signature(x="matrix"), 
  function(x, nrow, ncol)
    base::diag(x=x)
)
