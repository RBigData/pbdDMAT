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



# Create a diagonal distributed matrix
setMethod("diag", signature(x="vector"), 
  function(x, nrow, ncol, type="matrix", ..., bldim=.BLDIM, ICTXT=0){
    type <- match.arg(type, c("matrix", "ddmatrix"))
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    if (type=="ddmatrix")
      ret <- base.ddiagmk(x=x, nrow=nrow, ncol=ncol, bldim=bldim, ICTXT=ICTXT)
    else
      ret <- base::diag(x=x, nrow=nrow, ncol=ncol)
    
    return( ret )
  }
)

setMethod("diag", signature(x="matrix"), 
  function(x, nrow, ncol)
    base::diag(x=x)
)
