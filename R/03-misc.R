#---------------------------------------------
# redistribution tools
#---------------------------------------------

# dropper function, used in subsetting
base.dropper <- function(x, oldbldim, iorj, ij, ICTXT)
{
  blacs_ <- base.blacs(ICTXT)

  bldim <- oldbldim #x@bldim
  if (x@ICTXT != ICTXT){
# FIXME: would like to alter block dimension so that the data
# rebalances, but the code below causes a horrible null pointer
# problem for some cases.
#    if (ICTXT==1)
#      bldim <- c(dim(x)[1], ceiling(bldim[2] / blacs_$NPCOL))
#    if (ICTXT==2)
#      bldim <- c(ceiling(bldim[1] / blacs_$NPROW), dim(x)[2])
  
    newObj <- dmat.reblock(dx=x, bldim=bldim, ICTXT)
  }
  
  if (iorj=='i'){ # rows
    if (newObj@ldim[1L] == newObj@dim[1L]){
      new <- newObj@Data[ij, ]
      if (base::length(new)==0)
        new <- matrix(0.0)
      if (!is.matrix(new))
        dim(new) <- c(length(ij), newObj@ldim[2L])
#          new <- matrix(new, ncol=newObj@ldim[2])
      
      newObj@Data <- new
    }
  } else { # columns
    if (newObj@ldim[2L] == newObj@dim[2L]){
      new <- newObj@Data[, ij]
      if (base::length(new)==0)
        new <- matrix(0.0)
      if (!is.matrix(new))
        dim(new) <- c(newObj@ldim[1L], length(ij))
#          new <- matrix(new, nrow=newObj@ldim[1])
      
      newObj@Data <- new
    }
  }
  
  if (iorj=='i'){
    if (ij[1L] > 0)
      newObj@dim[1L] <- base::length(ij)
    else
      newObj@dim[1L] <- newObj@dim[1L] - base::length(ij)
  } else {
    if (ij[1L] > 0)
      newObj@dim[2L] <- base::length(ij)
    else
      newObj@dim[2L] <- newObj@dim[2L] - base::length(ij)
  }
  
  newObj@ldim <- dim(newObj@Data)
  
  return(newObj)
}

dropper <- base.dropper



# redistribute data from one BC type to another
dmat.reblock <- function(dx, bldim=dx@bldim, ICTXT=.ICTXT)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  
  dim <- dx@dim
  m <- dim[1]
  n <- dim[2]
  xattrs <- attributes(dx@Data)
  
  ldimB <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
  TldimB <- ldimB # true ldimB
  
  # lda's of 1 infuriate pdgemr2d
  mxx <- pbdMPI::allreduce(max(dx@ldim), op='max')
  mxb <- pbdMPI::allreduce(max(ldimB), op='max')
  
  if (all(dx@ldim==1))
    dx@ldim[1] <- mxx
  if (all(ldimB==1))
    ldimB[1] <- mxb
  
  if (pbdMPI::allreduce(dx@ldim[1], op='max')==1 && dx@dim[1]>1)
    dx@ldim[1] <- mxx
  if (pbdMPI::allreduce(ldimB[1], op='max')==1)
    ldimB[1] <- mxb
  
  descx <- base.descinit(dim=dim, bldim=dx@bldim, ldim=dx@ldim, ICTXT=dx@ICTXT)
  descy <- base.descinit(dim=dim, bldim=bldim, ldim=ldimB, ICTXT=ICTXT)
  
  dy <- new("ddmatrix", Data=matrix(0.0, 1, 1), dim=dim, ldim=TldimB, bldim=bldim, ICTXT=ICTXT)
  
  xblacs_ <- base.blacs(dx@ICTXT)
  if (xblacs_$MYROW==-1 || xblacs_$MYCOL==-1){
    descx[2] <- -1
  }
  
  blacs_ <- base.blacs(ICTXT=ICTXT)
  if (blacs_$MYROW==-1 || blacs_$MYCOL==-1){
    descy[2] <- -1
  }
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  
  ret <- base.rpdgemr2d(x=dx@Data, descx=descx, descy=descy)
  
  dy@Data <- ret
  
  
  if (length(xattrs) > 1){
    battrs <- union(attributes(dy@Data), xattrs[-1])
    names(battrs) <- names(xattrs)
    attributes(dy@Data) <- battrs
  }
  
  return( dy )
}

reblock <- dmat.reblock


dmat.redistribute <- function(dx, bldim=dx@bldim, ICTXT=.ICTXT)
{
#  if (dx@ICTXT != ICTXT)
  ret <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=ICTXT)
  
  return( ret )
}

redistribute <- dmat.redistribute


#---------------------------------------------
# *bind functions
#---------------------------------------------

base.rbind <- function(..., ICTXT=.ICTXT)
{
  args <- list(...)
  
  return( base.rbind2(args=args, ICTXT=ICTXT) )
}

base.rbind2 <- function(args, ICTXT=.ICTXT)
{ 
#  args <- list(...)
  
  oldctxt <- args[[1]]@ICTXT
  
  args <- lapply(args, 
    FUN=function(dx) dmat.redistribute(dx=dx, bldim=dx@bldim, ICTXT=1)
  )
  
  dim <- c(sum(sapply(args, function(x) dim(x)[1])), args[[1]]@dim[2])
  bldim <- args[[1]]@bldim
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=1, fixme=TRUE)
  
  Data <- lapply(args, submatrix)
  
  ret <- new("ddmatrix", Data=Reduce(base::rbind, Data), dim=dim, ldim=ldim, bldim=bldim, ICTXT=1)
  
  if (ICTXT!=1)
    ret <- dmat.redistribute(dx=ret, bldim=ret@bldim, ICTXT=ICTXT)
  
  return( ret )
}

base.cbind <- function(..., ICTXT=.ICTXT)
{
  args <- list(...)
  
  oldctxt <- args[[1]]@ICTXT
  
  args <- lapply(args, 
    FUN=function(dx) dmat.redistribute(dx=dx, bldim=dx@bldim, ICTXT=2)
  )
  
  dim <- c(args[[1]]@dim[1], sum(sapply(args, function(x) dim(x)[2])))
  bldim <- args[[1]]@bldim
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=2, fixme=TRUE)
  
  Data <- lapply(args, submatrix)
  
  ret <- new("ddmatrix", Data=Reduce(base::cbind, Data), dim=dim, ldim=ldim, bldim=bldim, ICTXT=2)
  
  if (ICTXT!=2)
    ret <- dmat.redistribute(dx=ret, bldim=ret@bldim, ICTXT=ICTXT)
  
  return( ret )
}


#---------------------------------------------
# rank (median)
#---------------------------------------------

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



#---------------------------------------------
# other
#---------------------------------------------

# checking compatibility between distributed matrices for use with
# scalapack/pblas. For internal use only.
base.checkem <- function(x, y, checks=1:3)
{
  # All dimension equal
  if (1 %in% checks)
    if (any(x@dim!=y@dim))
      comm.stop("Error: non-conformable distributed arrays")
  # Same BLACS context
  if (2 %in% checks)
    if (x@ICTXT != y@ICTXT)
      comm.stop("Error: Distributed matrices 'x' and 'y' must belong to the same BLACS context")
  # Same blocking dimension
  if (3 %in% checks)
    if (any(x@bldim != y@bldim))
      comm.stop("Distributed matrices 'x' and 'y' must have the same block dimension.")
}

checkem <- base.checkem


# print first few entries of the global matrix
base.firstfew <- function(dx, atmost=5)
{
  blacs_ <- base.blacs(dx@ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL
  NPROW <- blacs_$NPROW
  NPCOL <- blacs_$NPCOL

  MB <- dx@bldim[1]
  NB <- dx@bldim[2]

  if (prod(dx@dim) < atmost)
    atmost <- prod(dx@dim)

  dim <- c( min(dx@dim[1], atmost), min(dx@dim[2], atmost) )

  out <- numeric(atmost)
  ct <- 1
  for (j in 1:dim[2]-1){
    for (i in 1:dim[1]-1){
      l <- floor(i / (NPROW * MB))
      m <- floor(j / (NPCOL * NB))
      
      pr <- (0 + floor(i/MB)) %% NPROW
      pc <- (0 + floor(j/NB)) %% NPCOL
      
      if (MYROW==pr && MYCOL==pc){
        x <- 1 + i %% MB;
        y <- 1 + j %% NB;
        out[ct] <- dx@Data[x+MB*l,y+NB*m]

        ct <- ct+1
      }
      ct <- pbdMPI::allreduce(ct, op='max')
      if (ct == atmost+1)
         break
      barrier()
    }
    if (ct == atmost+1)
      break
  }
  barrier()
  out <- pbdMPI::allreduce(out, op='sum')
  return(out)
}

