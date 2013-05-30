### General redistribution
# redistribute data from one context and/or blocking factor to another
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
  
#  if (pbdMPI::allreduce(dx@ldim[1], op='max')==1 && dx@dim[1]>1)
#    dx@ldim[1] <- mxx
#  if (pbdMPI::allreduce(ldimB[1], op='max')==1)
#    ldimB[1] <- mxb
  
  descx <- base.descinit(dim=dim, bldim=dx@bldim, ldim=dx@ldim, ICTXT=dx@ICTXT)
  descy <- base.descinit(dim=dim, bldim=bldim, ldim=ldimB, ICTXT=ICTXT)
  
  dy <- new("ddmatrix", Data=matrix(0.0, 1, 1), dim=dim, ldim=TldimB, bldim=bldim, ICTXT=ICTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  ret <- base.rpdgemr2d(x=dx@Data, descx=descx, descy=descy)
  
  dy@Data <- ret
  
  if (length(xattrs) > 1){
    xattrs$dim <- dy@ldim
    attributes(dy@Data) <- xattrs
  }
  
  
  return( dy )
}

reblock <- dmat.reblock
dmat.redistribute <- dmat.reblock
redistribute <- dmat.redistribute




### Simple interfaces

# cyclic
dmat.as.rowcyclic <- function(dx, bldim=.BLDIM)
{
  if (dx@ICTXT == 2 && all(dx@bldim == bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=2L)
  
  return( ret )
}


dmat.as.colcyclic <- function(dx, bldim=.BLDIM)
{
  if (dx@ICTXT == 1 && all(dx@bldim == bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=1L)
  
  return( ret )
}

as.rowcyclic <- dmat.as.rowcyclic
as.colcyclic <- dmat.as.colcyclic


# block-cyclic
dmat.as.blockcyclic <- function(dx, bldim=.BLDIM)
{
  if (dx@ICTXT == 0L && all(dx@bldim == bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=0L)
}

as.blockcyclic <- dmat.as.blockcyclic


# block
dmat.as.block <- function(dx, square.bldim=TRUE)
{
  
  blacs_ <- base.blacs(ICTXT=dx@ICTXT)
  procs <- c(blacs_$NPROW, blacs_$NPCOL)
  
  if (square.bldim)
    new.bldim <- rep(max(sapply(1L:2L, function(i) ceiling(dx@dim[i]/procs[i]))), 2L)
  else
    new.bldim <- sapply(1L:2L, function(i) ceiling(dx@dim[i]/procs[i]))
  
  if (all(dx@bldim == new.bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=new.bldim, ICTXT=dx@ICTXT)
}

as.block <- dmat.as.block


dmat.as.rowblock <- function(dx)
{
  new.bldim <- rep(ceiling(dx@dim[1L]/comm.size()), 2L)
  
  if (dx@ICTXT == 2 && all(dx@bldim == new.bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=new.bldim, ICTXT=2L)
}

dmat.as.colblock <- function(dx)
{
  new.bldim <- rep(ceiling(dx@dim[2L]/comm.size()), 2L)
  
  if (dx@ICTXT == 1 && all(dx@bldim == new.bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=new.bldim, ICTXT=1L)
}

as.rowblock <- dmat.as.rowblock
as.colblock <- dmat.as.colblock

