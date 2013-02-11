dmat.clmn <- function(x, na.rm=FALSE)
{
  Data <- matrix(colSums(x@Data) / x@dim[1L], nrow=1)
  ldim <- dim(Data)
  nprows <- base.blacs(ICTXT=x@ICTXT)$NPROW
  
  if (!is.double(Data))
    storage.mode(Data) <- "double"
  
  out <- .Call("R_dgsum2d1", as.integer(x@ICTXT), 'Col', as.integer(1L),
#  as.integer(nprows), 
                as.integer(x@ldim[2L]), Data, as.integer(1), PACKAGE="pbdBASE")
  
  return( out )
}


# ------------------
# Scale
# ------------------

setMethod("scale", signature(x="ddmatrix", center="ANY", scale="ANY"),
function(x, center=TRUE, scale=TRUE)
  {
    # Center
    if (is.ddmatrix(center)){
    
    }
    else if (is.matrix(center) || is.vector(center)){
    
    }
    else if (is.logical(center) && center){
      cntr <- dmat.clmn(x, na.rm=FALSE)
      x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
    }
    else {
      comm.print("ERROR : invalid argument for 'center'")
      stop("")
    }
    
    return( ret )
  }
)


# ------------------
# cov
# ------------------

setMethod("cov", signature(x="ddmatrix"),
function (x, y = NULL, use = "everything", method = "pearson") 
  {
    yexists <- !is.null(y)

    if (yexists){
      if (!is.ddmatrix(y)){
        comm.print("Error : 'y' must be a distributed matrix")
        stop("")
      } else if (x@dim[1] != y@dim[1]){
          comm.print("Error : incompatible dimensions")
          stop("")
        }
    }
    
    if (use=="all.obs"){
      anyna <- FALSE
      if (any(is.na(x)))
        anyna <- TRUE
      if (yexists)
        if (any(is.na(y)))
          anyna <- TRUE
      
      if (anyna){
        comm.print("Error : missing observations in cov")
        stop("")
      }
    }
    if (use=="complete.obs"){
      if (yexists){
        narows <- unique(which(is.na(rowSums(x) + rowSums(y))))
        lnarows <- length(narows)
        if (lnarows > 0) {
          if (lnarows < x@dim[1]){
            x <- x[-narows, ]
            y <- y[-narows, ]
          } else {
            comm.print("Error : no complete element pairs")
            stop("")
          }
        }
      } else
        x <- na.exclude(x)
    } else if (use=="na.or.complete"){
      comm.print("Error : na.or.complete not yet implemented")
      stop("")
    } else if (use=="pairwise.complete.obs"){
      comm.print("Error : pairwise.complete.obs not yet implemented")
      stop("")
    } else if (use!="everything"){
      comm.print("Error : invalid 'use' argument")
      stop("")
    }
    
    method <- match.arg(method)
    if (method == "pearson") {
#########################################################
      cntr <- dmat.clmn(x, na.rm=FALSE)
      x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
#      x <- scale(x, scale=FALSE)
#########################################################
      if (is.null(y))
        ret <- crossprod(x=x) / max(1, nrow(x) - 1)
      else {
        scale(y, scale=FALSE)
        ret <- base.rpdgemm(transx='T', transy='N', x=x, y=y, outbldim=x@bldim) / (nrow(x) - 1)
      }
    }
    else {
      comm.print("Error : Other methods not yet implemented")
      stop("")
    }
    
    if (x@dim[1] == 1)      
      ret[, ] <- NA
    
    return( ret )
  }
)


# Much of this wrapper taken from core R's var function
setMethod("var", signature(x="ddmatrix"),
function (x, y = NULL, na.rm = FALSE, use) 
  {
    if (missing(use)) {
      if (na.rm) 
        use <- "na.or.complete"
      else 
        use <- "everything"
    }
    
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"))
    if (is.na(na.method)) 
        stop("invalid 'use' argument")
    
    ret <- cov(x, y, na.method, FALSE)
    
    return( ret )
  }
)



setMethod("sd", signature(x="ddmatrix"),
function (x, na.rm = FALSE, reduce = FALSE, proc.dest="all") 
  {
    if (na.rm)
      x <- na.exclude(x)
    
    sdv <- .Call("R_DDMATVAR", 
                  x@Data, as.integer(x@dim[1]), 
                  as.integer(x@ldim[1]), as.integer(x@ldim[2]), 
                  as.integer(x@ICTXT),
                  PACKAGE="pbdBASE")
    
    sdv <- matrix(sqrt(sdv), nrow=1)
    
    ret <- new("ddmatrix", 
               Data=sdv, dim=c(1, x@dim[2]), 
               ldim=dim(sdv), bldim=x@bldim, ICTXT=x@ICTXT)
    
    if (reduce)
      ret <- as.vector(x=ret, proc.dest=proc.dest)
    
    return( ret )
  }
)

setMethod("sd", signature(x="ANY"), 
  function(x, na.rm = FALSE) 
    stats::sd(x=x, na.rm=na.rm)
)





