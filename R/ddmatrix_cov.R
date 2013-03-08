# ------------------
# cov
# ------------------

setMethod("cov", signature(x="ddmatrix"),
function (x, y = NULL, use = "everything", method = "pearson") 
  {
    yexists <- !is.null(y)
    
    if (yexists){
      if (!is.ddmatrix(y))
        comm.stop("Error : 'y' must be a distributed matrix")
      else if (x@dim[1] != y@dim[1])
        comm.stop("Error : incompatible dimensions")
    }
    
    if (use=="all.obs"){
      anyna <- FALSE
      if (any(is.na(x)))
        anyna <- TRUE
      if (yexists)
        if (any(is.na(y)))
          anyna <- TRUE
      
      if (anyna)
        comm.stop("Error : missing observations in cov")
    }
    if (use=="complete.obs"){
      if (yexists){
        narows <- unique(which(is.na(rowSums(x) + rowSums(y))))
        lnarows <- length(narows)
        if (lnarows > 0) {
          if (lnarows < x@dim[1]){
            x <- x[-narows, ]
            y <- y[-narows, ]
          } 
          else 
            comm.stop("Error : no complete element pairs")
        }
      } 
      else
        x <- na.exclude(x)
    } 
    else if (use=="na.or.complete")
      comm.stop("Error : na.or.complete not yet implemented")
    else if (use=="pairwise.complete.obs")
      comm.stop("Error : pairwise.complete.obs not yet implemented")
    else if (use!="everything")
      comm.stop("Error : invalid 'use' argument")
    
    method <- match.arg(method)
    if (method == "pearson") {
#########################################################
      cntr <- dmat.clmn(x, na.rm=FALSE)
      if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
        x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
#      x <- scale(x, scale=FALSE)
#########################################################
      if (is.null(y))
        ret <- crossprod(x=x) / max(1, nrow(x) - 1)
      else {
        scale(x=y, center=TRUE, scale=FALSE)
        ret <- t(x) %*% y / (nrow(x) - 1)
      }
    }
    else 
      comm.stop("Error : Other methods not yet implemented")
    
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
        comm.stop("invalid 'use' argument")
    
    ret <- cov(x, y, na.method, FALSE)
    
    return( ret )
  }
)


# ------------------
# sd
# ------------------

setMethod("sd", signature(x="ddmatrix"),
function (x, na.rm = FALSE, reduce = FALSE, proc.dest="all") 
  {
    if (na.rm)
      x <- na.exclude(x)
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    sdv <- base.pdclvar(x=x@Data, descx=descx)
    
    sdv <- sqrt(sdv)
    dim(sdv) <- c(1L, base::length(sdv))
#    sdv <- matrix(sqrt(sdv), nrow=1)
    
    ret <- new("ddmatrix", Data=sdv, dim=c(1, x@dim[2]), ldim=dim(sdv), bldim=x@bldim, ICTXT=x@ICTXT)
    
    if (reduce)
      ret <- as.vector(x=ret, proc.dest=proc.dest)
    
    return( ret )
  }
)

setMethod("sd", signature(x="ANY"), 
  function(x, na.rm = FALSE) 
    stats::sd(x=x, na.rm=na.rm)
)

# ------------------
# cor
# ------------------

# experimental
setMethod("cor", signature(x="ddmatrix"),
function (x, y = NULL, use = "everything", method = "pearson") 
  {
    if (method != "pearson")
      comm.stop("Not yet implemented")
    
    if (use != "everything")
      comm.stop("Not yet implemented")
    
    x <- scale(x=x, center=TRUE, scale=TRUE)
    
    if (!is.null(y)){
      yscaled <- scale(x=y, center=TRUE, scale=TRUE)
      ret <- t(x) %*% y / (nrow(x) - 1)
    }
    else {
      ret <- crossprod(x=x) / max(1, nrow(x) - 1)
    }
    
    return( ret )
  }
)



setMethod("cov2cor", signature(V="ddmatrix"),
function(V)
  {
    d <- sqrt(1/diag(V))
    
    r <- V@Data
    descv <- base.descinit(dim=V@dim, bldim=V@bldim, ldim=V@ldim, ICTXT=V@ICTXT)
    
    r <- base.pdsweep(x=r, descx=descv, vec=d, MARGIN=1L, FUN="*")
    r <- base.pdsweep(x=r, descx=descv, vec=d, MARGIN=2L, FUN="*")
    
    V@Data <- r
    
    return( V )
  }
)



