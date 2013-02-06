# Large parts of these R wrappers taken from core R prcomp.default and 
# scale.default functions
setMethod("prcomp", signature(x="ddmatrix"),
  function(x, retx=TRUE, center=TRUE, scale.=FALSE, tol=NULL) 
  {
      x <- scale(x, center = center, scale = scale.)
      cen <- attr(x@Data, "scaled:center")
      sc <- attr(x@Data , "scaled:scale")
      if (any(sc == 0)) {
          comm.print("cannot rescale a constant/zero column to unit variance")
          stop("")
         }
      s <- svd(x, nu=0)
      s$d <- s$d/sqrt(max(1, nrow(x) - 1))
      if (!is.null(tol)) {
          rank <- max(sum(s$d > (s$d[1L] * tol)), 1)
          if (rank < ncol(x)) {
              s$v <- s$v[, 1L:rank]
              s$d <- s$d[1L:rank]
          }
      }
      r <- list(sdev = s$d, rotation = s$v, center = if (is.null(cen)) FALSE else cen, 
          scale = if (is.null(sc)) FALSE else sc)
      if (retx) 
          r$x <- x %*% s$v
      class(r) <- "prcomp"
      
      return(r)
  }
)

setMethod("scale", signature(x="ddmatrix"),
  function(x, center=TRUE, scale=TRUE) 
  {
    if (!is.logical(center) || !is.logical(scale)){
      comm.print("argument 'scale' must be logical for a distributed matrix")
      stop("")
    }
    
    if (x@dim[1L] == 1){ # REALLY annoying special cases
      if (center){
        center <- as.vector(x)
        if (scale){
          scale <- rep(0, length=x@dim[2L])
          x@Data <- matrix(rep(NaN, length=x@ldim[2L]), nrow=1)
        }
        else {
          x@Data <- matrix(rep(0, length=x@ldim[2L]), nrow=1)
        }
      }
      else {
        if (scale){
          scale <- as.vector(x)
          x@Data <- matrix(rep(0, length=x@ldim[2L]), nrow=1)
        }
      }
    }
    else {
      if (center) {
        center <- as.vector(colMeans(x, na.rm = TRUE))
        x <- base.pdsweep(dx=x, vec=center, MARGIN=2L, FUN="-")
      }
      if (scale) {
        scale <- sqrt(as.vector(colSums(x^2))/max(1, nrow(x) - 1L))
        x <- base.pdsweep(dx=x, vec=scale, MARGIN=2L, FUN="/")
      }
    }
    
    if (is.numeric(center)) 
      attr(x@Data, "scaled:center") <- center
    if (is.numeric(scale)) 
      attr(x@Data, "scaled:scale") <- scale
    
    return( x )
  }
)

