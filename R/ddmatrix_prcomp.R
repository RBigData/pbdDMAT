# R wrappers taken mostly unmodified from core R 
# prcomp.default and scale.default functions
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
#              s$v <- s$v[, 1L:rank, drop = FALSE]
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
      if (is.logical(center)) {
          if (center) {
              center <- as.vector(as.matrix(colMeans(x, na.rm = TRUE)))
              x <- sweep(x, 2L, center, check.margin = FALSE)
          }
      }
      else {
        comm.print("argument 'scale' must be logical for a distributed matrix")
        stop("")
      }
      if (is.logical(scale)) {
          if (scale) {
              scale <- sqrt(as.vector(as.matrix(colSums(x^2)))/max(1, nrow(x) - 1L))
              x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
          }
      }
      else {
        comm.print("argument 'scale' must be logical for a distributed matrix")
        stop("")
      }
      if (is.numeric(center)) 
          attr(x@Data, "scaled:center") <- center
      if (is.numeric(scale)) 
          attr(x@Data, "scaled:scale") <- scale
      x
  }
)


# this is extremely ad hoc and not ready for the big time; only meant to be
# used internally at the moment
setMethod("sweep", signature(x="ddmatrix"), 
  dmat.sweep
)
