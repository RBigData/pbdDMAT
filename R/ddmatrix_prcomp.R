# Large parts of these R wrappers taken from core R prcomp.default and 
# scale.default functions
setMethod("prcomp", signature(x="ddmatrix"),
  function(x, retx=TRUE, center=TRUE, scale.=FALSE, tol=NULL) 
  {
      x <- scale(x, center = center, scale = scale.)
      cen <- attr(x@Data, "scaled:center")
      sc <- attr(x@Data , "scaled:scale")
      if (any(sc == 0))
          comm.stop("cannot rescale a constant/zero column to unit variance")
      
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


