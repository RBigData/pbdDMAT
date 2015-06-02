### S4 methods

# Misc
setGeneric(name = "all.equal", useAsDefault = base::all.equal, package="pbdDMAT")
setGeneric(name = "isSymmetric", useAsDefault = base::isSymmetric, package="pbdDMAT")
#setGeneric(name = "polyroot", useAsDefault = base::polyroot, package="pbdDMAT")



# Stats
setGeneric(name = "scale", useAsDefault = base::scale, package="pbdDMAT")
setGeneric(name = "var", useAsDefault = stats::var, package="pbdDMAT")
setGeneric(name = "cov", useAsDefault = stats::cov, package="pbdDMAT")
setGeneric(name = "cor", useAsDefault = stats::cor, package="pbdDMAT")
setGeneric(name = "cov2cor", useAsDefault = stats::cov2cor, package="pbdDMAT")
setGeneric(name = "prcomp", useAsDefault = stats::prcomp, package="pbdDMAT")
setGeneric(name = "scale", useAsDefault = base::scale, package="pbdDMAT")
setGeneric(name = "sweep", useAsDefault = base::sweep, package="pbdDMAT")

setGeneric(name="sd", 
  function(x, ...)
    standardGeneric("sd"),
  package="pbdDMAT"
)



# Reductions
setGeneric(name = "mean", useAsDefault = base::mean, package="pbdDMAT")
setGeneric(name = "median", useAsDefault = stats::median, package="pbdDMAT")
setGeneric(name = "rowSums", useAsDefault = base::rowSums, package="pbdDMAT")
setGeneric(name = "colSums", useAsDefault = base::colSums, package="pbdDMAT")
setGeneric(name = "rowMeans", useAsDefault = base::rowMeans, package="pbdDMAT")
setGeneric(name = "colMeans", useAsDefault = base::colMeans, package="pbdDMAT")

setGeneric(name="rowMin", 
  function(x, ...)
    standardGeneric("rowMin"),
  package="pbdDMAT"
)

setGeneric(name="rowMax", 
  function(x, ...)
    standardGeneric("rowMax"),
  package="pbdDMAT"
)

setGeneric(name="colMin", 
  function(x, ...)
    standardGeneric("colMin"),
  package="pbdDMAT"
)

setGeneric(name="colMax", 
  function(x, ...)
    standardGeneric("colMax"),
  package="pbdDMAT"
)


# Algebra
setGeneric(name = "t", useAsDefault = base::t, package="pbdDMAT")
setGeneric(name = "crossprod", useAsDefault = base::crossprod, package="pbdDMAT")
setGeneric(name = "tcrossprod", useAsDefault = base::tcrossprod, package="pbdDMAT")
setGeneric(name = "solve", useAsDefault = base::solve, package="pbdDMAT")
setGeneric(name = "chol", useAsDefault = base::chol, package="pbdDMAT")
setGeneric(name = "chol2inv", useAsDefault = base::chol2inv, package="pbdDMAT")
setGeneric(name = "norm", useAsDefault = base::norm, package="pbdDMAT")
setGeneric(name = "rcond", useAsDefault = base::rcond, package="pbdDMAT")
#setGeneric(name = "lu", useAsDefault = Matrix::lu, package="pbdDMAT")
setGeneric(name="lu", function(x, ...) standardGeneric("lu"), package="pbdDMAT")
setGeneric(name = "eigen", useAsDefault = base::eigen, package="pbdDMAT")

# Games to satisfy codetools' global variable checking --- they don't always play nice with S4
setGeneric(name="La.svd", 
  function(x, ...)
    standardGeneric("La.svd"),
  package="pbdDMAT"
)

setMethod("La.svd", signature(x="ANY"), 
  function(x, nu = min(n, p), nv = min(n, p)){
    n <- nrow(x)
    p <- ncol(x)
    
    base::La.svd(x=x, nu=nu, nv=nv)
  }
)

setGeneric(name="svd", 
  function(x, ...)
    standardGeneric("svd"),
  package="pbdDMAT"
)

setMethod("svd", signature(x="ANY"), 
  function(x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    base::svd(x=x, nu=nu, nv=nv, LINPACK=LINPACK)
  }
)

setGeneric(name="expm", 
  function(x, y, ...) 
    standardGeneric("expm"), 
  package="pbdDMAT"
)




# QR's
setGeneric(name="qr", 
  function(x, ...)
    standardGeneric("qr"),
  package="pbdDMAT"
)

setGeneric(name="qr.Q", 
  function(x, ...)
    standardGeneric("qr.Q"),
  package="pbdDMAT"
)

setGeneric(name="qr.R", 
  function(x, ...)
    standardGeneric("qr.R"),
  package="pbdDMAT"
)

setGeneric(name="qr.qy", 
  function(x, ...)
    standardGeneric("qr.qy"),
  package="pbdDMAT"
)

setGeneric(name="qr.qty", 
  function(x, ...)
    standardGeneric("qr.qty"),
  package="pbdDMAT"
)



### S4 methods for new things
setGeneric(name="dmat", 
  function(data, ...)
    standardGeneric("dmat"),
  package="pbdDMAT"
)



setGeneric("submatrix<-", 
  function(x, value)
    standardGeneric("submatrix<-"),
  package="pbdDMAT"
)

setGeneric(name="llen", 
  function(x, ...) 
    standardGeneric("llen"), 
  package="pbdDMAT"
)




setGeneric(name="sparsity", 
  function(x, count="zero", out="count", tol=.Machine$double.eps) 
    standardGeneric("sparsity"), 
  package="pbdDMAT"
)



setGeneric(name="as.dmat", 
  function(x, ...) 
    standardGeneric("as.dmat"), 
  package="pbdDMAT"
)

setGeneric(name="as.dsmatrix", 
  function(x, ...) 
    standardGeneric("as.dsmatrix"), 
  package="pbdDMAT"
)

setGeneric(name="as.dsvector", 
  function(x, ...) 
    standardGeneric("as.dsvector"), 
  package="pbdDMAT"
)


