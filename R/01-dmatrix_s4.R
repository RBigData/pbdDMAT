### S4 methods

# Misc
setGeneric(name = "print", useAsDefault = base::print, package="pbdDMAT")
setGeneric(name = "nrow", useAsDefault = base::nrow, package="pbdDMAT")
setGeneric(name = "ncol", useAsDefault = base::ncol, package="pbdDMAT")
setGeneric(name = "NROW", useAsDefault = base::NROW, package="pbdDMAT")
setGeneric(name = "NCOL", useAsDefault = base::NCOL, package="pbdDMAT")
setGeneric(name = "as.matrix", useAsDefault = base::as.matrix, package="pbdDMAT")
setGeneric(name = "na.exclude", useAsDefault = stats::na.exclude, package="pbdDMAT")
setGeneric(name = "all.equal", useAsDefault = base::all.equal, package="pbdDMAT")
setGeneric(name = "summary", useAsDefault = base::summary, package="pbdDMAT")

setGeneric(name="as.vector", 
  function(x, ...)
    standardGeneric("as.vector"),
  package="pbdDMAT"
)

setGeneric(name="rbind", 
  function(..., ICTXT=0, deparse.level=1)
    standardGeneric("rbind"),
  package="pbdDMAT"
)

setGeneric(name="cbind", 
  function(..., ICTXT=0, deparse.level=1)
    standardGeneric("cbind"),
  package="pbdDMAT"
)

setGeneric(name = "apply", useAsDefault = base::apply, package="pbdDMAT")

setGeneric(name="diag",
  function(x, ...)
    standardGeneric("diag"),
  package="pbdDMAT"
)


# Stats
setGeneric(name = "scale", useAsDefault = base::scale, package="pbdDMAT")
setGeneric(name = "var", useAsDefault = stats::var, package="pbdDMAT")
setGeneric(name = "cov", useAsDefault = stats::cov, package="pbdDMAT")
setGeneric(name = "cor", useAsDefault = stats::cor, package="pbdDMAT")
setGeneric(name = "cov2cor", useAsDefault = stats::cov2cor, package="pbdDMAT")
setGeneric(name = "prcomp", useAsDefault = stats::prcomp, package="pbdDMAT")
setGeneric(name = "scale", useAsDefault = base::scale, package="pbdDMAT")
setGeneric(name = "sweep", useAsDefault = base::sweep, package="pbdDMAT")
setGeneric(name = "lm.fit", useAsDefault = stats::lm.fit, package="pbdDMAT")

setGeneric(name="sd", 
  function(x, ...)
    standardGeneric("sd"),
  package="pbdDMAT"
)



# Reductions
#setGeneric(name = "diag", useAsDefault = base::diag, package="pbdDMAT")
setGeneric(name = "mean", useAsDefault = base::mean, package="pbdDMAT")
setGeneric(name = "median", useAsDefault = stats::median, package="pbdDMAT")
setGeneric(name = "rowSums", useAsDefault = base::rowSums, package="pbdDMAT")
setGeneric(name = "colSums", useAsDefault = base::colSums, package="pbdDMAT")
setGeneric(name = "rowMeans", useAsDefault = base::rowMeans, package="pbdDMAT")
setGeneric(name = "colMeans", useAsDefault = base::colMeans, package="pbdDMAT")



# Algebra
setGeneric(name = "t", useAsDefault = base::t, package="pbdDMAT")
setGeneric(name = "crossprod", useAsDefault = base::crossprod, package="pbdDMAT")
setGeneric(name = "tcrossprod", useAsDefault = base::tcrossprod, package="pbdDMAT")
setGeneric(name = "solve", useAsDefault = base::solve, package="pbdDMAT")
setGeneric(name = "chol", useAsDefault = base::chol, package="pbdDMAT")
setGeneric(name = "lu", useAsDefault = Matrix::lu, package="pbdDMAT")
setGeneric(name = "chol2inv", useAsDefault = base::chol2inv, package="pbdDMAT")

setGeneric(name = "norm", useAsDefault = base::norm, package="pbdDMAT")
setGeneric(name = "rcond", useAsDefault = base::rcond, package="pbdDMAT")



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
setGeneric(name="ddmatrix", 
  function(data, ...) 
    standardGeneric("ddmatrix"), 
  package="pbdDMAT"
)

setGeneric(name="as.ddmatrix", 
  function(x, ...) 
    standardGeneric("as.ddmatrix"), 
  package="pbdDMAT"
)

setGeneric(name="submatrix", 
  function(x, ...) 
    standardGeneric("submatrix"), 
  package="pbdDMAT"
)

setGeneric("submatrix<-", 
  function(x, value)
    standardGeneric("submatrix<-"),
  package="pbdDMAT"
)

setGeneric(name="ldim", 
  function(x, ...) 
    standardGeneric("ldim"), 
  package="pbdDMAT"
)

setGeneric(name="bldim", 
  function(x, ...) 
    standardGeneric("bldim"), 
  package="pbdDMAT"
)

setGeneric(name="ictxt", 
  function(x, ...) 
    standardGeneric("ictxt"), 
  package="pbdDMAT"
)

