### S4 methods

# Misc
setGeneric(name = "isSymmetric", useAsDefault = base::isSymmetric, package="pbdDMAT")
#setGeneric(name = "polyroot", useAsDefault = base::polyroot, package="pbdDMAT")



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


