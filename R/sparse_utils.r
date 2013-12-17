sparse_count_zeros <- function(x, tol=.Machine$double.eps)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  .Call("R_sparse_count_zeros", x, tol)
}

setMethod("sparsity", signature(x="matrix"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    count <- match.arg(tolower(count), c("zero", "other"))
    out <- match.arg(tolower(out), c("count", "proportion", "percent"))
    
    ret <- sparse_count_zeros(x=x, tol=tol)
    
    dim <- prod(dim(x))
    
    if (count == "other")
      ret <- dim - ret
    
    if (out == "proportion")
      ret <- ret/dim
    else if (out == "percent")
      ret <- ret/dim*100
    
    return( ret )
  }
)

setMethod("sparsity", signature(x="vector"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    if (!is.numeric(x))
      comm.stop("argument 'x' must be a numeric vector")
    
    dim(x) <- c(length(x), 1L)
    ret <- sparsity(x=x, count=count, out=out, tol=tol)
    
    return( ret )
  }
)

setMethod("sparsity", signature(x="ddmatrix"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    ret <- sparsity(x=x@Data, count=count, out="count", tol=tol)
    
    ret <- allreduce(ret)
    
    dim <- prod(x@dim)
    
    if (count == "other")
      ret <- dim - ret
    
    if (out == "proportion")
      ret <- ret/dim
    else if (out == "percent")
      ret <- ret/dim*100
    
    return( ret )
  }
)


