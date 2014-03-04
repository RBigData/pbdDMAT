sparse_count_zeros <- function(x, tol=.Machine$double.eps)
{
  if (is.logical(x))
    storage.mode(x) <- "integer"
  else if (!is.integer(x) && !is.double(x))
    storage.mode(x) <- "double"
  
  if (is.integer(x))
    ret <- .Call("R_int_sparse_count_zeros", x)
  else
    ret <- .Call("R_sparse_count_zeros", x, tol)
  
  return( ret )
}



check_sparsity_inputs <- function(count, out, tol)
{
  match.arg(tolower(count), c("zero", "other"))
  match.arg(tolower(out), c("count", "proportion", "percent"))
}



calc_sparsity_return <- function(n, dim, count, out)
{
  count <- match.arg(tolower(count), c("zero", "other"))
  out <- match.arg(tolower(out), c("count", "proportion", "percent"))
  
  if (count == "other")
    n <- dim - n
  
  if (out == "count")
    ret <- n
  else if (out == "proportion")
    ret <- n/dim
  else if (out == "percent")
    ret <- n/dim*100
  
  return( ret )
}



setMethod("sparsity", signature(x="matrix"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    check_sparsity_inputs(count=count, out=out, tol=tol)
    n <- sparse_count_zeros(x=x, tol=tol)
    
    dim <- prod(dim(x))
    ret <- calc_sparsity_return(n=n, dim=dim, count=count, out=out)
    
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



setMethod("sparsity", signature(x="dmat"), 
  function(x, count="zero", out="count", tol=.Machine$double.eps)
  {
    n <- sparsity(x=x@Data, count="zero", out="count", tol=tol)
    
    n <- pbdMPI::allreduce(n)
    
    dim <- prod(x@dim)
    ret <- calc_sparsity_return(n=n, dim=dim, count=count, out=out)
    
    return( ret )
  }
)

