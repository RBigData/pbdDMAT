setMethod("t", signature(x="dsmatrix"),
  function(x)
  {
    m <- x@dim[2L]
    n <- x@dim[1L]
    
    new_dim <- c(m, n)
    new_ldim <- c(dmat_ldim(m), n)
    
    out <- petsc_mattranspose(dim=x@dim, ldim=x@ldim, Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind)
    
    if (new_ldim[1L] == 0)
      Data <- matrix(0.0, 1L, 1L)
    else
      Data <- out$Data
    
    ret <- new("dsmatrix", Data=Data, dim=new_dim, ldim=new_ldim, row_ptr=out$row_ptr, col_ind=out$col_ind)
    
    return( ret )
  }
)

