### TODO put in appropriate as.___.r





# -----------------------------------------------------------
# x = dsmatrix
# -----------------------------------------------------------

setMethod("as.dmat", signature(x="dsmatrix"),
  function(x)
  {
    Data <- convert_csr_to_dense(dim=x@ldim, Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind)
    llb <- new("dmat", Data=Data, dim=x@dim, ldim=x@ldim, storage="llb")
    
    return( llb )
  }
)



setMethod("as.dsvector", signature(x="dsmatrix"),
  function(x)
  {
    if (x@dim[2L] != 1)
      comm.stop("not yet supported")
    
    y <- new("dsvector", length=x@dim[1L], llength=x@ldim[1L], Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind, storage=x@storage)
    
    return( y )
  }
)



setMethod("as.matrix", signature(x="dsmatrix"),
  function(x)
  {
    y <- as.matrix(as.dmat(x))
    
    return( y )
  }
)




# -----------------------------------------------------------
# x = dmat
# -----------------------------------------------------------

setMethod("as.dsmatrix", signature(x="matrix"),
  function(x)
    as.dsmatrix(as.dmat(x))
)



setMethod("as.dsmatrix", signature(x="dmat"),
  function(x)
  {
    l <- convert_dense_to_csr(x@Data)
    sparse <- new("dsmatrix", Data=l$Data, dim=x@dim, ldim=x@ldim, row_ptr=l$row_ptr, col_ind=l$col_ind, storage="csr")
    
    return( sparse )
  }
)



setMethod("as.matrix", signature(x="dmat"),
  function(x)
  {
    mat <- matrix(0.0, x@dim[1L], x@dim[2L])
    
    dim <- x@dim
    nrows <- dim[1L]
    
    nrows.local <- dmat_ldim(nrows)
    ldim <- c(nrows.local, dim[2L])
    
    start <- dmat_index(nrows)
    end <- start + nrows.local - 1L
    
    if (ldim[1L] > 0)
      mat[start:end, ] <- x@Data
    
    # FIXME make this bcast later, too lazy atm
    mat <- allreduce(mat)
    
    return( mat )
  }
)



# -----------------------------------------------------------
# x = matrix
# -----------------------------------------------------------



setMethod("as.dmat", signature(x="matrix"),
  function(x)
  {
    dim <- dim(x)
    nrows <- dim[1L]
    
    nrows.local <- dmat_ldim(nrows)
    ldim <- c(nrows.local, dim[2L])
    
    start <- dmat_index(nrows)
    end <- start + nrows.local - 1L
    
    if (nrows.local == 0)
      Data <- matrix(0.0, 1L, 1L)
    else
      Data <- x[start:end, ]
    
    dmat <- new("dmat", Data=Data, dim=dim, ldim=ldim)
    
    return( dmat )
  }
)


