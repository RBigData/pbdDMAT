convert_dense_to_csr <- function(x)
{
    if (!is.double(x))
        storage.mode(x) <- "double"
    
    .Call("sbase_convert_dense_to_csr", x, PACKAGE="pbdDMAT")
}



convert_csr_to_dense <- function(dim, data, row_ptr, col_ind)
{
    if (!is.double(data))
      storage.mode(data) <- "double"
    
    if (!is.integer(dim))
      storage.mode(dim) <- "integer"
    
    if (!is.integer(row_ptr))
      storage.mode(row_ptr) <- "integer"
    
    if (!is.integer(col_ind))
      storage.mode(col_ind) <- "integer"
    
    .Call("sbase_convert_csr_to_dense", dim, data, row_ptr, col_ind, PACKAGE="pbdDMAT")
}
