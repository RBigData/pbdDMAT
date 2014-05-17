convert_dense_to_csr <- function(x)
{
    if (!is.double(x))
        storage.mode(x) <- "double"
    
    .Call("convert_dense_to_csr", x, PACKAGE="pbdDMAT")
}



convert_csr_to_dense <- function(dim, Data, row_ptr, col_ind)
{
    if (!is.double(Data))
      storage.mode(Data) <- "double"
    
    if (!is.integer(dim))
      storage.mode(dim) <- "integer"
    
    if (!is.integer(row_ptr))
      storage.mode(row_ptr) <- "integer"
    
    if (!is.integer(col_ind))
      storage.mode(col_ind) <- "integer"
    
    .Call("convert_csr_to_dense", dim, Data, row_ptr, col_ind, PACKAGE="pbdDMAT")
}
