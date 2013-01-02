# kappa.ddmatrix source lifted heavily from R's kappa.default function
kappa.ddmatrix <- function (z, exact = FALSE, norm = NULL, method = c("qr", "direct"), ...) 
{
  method <- match.arg(method)
  
  if (!is.null(norm)) 
    norm <- match.arg(norm, c("2", "1", "O", "I"))
  else 
    norm <- "2"
  
  if (exact && norm == "2") {
    s <- svd(z, nu = 0, nv = 0)$d
    ret <- max(s)/min(s[s > 0])
  } else {
    if (exact)
      warning(gettextf("norm '%s' currently always uses exact = FALSE", norm))
    if (norm=="2")
      norm <- "O"
    
    d <- z@dim
    
    if (method == "qr" || d[1L] != d[2L])
      z <- qr.R(qr(if (d[1L] < d[2L]) t(z) else z), complete=TRUE)
    
    if (norm=="I")
      ret <- 1/base.rpdtrcon(x=z, type=norm, uplo="L")
    else
      ret <- 1/base.rpdtrcon(x=z, type=norm, uplo="U")
  }
  
  return( ret )
}



setMethod("rcond", signature(x="ddmatrix"),
  function (x, norm = c("O", "I", "1"), triangular = FALSE, ...) 
  {
    norm <- match.arg(norm)
    d <- x@dim
    
    if (d[1L] != d[2L]){
      x <- qr.R(qr(if (d[1L] < d[2L]) t(x) else x))
      triangular <- TRUE
    }
    if (triangular) 
      ret <- base.rpdtrcon(x=x, type=norm, uplo="U")
    else 
      ret <- base.rpdgecon(x=x, type=norm)
    
    return( ret )
  }
)


