# kappa* sources lifted heavily from R's kappa.default function

###kappa.ddmatrix <- function (z, exact = FALSE, norm = NULL, method = c("qr", "direct"), ...) 
###{
###  method <- match.arg(method)
###  
###  if (!is.null(norm))
###    norm <- match.arg(norm, c("2", "1", "O", "I"))
###  else 
###    norm <- "2"
###  
###  if (exact && norm == "2") {
###    s <- svd(z, nu = 0, nv = 0)$d
###    ret <- max(s)/min(s[s > 0])
###  } else {
###    if (exact)
###      warning(gettextf("norm '%s' currently always uses exact = FALSE", norm))
###    if (norm=="2")
###      norm <- "O"
###    
###    d <- z@dim
###    
###    if (method == "qr" || d[1L] != d[2L])
###      z <- qr.R(qr(if (d[1L] < d[2L]) t(z) else z), complete=TRUE)
###    
###    if (norm=="I")
###      ret <- 1/base.rpdtrcon(x=z, type=norm, uplo="L")
###    else
###      ret <- 1/base.rpdtrcon(x=z, type=norm, uplo="U")
###  }
###  
###  return( ret )
###}

.kappa_tri2 <- function (z, exact = FALSE, norm = NULL, ...) 
{
  if (exact) {
    stopifnot(is.null(norm) || identical("2", norm))
    kappa.default(z, exact = TRUE)
  } else {
    p <- as.integer(nrow(z))
    if (is.na(p)) 
        stop("invalid nrow(x)")
    if (p != ncol(z)) 
        stop("triangular matrix should be square")
    if (is.null(norm)) 
        norm <- "1"
    else 
        1/base.rpdtrcon(x=z, type=norm, uplo="U")
  }
}

kappa.qr2 <- function (z, ...) 
{
    R <- qr.R(z, complete=FALSE)
    .kappa_tri2(R, ...)
}

kappa.ddmatrix <- function (z, exact = FALSE, norm = NULL, method = c("qr", "direct"), ...) 
{
  method <- match.arg(method)
  norm <- if (!is.null(norm)) 
            match.arg(norm, c("2", "1", "O", "I"))
          else 
            "2"
  if (exact && norm == "2") {
    s <- svd(z, nu = 0, nv = 0)$d
    max(s)/min(s[s > 0])
  } else {
    if (exact) 
      warning(gettextf("norm '%s' currently always uses exact = FALSE", norm))
    if (norm=="2")
      norm <- "O"
    d <- dim(z)
    if (method == "qr" || d[1L] != d[2L]) 
      kappa.qr2(qr(if (d[1L] < d[2L]) t(z) else z), exact = FALSE, norm = norm, ...)
    else 
      .kappa_tri2(z, exact = FALSE, norm = norm, ...)
  }
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


