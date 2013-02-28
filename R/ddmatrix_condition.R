# kappa* sources lifted heavily from R's kappa.default function

.kappa_tri2 <- function (z, exact = FALSE, norm = NULL, ...) 
{
  if (exact) {
    stopifnot(is.null(norm) || identical("2", norm))
    kappa.default(z, exact = TRUE)
  } else {
    p <- as.integer(nrow(z))
    if (is.na(p)) 
        comm.stop("invalid nrow(x)")
    if (p != ncol(z)) 
        comm.stop("triangular matrix should be square")
    if (is.null(norm)) 
        norm <- "1"
    else {
      desca <- base.descinit(dim=z@dim, bldim=z@bldim, ldim=z@ldim, ICTXT=z@ICTXT)
      n <- z@dim[2L]
      
      1/base.rpdtrcon(norm=norm, uplo='U', diag='N', n=n, a=z@Data, desca=desca)
    }
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
      comm.warning(gettextf("norm '%s' currently always uses exact = FALSE", norm))
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
    if (triangular) {
      desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
      n <- x@dim[2L]
      
      ret <- base.rpdtrcon(norm=norm, uplo='U', diag='N', n=n, a=x@Data, desca=desca)
    }
    else {
      desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
      
      m <- x@dim[1L]
      n <- x@dim[2L]
      
      ret <- base.rpdgecon(norm=norm, m=m, n=n, a=x@Data, desca=desca)
    }
    
    return( ret )
  }
)


