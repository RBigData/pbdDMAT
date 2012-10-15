setMethod("cov", signature(x="ddmatrix"),
function (x, y = NULL, use = "everything", method = "pearson") 
  {
    yexists <- !is.null(y)

    if (yexists){
      if (!is.ddmatrix(y)){
        comm.print("Error : 'y' must be a distributed matrix")
        stop("")
      } else if (x@dim[1] != y@dim[1]){
          comm.print("Error : incompatible dimensions")
          stop("")
        }
    }

    if (use=="all.obs"){
      anyna <- FALSE
      if (any(is.na(x)))
        anyna <- TRUE
      if (yexists)
        if (any(is.na(y)))
          anyna <- TRUE
      
      if (anyna){
        comm.print("Error : missing observations in cov")
        stop("")
      }
    }
    if (use=="complete.obs"){
      if (yexists){
        narows <- unique(which(is.na(rowSums(x) + rowSums(y))))
        lnarows <- length(narows)
        if (lnarows > 0) {
          if (lnarows < x@dim[1]){
            x <- x[-narows, ]
            y <- y[-narows, ]
          } else {
            comm.print("Error : no complete element pairs")
            stop("")
          }
        }
      } else
        x <- na.exclude(x)
    } else if (use=="na.or.complete"){
      comm.print("Error : na.or.complete not yet implemented")
      stop("")
    } else if (use=="pairwise.complete.obs"){
      comm.print("Error : pairwise.complete.obs not yet implemented")
      stop("")
    } else if (use!="everything"){
      comm.print("Error : invalid 'use' argument")
      stop("")
    }
    
    method <- match.arg(method)
    if (method == "pearson") {
        x <- scale(x, scale=FALSE)
        if (is.null(y))
            ret <- base.pdgemm(base.pdtran(x), x) / (nrow(x) - 1)
        else {
            scale(y, scale=FALSE)
            ret <- base.pdgemm(base.pdtran(x), y) / (nrow(x) - 1)
        }
    }
    else {
      comm.print("Error : Other methods not yet implemented")
      stop("")
    }

    if (x@dim[1] == 1)      
      ret[, ] <- NA

    return( ret )
  }
)

