# ---------------------------------------------------------
# lm.fit
# ---------------------------------------------------------

setMethod("lm.fit", signature(x="ddmatrix", y="ddmatrix"), 
  function (x, y, tol = 1e-07, singular.ok=TRUE)
  {
    # checks
    base.checkem(x=x, y=y, checks=2L:3L)
    if (x@dim[1L] != y@dim[1L])
      comm.stop("Error : incompatible dimensions")
    
    oldctxt <- y@ICTXT
    
    # Matrix descriptors
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=oldctxt)
    descb <- base.descinit(dim=y@dim, bldim=y@bldim, ldim=y@ldim, ICTXT=oldctxt)
    
    m <- desca[3L]
    n <- desca[4L]
    nrhs <- descb[4L]
    
    # fit the model
    out <- base.rpdgels(tol=tol, m=m, n=n, nrhs=nrhs, a=x@Data, desca=desca, b=y@Data, descb=descb)
    
    if (!singular.ok && out$RANK < x@dim[2L]) 
      comm.stop("singular fit encountered")
    
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- out$A
    if (base.ownany(dim=y@dim, bldim=y@bldim, ICTXT=y@ICTXT))
      y@Data <- out$B
    
    
    eff <- new("ddmatrix", Data=out$EFF, dim=y@dim, 
               ldim=y@bldim, bldim=y@bldim, ICTXT=y@ICTXT)
    fitted.values <- new("ddmatrix", Data=out$FT, dim=y@dim,
                         ldim=y@ldim, bldim=y@bldim, ICTXT=y@ICTXT)
    residuals <- new("ddmatrix", Data=out$RSD, dim=y@dim,
                     ldim=y@ldim, bldim=y@bldim, ICTXT=y@ICTXT)
    
    # rearranging solution in the overdetermined and/or rank deficient case
    temp <- 1L:n # indexing of coefficients
    if (m >= n)
    {
      y <- y[temp, , ICTXT=1L]
    } 
    else 
    {
      cdim <- c(n-y@dim[1L], y@dim[2L])
      cldim <- base.numroc(dim=cdim, bldim=y@bldim, ICTXT=y@ICTXT, fixme=TRUE)
      c <- new("ddmatrix", Data=matrix(as.double(NA), nrow=cldim[1L], ncol=cldim[2L]), dim=cdim, ldim=cldim, bldim=y@bldim, ICTXT=y@ICTXT)
      y <- base.rbind(y, c, ICTXT=1L)
    }
    
    # convert IPIV to global vector if it isn't already
    if (base.blacs(ICTXT=x@ICTXT)$NPCOL > 1L)
    {
      dim(out$IPIV) <- c(1L, length(out$IPIV))
      c <- new("ddmatrix", Data=out$IPIV,
                dim=c(1L, y@dim[1L]), ldim=c(1L, y@ldim[1L]), 
                bldim=y@bldim, ICTXT=oldctxt)
      pivot <- as.vector(c)
    } 
    else 
    {
      pivot <- out$IPIV
    }
    
    
    if (out$RANK < n)
    {
  #    vec <- as.ddmatrix(matrix(NA, nrow=1, ncol=nrhs), bldim=y@bldim)
      if (m >= n)
        y[(out$RANK+1L):n, , ICTXT=y@ICTXT] <- as.double(NA)
      else 
      {
        if (out$RANK < m)
          y[(out$RANK+1L):m, , ICTXT=y@ICTXT] <- as.double(NA)
      }
      if (any(pivot - temp != 0L))
      {
        perm <- sapply(temp, function(i) temp[which(i==pivot)])
        y <- dmat.redistribute(dx=y, bldim=y@bldim, ICTXT=2L)
        y <- y[perm, , ICTXT=oldctxt]
      } 
      else 
      {
        y <- dmat.redistribute(dx=y, bldim=y@bldim, ICTXT=oldctxt)
      }
    } 
    else 
    {
      y <- dmat.redistribute(dx=y, bldim=y@bldim, ICTXT=oldctxt)
    }
    
    # rownames
  #  if (base.ownany(dim=y@dim, bldim=y@bldim, ICTXT=y@CTXT)){
  #    coords <- sapply(temp, function(i) base.g2l_coord(ind=i, dim=y@dim, bldim=y@bldim, ICTXT=y@CTXT)[5])
  #    mycoords <- coords[which(!is.na(coords))]
  #    
  #    rownames(y@Data) <- paste("x", mycoords, sep="")
  #  }
    
    qr <- list(qr=x, tau=out$TAU, pivot=pivot, tol=tol, rank=out$RANK)
    
    attr(qr, "class") <- "qr"
    
    ret <- list(coefficients=y, residuals=residuals, effects=eff, 
                rank=out$RANK, fitted.values=fitted.values, assign=attr(x@Data, "assign"),
                qr=qr, df.residual=(x@dim[1] - out$RANK))
    
    return( ret )
  }
)


