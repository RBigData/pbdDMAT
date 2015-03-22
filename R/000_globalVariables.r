### Do not delete this file and file name!!
### This file should be loaded before all other *.r files.

### This is to avoid the fale positive messages from R CMD check.
###   "no visible binding for global variable"
### Suggested by Prof Brian Ripley
### ?globalVariables

utils::globalVariables(c(".BLDIM", ".ICTXT", ".conflicts.OK"))


#' Some default parameters for pbdDMAT.
#' 
#' This set of controls is used to provide default values in this package.
#' 
#' The default blocking \code{.BLDIM} is \code{c(4,4)}, which results in a 4 by
#' 4 blocking dimension for distributed matrices.  Any time a function takes
#' the \code{bldim=} argument, it will default to this value unless the user
#' specifies an alternative.
#' 
#' The default ICTXT is 0.  This is the full 2-dimensional processor grid.
#' 
#' @name pbdDMAT Control
#' @rdname control
#' @docType data
#' @format Objects contain several parameters for communicators and methods.
NULL

#' @rdname control
.BLDIM <- c(4, 4)

#' @rdname control
.ICTXT <- 0

.conflicts.OK <- TRUE
