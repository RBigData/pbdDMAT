\name{ddmatrix-class}
\docType{class}
\alias{ddmatrix}
\alias{ddmatrix-class}
\title{Class ddmatrix}
\description{ Distributed matrix class. }
\section{Creating Objects}{
  \code{ new('ddmatrix', Data = ..., dim = ..., ldim = ..., 
  bldim = ..., CTXT = ...)}
}
\section{Slots}{
\describe{
\item{\code{Data}:}{LOCAL: Object of class \code{matrix}}
\item{\code{dim}:}{GLOBAL: Object of class \code{numeric}}
\item{\code{ldim}:}{LOCAL: Object of class \code{numeric}}
\item{\code{bldim}:}{GLOBAL: Object of class \code{numeric}}
\item{\code{CTXT}:}{GLOBAL: Object of class \code{numeric}}
}
}
\section{Prototype}{
\describe{
\item{matrix Data}{matrix(0)}
\item{numeric dim}{c(1,1)}
\item{numeric ldim}{c(1,1)}
\item{numeric bldim}{c(1,1)}
\item{numeric CTXT}{0}
}
}
\section{Details}{
  \code{ddmatrix} is the container for ScaLAPACK-friendly parallel
  block-cyclically distributed matrices. The class object is 
  instantiated in SPMD fashion, whereby each process owns a piece
  of the "whole" matrix (which no single R process need ever own in 
  its entirety), and each process stores its piece of the whole into a
  container with a name common to all processes.
  
  The \code{Data} slot contains the data (submatrix) belonging to
  that process. Accessible via \code{\link[=SlotAccessors]{submatrix()}}. 
  Values in the 
  \code{Data} slot will vary from process to process.
  
  The \code{dim} slot contains the global dimension; the dimension
  of the full matrix. Accessible via 
  \code{\link[=SlotAccessors]{dim()}}. The \code{dim} 
  slot is global, i.e. each process stores the same information in 
  this slot.
  
  The \code{ldim} slot contains the local dimension; here, 
  \code{all(ldim == dim(Data))}. Accessible via 
  \code{\link[=SlotAccessors]{ldim()}}. Values 
  in the \code{Data} slot will vary from process to process.
  
  The \code{bldim} slot contains the blocking factor for the block-
  cyclic distribution of the data. It consists of two numbers, 
  namely the row and column blocking, respectively. Accessible via 
  \code{\link[=SlotAccessors]{bldim()}}. 
  The \code{bldim} slot is global, i.e. each process 
  stores the same information in this slot.
  
  The \code{CTXT} slot contains the BLACS context onto which the 
  matrix information is stored. This is mostly for internal 
  bookkeeping, though advanced users might be able to effectively
  leverage differing BLACS contexts for performance improvements. 
  For details, see \code{\link{InitGrid}}. Accessible via
  \code{\link[=SlotAccessors]{ctxt()}}. 
  The \code{CTXT} slot is global, i.e. each process 
  stores the same information in this slot.
  
  A very important piece of information is that every process must 
  own something in the Data slot. This is essentially a ScaLAPACK 
  "problem", but one that is not particularly hard to avoid so long
  as you are aware that it exists. A submatrix of
  \code{matrix(0, nrow=1, ncol=1)} is used if the matrix should not 
  actually, technically, own part of the whole global matrix. You can 
  easily still see if the stored submatrix is indeed part of the global
  matrix or just a placeholder with the 
  \code{\link[=SlotAccessors]{ownany()}} function, which
  is just a wrapper on \code{numroc()} with argument \code{fixme=FALSE}.
}
\seealso{
  \code{\link{InitGrid}, \link{SlotAccessors}}
}
\keyword{Classes}
