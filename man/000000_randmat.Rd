\name{ DistributedMatrixGeneration }
\alias{Hconst}
\alias{Hnorm}
\alias{Hunif}
\alias{Hconst.local}
\alias{Hnorm.local}
\alias{Hunif.local}
\title{ Distributed Matrix Generation }
\description{
  Routines for generating random distributed matrices.
}
\usage{
  Hconst(dim, bldim, const=0, ICTXT=0)
  Hunif(dim, bldim, min=0, max=1, ICTXT=0)
  Hnorm(dim, bldim, mean=0, sd=1, ICTXT=0)
  Hconst.local(ldim, bldim, const=0, ICTXT=0)
  Hunif.local(ldim, bldim, min=0, max=1, ICTXT=0)
  Hnorm.local(ldim, bldim, mean=0, sd=1, ICTXT=0)
}
\arguments{
  \item{dim}{global dimension.}
  \item{ldim}{local dimension.}
  \item{bldim}{blocking dimension.}
  \item{ICTXT}{BLACS context number.}
  \item{const}{numerical value for generating a constant \code{ddmatrix}.}
  \item{min, max}{numerical value for least/greatest possible value in sampling 
  from the uniform distribution.}
  \item{mean, sd}{mean and standard deviation when sampling from a 
  normal distribution.}
s}
\details{
  Each functions generates a constant or random distributed matrix
  without performing a data redistribution.  That is, the data that is 
  needed by each processor is generated locally strictly on an "as 
  needed" basis.
  
  The functions predicated with \code{.local} generate data with a 
  fixed local dimension, i.e., each processor gets an identical amount
  of data.  Likewise, the remaining functions generate a fixed global
  amount of data, and each processor may or may not have an identical 
  amount of local data.
  
  To ensure good random number generation, you should only consider 
  using these function with the \code{comm.set.seed()} function from 
  pbdBASE which relies on rlecuyer.
}
\keyword{Data Generation}




















