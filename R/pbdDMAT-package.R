#' Distributed Matrix Methods
#' 
#' A package for dense distributed matrix computations. Includes the use of
#' PBLAS and ScaLAPACK libraries via pbdSLAP, communicating over MPI via the
#' BLACS library and pbdMPI.
#' 
#' \tabular{ll}{ 
#'    Package: \tab pbdDMAT \cr 
#'    Type: \tab Package \cr 
#'    License: \tab GPL \cr 
#'    LazyLoad: \tab yes \cr 
#' } 
#' 
#' This package requires an MPI library (OpenMPI, MPICH2, or LAM/MPI).
#' 
#' @useDynLib pbdDMAT,
#'   R_int_sparse_count_zeros, R_sparse_count_zeros, 
#'   R_convert_dense_to_csr, R_convert_csr_to_dense
#' 
#' @import methods pbdSLAP pbdBASE
#' @importFrom pbdMPI comm.cat comm.rank comm.print comm.size 
#'    comm.stop comm.warning allreduce barrier
#' @importFrom utils head tail
#' 
#' @name pbdDMAT-package
#' @docType package
#' @author Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen Chen, George
#' Ostrouchov, and Pragneshkumar Patel, with contributions from R Core team
#' (some wrappers taken from the base and stats packages).
#' @references Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords Package
NULL








#' Inverse from Choleski (or QR) Decomposition
#' 
#' \code{qr()} takes the QR decomposition.
#' 
#' The function returns the inverse of a choleski factored matrix, or the
#' inverse of \code{crossprod(x)} if \code{qr.R(qr(x))} is passed.
#' 
#' @name chol2inv
#' @aliases chol2inv chol2inv-method chol2inv,ddmatrix-method chol2inv
#' @docType methods
#' @param x numeric distributed matrices for
#' @param size number of columns of \code{x} containing the Choleski
#' factorization.
#' @return A numeric distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"ANY\")")}{} }
#' @seealso \code{\link{lm.fit}}
#' @keywords Methods Linear Algebra
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 3, 3)
#' 
#' R <- qr.R(qr(x))
#' xtx.inv <- chol2inv(R)
#' 
#' id <- xtx.inv %*% crossprod(x)
#' 
#' print(id)
#' 
#' finalize()
#' }
#' 
NULL





#' Compute or estimate the Condition Number of a Distributed Matrix
#' 
#' Computes or estimates the condition number.
#' 
#' 
#' @name kappa
#' @aliases ConditionNumbers kappa.ddmatrix kappa rcond-method
#' rcond,ddmatrix-method rcond
#' @docType methods
#' @param x,z numeric distributed matrices.
#' @param exact logical. Determines whether exact condition number or
#' approximation should be computed.
#' @param norm character. Determines which matrix norm is to be used.
#' @param method character. Determines the method use in computing condition
#' number.
#' @param triangular logical. If true, only the lower triangle is used.
#' @param ... Extra arguments.
#' @return Returns a number.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(z = \"ddmatrix\")")}{} }
#' @seealso \code{\link{Norm}}
#' @keywords Methods Linear Algebra ConditionNumbers
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 10, 10)
#' 
#' cnm <- rcond(x)
#' 
#' comm.print(cnm)
#' 
#' finalize()
#' }
#' 
NULL





#' Compute the Norm of a Distributed Matrix
#' 
#' Computes the norm.
#' 
#' 
#' @name Norm
#' @aliases Norm norm norm-method norm,ddmatrix-method norm
#' @docType methods
#' @param x numeric distributed matrices.
#' @param type character. Determines which matrix norm is to be used.
#' @return Returns a number.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{ConditionNumbers}}
#' @keywords Methods Linear Algebra Norm
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 10, 10)
#' 
#' nrm <- norm(x)
#' 
#' comm.print(nrm)
#' 
#' finalize()
#' }
#' 
NULL





#' Matrix Exponentiation
#' 
#' Routines for matrix exponentiation.
#' 
#' Formally, the exponential of a square matrix \code{X} is a power series:
#' 
#' \eqn{expm(X) = id + X/1! + X^2/2! + X^3/3! + \dots}
#' 
#' where the powers on the matrix correspond to matrix-matrix multiplications.
#' 
#' \code{expm()} directly computes the matrix exponential of a distributed,
#' dense matrix.  The implementation uses Pade' approximations and a
#' scaling-and-squaring technique (see references).
#' 
#' @name Expm
#' @aliases expm-method expm,ddmatrix-method expm,matrix-method expm
#' @docType methods
#' @param x A numeric matrix or a numeric distributed matrix.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"matrix\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{Reductions}, \link{MatMult},
#' \link{LinAlg}}
#' @references Matrix exponentiation using Pade' approximations and scaling and
#' squaring from: "New Scaling and Squaring Algorithm for the Matrix
#' Exponential" Awad H. Al-Mohy and Nicholas J. Higham, August 2009
#' @keywords Methods Linear Algebra
NULL





#' Sparsity of Matrix Objects
#' 
#' Determine the sparsity of a matrix, distributed, dense, or otherwise.
#' 
#' The sparsity count of a matrix is returned.
#' 
#' @name Sparsity
#' @aliases sparsity-method sparsity,vector-method sparsity,matrix-method
#' sparsity,dmat-method sparsity
#' @docType methods
#' @param x numeric matrix
#' @param count character; options are "zero" and "other". The former counts
#' the number of zeros, while the latter counts the number of non-zeros
#' ('other' elements).
#' @param out character; options are "count", "proportion", and "percent". This
#' determines whether a pure count, proportion of \code{count} elements in the
#' matrix, or percentage of \code{count} elements in the matrix.
#' @param tol numeric; the tolerance for numerical zero. This is ignored if the
#' input data is integer/logical.
#' @section Methods: \describe{ \item{list("signature(x = \"vector\")")}{}
#' \item{list("signature(x = \"matrix\")")}{} \item{list("signature(x =
#' \"dmat\")")}{} }
#' @keywords Methods,Sparse
NULL


