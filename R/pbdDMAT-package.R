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








#' Arithmetic Reductions: Sums, Means, and Prods
#' 
#' Arithmetic reductions for distributed matrices.
#' 
#' Performs the reduction operation on a distributed matrix.
#' 
#' There are four legitimately new operations, namely \code{rowMin()},
#' \code{rowMax()}, \code{colMin()}, and \code{colMax()}.  These
#' implementations are not really necessary in R because one can easily (and
#' reasonably efficiently) do something like
#' 
#' \code{apply(X=x, MARGIN=1L, FUN=min, na.rm=TRUE)}
#' 
#' But \code{apply()} on a \code{ddmatrix} is \emph{very} costly, and should be
#' used sparingly.
#' 
#' @name Reductions
#' @aliases Reductions sum-method sum,ddmatrix-method sum mean-method
#' mean,ddmatrix-method mean median-method median,ddmatrix-method median
#' prod-method prod,ddmatrix-method prod rowSums-method rowSums,ddmatrix-method
#' rowSums colSums-method colSums,ddmatrix-method colSums rowMeans-method
#' rowMeans,ddmatrix-method rowMeans colMeans-method colMeans,ddmatrix-method
#' colMeans min-method min,ddmatrix-method min max-method max,ddmatrix-method
#' max rowMin-method rowMin,ddmatrix-method rowMin rowMax-method
#' rowMax,ddmatrix-method rowMax colMin-method colMin,ddmatrix-method colMin
#' colMax-method colMax,ddmatrix-method colMax
#' @docType methods
#' @param x numeric distributed matrix
#' @param na.rm logical. Should missing (including \code{NaN}) be removed?
#' @param ... additional arguments
#' @return Returns a global numeric vector.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{Arithmetic}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- sum(colMeans(x))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Arithmetic Operators
#' 
#' Binary operations for distributed matrix/distributed matrix and distributed
#' matrix/vector operations.
#' 
#' If \code{x} and \code{y} are distributed matrices, then they must be
#' conformable, on the same BLACS context, and have the same blocking
#' dimension.
#' 
#' @name Arithmetic
#' @aliases Arithmetic +-method +,ddmatrix,ddmatrix-method
#' +,ddmatrix,numeric-method +,numeric,ddmatrix-method + --method
#' -,ddmatrix,ddmatrix-method -,ddmatrix,missing-method
#' -,ddmatrix,numeric-method -,numeric,ddmatrix-method - *-method
#' *,ddmatrix,ddmatrix-method *,ddmatrix,numeric-method
#' *,numeric,ddmatrix-method * /-method /,ddmatrix,ddmatrix-method
#' /,ddmatrix,numeric-method /,numeric,ddmatrix-method / ^-method
#' ^,ddmatrix,ddmatrix-method ^,ddmatrix,numeric-method
#' ^,numeric,ddmatrix-method ^ %%-method %%,ddmatrix,ddmatrix-method
#' %%,ddmatrix,numeric-method %%,numeric,ddmatrix-method %% %/%-method
#' %/%,ddmatrix,ddmatrix-method %/%,ddmatrix,numeric-method
#' %/%,numeric,ddmatrix-method %/%
#' @docType methods
#' @param x,y numeric distributed matrices or numeric vectors
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\", y =
#' \"ddmatrix\")")}{} \item{list("signature(x = \"numeric\", y =
#' \"ddmatrix\")")}{} \item{list("signature(x = \"ddmatrix\", y =
#' \"numeric\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{LinAlg}, \link{MatMult}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- (2*x) - x^(.5)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Matrix Multiplication
#' 
#' Multiplies two distributed matrices, if they are conformable.
#' 
#' \code{x} and \code{y} must be conformable, on the same BLACS context, but
#' they need not be blocked with the same blocking dimension. The return will
#' default to the blocking dimension of \code{x}.
#' 
#' If you need to use \code{x} and \code{y} with differing blocking dimensions
#' and you want the return to have blocking different from that of \code{x},
#' then use the function \code{base.rpdgemm()}.
#' 
#' The \code{crossprod()} and \code{tcrossprod()} functions behave exactly as
#' their R counterparts.
#' 
#' @name MatMult
#' @aliases MatMult %*%-method %*%,ddmatrix,ddmatrix-method %*%
#' crossprod-method crossprod,ddmatrix-method crossprod,ddmatrix,ANY-method
#' crossprod tcrossprod-method tcrossprod,ddmatrix-method
#' tcrossprod,ddmatrix,ANY-method tcrossprod
#' @docType methods
#' @param x,y numeric distributed matrices
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\", y =
#' \"ddmatrix\")")}{} \item{list("signature(x = \"ddmatrix\", y = \"ANY\")")}{}
#' }
#' @seealso \code{\link{Arithmetic}, \link{LinAlg}, \link{MatMult}}
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
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- x %*% x
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Miscellaneous Mathematical Functions
#' 
#' Binary operations for distributed matrix/distributed matrix and distributed
#' matrix/vector operations.
#' 
#' Performs the miscellaneous mathematical calculation on a distributed matrix.
#' 
#' @name MiscMath
#' @aliases MiscMath abs-method abs,ddmatrix-method abs sqrt-method
#' sqrt,ddmatrix-method sqrt exp-method exp,ddmatrix-method exp log-method
#' log,ddmatrix-method log log2-method log2,ddmatrix-method log2 log10-method
#' log10,ddmatrix-method log10 log1p-method log1p,ddmatrix-method log1p
#' sin-method sin,ddmatrix-method sin cos-method cos,ddmatrix-method cos
#' tan-method tan,ddmatrix-method tan asin-method asin,ddmatrix-method asin
#' acos-method acos,ddmatrix-method acos atan-method atan,ddmatrix-method atan
#' sinh-method sinh,ddmatrix-method sinh cosh-method cosh,ddmatrix-method cosh
#' tanh-method tanh,ddmatrix-method tanh
#' @docType methods
#' @param x numeric distributed matrix
#' @param base a positive number: the base with respect to which logarithms are
#' computed. Defaults to e='exp(1)'.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{Reductions}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- sqrt(abs(log(x/10)))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Rounding of Numbers
#' 
#' Extensions of R rounding functions for distributed matrices.
#' 
#' Rounding to a negative number of digits means rounding to a power of ten, so
#' for example \code{round(x, digits = -2)} rounds to the nearest hundred.
#' 
#' @name Round
#' @aliases Round round-method round,ddmatrix-method round ceiling-method
#' ceiling,ddmatrix-method ceiling floor-method floor,ddmatrix-method floor
#' @docType methods
#' @param x numeric distributed matrix
#' @param digits integer indicating the number of decimal places
#' (\code{round()}) or significant digits (\code{signif()}) to be used.
#' Negative values are allowed (see 'Details').
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{MiscMath}, \link{NAs}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- ceiling(x/3)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Linear Algebra Functions
#' 
#' Linear alegbra functions for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @name LinAlg
#' @aliases LinAlg t-method t,ddmatrix-method t isSymmetric-method
#' isSymmetric,ddmatrix-method isSymmetric solve-method
#' solve,ddmatrix,ddmatrix-method solve,ddmatrix,ANY-method solve La.svd-method
#' La.svd,ddmatrix-method La.svd svd-method svd,ddmatrix-method svd
#' eigen-method eigen,ddmatrix-method eigen chol-method chol,ddmatrix-method
#' chol lu-method lu,ddmatrix-method lu
#' @docType methods
#' @param object,x,a,b numeric distributed matrices.  If applicable, \code{a}
#' and \code{b} must be on the same BLACS context and have the same blocking
#' dimension.
#' @param tol precision tolerance.
#' @param ... additional arguments.
#' @param nu number of left singular vectors to return when calculating
#' singular values.
#' @param nv number of right singular vectors to return when calculating
#' singular values.
#' @param symmetric logical, if \code{TRUE} then the matrix is assumed to be
#' symmetric and only the lower triangle is used.  Otherwise \code{x} is
#' inspected for symmetry.
#' @param only.values logical, if \code{TRUE} then only the eigenvalues are
#' returned.  Otherwise both eigenvalues and eigenvectors are returned.
#' @return \code{t()} returns the transposed matrix.
#' 
#' \code{solve()} solves systems and performs matrix inversion when argument
#' \code{b=} is missing.
#' 
#' \code{La.svd()} performs singular value decomposition, and returns the
#' transpose of right singular vectors if any are requested. Singular values
#' are stored as a global R vector. Left and right singular vectors are unique
#' up to sign. Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to
#' what the left/right singular vectors are, but the disagreement is always
#' only up to sign.
#' 
#' \code{svd()} performs singular value decomposition. Differs from
#' \code{La.svd()} in that the right singular vectors, if requested, are
#' returned non-transposed. Singular values are stored as a global R vector.
#' Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to what the
#' left/right singular vectors are, but the disagreement is always only up to
#' sign.
#' 
#' \code{eigen()} computes the eigenvalues, and eigenvectors if requested.  As
#' with \code{svd()}, eigenvalues are stored in a global R vector.
#' 
#' \code{chol()} performs Cholesky factorization.
#' 
#' \code{lu()} performs LU factorization.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(a = \"ddmatrix\")")}{} \item{list("signature(b =
#' \"ddmatrix\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{Reductions}, \link{MatMult},
#' \link{MiscMath}}
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
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- solve(t(A) %*% A)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' QR Decomposition Methods
#' 
#' \code{qr()} takes the QR decomposition.
#' 
#' \code{qr.Q()} recovers Q from the output of \code{qr()}.
#' 
#' \code{qr.R()} recovers R from the output of \code{qr()}.
#' 
#' \code{qr.qy()} multiplies \code{y} by Q.
#' 
#' \code{qr.qty()} multiplies \code{y} by the transpose of Q.
#' 
#' Functions for forming a QR decomposition and for using the outputs of these
#' numerical QR routines.
#' 
#' @name QR Decomposition
#' @aliases QR qr-method qr,ddmatrix-method qr qr.Q-method qr.Q,ANY-method qr.Q
#' qr.R-method qr.R,ANY-method qr.R qr.qy-method qr.qy,ANY-method qr.qy
#' qr.qty-method qr.qty,ANY-method qr.qty
#' @docType methods
#' @param x,y numeric distributed matrices for \code{qr()}. Otherwise, \code{x}
#' is a list, namely the return from \code{qr()}.
#' @param tol logical value, determines whether or not columns are zero
#' centered.
#' @param complete logical expression of length 1.  Indicates whether an
#' arbitrary orthogonal completion of the Q or X matrices is to be made, or
#' whether the R matrix is to be completed by binding zero-value rows beneath
#' the square upper triangle.
#' @param Dvec Not implemented for objects of class \code{ddmatrix}.  vector
#' (not matrix) of diagonal values.  Each column of the returned Q will be
#' multiplied by the corresponding diagonal value.  Defaults to all 1's.
#' @return \code{qr()} returns a list consisting of: \code{qr} - \code{rank} -
#' calculated numerical rank, \code{tau} - \code{pivot} - \code{"class"} -
#' attribute "qr".
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
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' Q <- qr.Q(qr(x))
#' print(Q)
#' 
#' finalize()
#' }
#' 
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





#' Variance, Covariance, and Correlation
#' 
#' \code{sd()} forms the vector of column standard deviations.  \code{cov()}
#' and \code{var()} form the variance-covariance matrix.  \code{cor()} forms
#' the correlation matrix.  \code{cov2cor()} scales a covariance matrix into a
#' correlation matrix.
#' 
#' \code{sd()} will compute the standard deviations of the columns, equivalent
#' to calling \code{apply(x, MARGIN=2, FUN=sd)} (which will work for
#' distributed matrices, by the way). However, this should be much faster and
#' use less memory than \code{apply()}.  If \code{reduce=FALSE} then the return
#' is a distributed matrix consisting of one (global) row; otherwise, an
#' \code{R} vector is returned, with ownership of this vector determined by
#' \code{proc.dest}.
#' 
#' \code{cov()} forms the variance-covariance matrix. Only
#' \code{method="pearson"} is implemented at this time.
#' 
#' \code{var()} is a shallow wrapper for \code{cov()} in the case of a
#' distributed matrix.
#' 
#' \code{cov2cor()} scales a covariance matrix into a correlation matrix.
#' 
#' @name Variance/Covariance
#' @aliases Variance/Covariance sd-method sd,ddmatrix-method sd,ANY-method sd
#' cov-method cov,ddmatrix-method cov var-method var,ddmatrix-method var
#' cor-method cor,ddmatrix-method cor cov2cor-method cov2cor,ddmatrix-method
#' cov2cor
#' @docType methods
#' @param x,y,V numeric distributed matrices.
#' @param na.rm logical, determines whether or not \code{NA}'s should be dealth
#' with.
#' @param reduce logical or string. See details
#' @param proc.dest Destination process (or 'all') if a reduction occurs
#' @param use character indicating how missing values should be treated.
#' Acceptable values are the same as \code{R}'s, namely "everything",
#' "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method character argument indicating which method should be used to
#' calculate covariances. Currently only "spearman" is available for
#' \code{ddmatrix}.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(V = \"ddmatrix\")")}{} }
#' @author R Core Team, Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen
#' Chen, George Ostrouchov, and Pragneshkumar Patel.
#' @seealso \code{\link{prcomp}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix("rnorm", nrow=3, ncol=3)
#' 
#' cv <- cov(x)
#' print(cv)
#' 
#' finalize()
#' }
#' 
NULL

