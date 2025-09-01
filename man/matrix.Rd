\name{RcppArmadillo-matrix-functions}
\alias{solve_cpp}
\alias{spsolve_cpp}
\alias{kronecker_cpp}
\alias{spkronecker_cpp}
\title{C++ matrix functions}
\description{
Inverse matrix and kronecker products functions coded in C++ and brought to R with RcppArmadillo for faster implementation.
}
\usage{
solve_cpp(M, symmetric)
spsolve_cpp(M)
kronecker_cpp(M, Q)
spkronecker_cpp(M, Q)
}
\arguments{
  \item{M}{a numeric matrix (invertible if used with solve functions)}
  \item{symmetric}{TRUE if \code{M} is symmetric}
  \item{Q}{a numeric matrix}
}
\value{
  \code{solve_cpp(M, symmetric)} returns a numeric matrix computed as the
  inverse of \code{M} based on Rcpp. Can use \code{symmetric} == TRUE if matrix is symmetric positive definite to use a faster implementation.

  \code{spsolve_cpp(M)} same as \code{solve_cpp}, but useful for sparse matrices.

  \code{kronecker_cpp(M, Q)} kronecker product between two compatible matrices \code{M} and \code{Q}.

  \code{spkronecker_cpp(M, Q)} same as \code{kronecker_cpp}, but useful for sparse matrices.

}
\details{
  These are matrix related function coded in C++ with the Armadillo library and brought to R with RcppArmadillo to have fast functions which are used to compute matrix inverse and kronecker products. Very useful in the implementation of this package MCMC algorithms which requires several Kronecker products and matrix inversions.
}
\references{
  See the documentation for Armadillo, reference: Conrad Sanderson and Ryan Curtin. Armadillo: An Efficient Framework for Numerical Linear Algebra. International Conference on Computer and Automation Engineering, 2025. https://doi.org/10.1109/ICCAE64891.2025.10980539.
}
\author{Rayane Dounas}
