# sparseinv: An R Software package for computing the sparse inverse subset with the Takahashi equations with large datasets.
# Copyright (c) 2017 Andrew Zammit-Mangion
# Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' sparseinv
#'
#' This package creates a wrapper for the SuiteSparse routines in C that use the Takahashi equations to compute the elements of the inverse of a sparse matrix at locations where the (permuted) Cholesky factor is structurally non-zero. The resulting matrix is known as a sparse inverse subset. Some helper functions (like the permuted Cholesky factorisation) are also implemented. Support for spam matrices is currently limited and will be implemented in the future.
#' @name sparseinv-package
#' @docType package
#' @useDynLib sparseinv, .registration=TRUE
#' @import Matrix
#' @import Rcpp
#' @importFrom spam as.dgCMatrix.spam as.spam.dgCMatrix
#' @importFrom methods as is
NULL
