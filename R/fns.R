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

#' @title Takahashi equations
#' @description Computes the sparse inverse subset of a sparse matrix \code{Q} using the Takahashi equations.
#'
#' @details This function first computes the Cholesky factor of \code{Q}. The fill-in reduction permutation is the approximate minimum degree permutation (amd) of Timothy Davis' SuiteSparse package configured to be slightly more aggressive than that in the \code{Matrix} package. The function then uses the Takahashi equations to compute the variances at the non-zero locations in the Cholesky factor from the factor itself. The equations themselves are implemented in C using the SparseSuite package of Timothy Davis.
#'
#' @param Q precision matrix of class \code{matrix}, \code{Matrix} (column-compressed, i.e., \code{dgCMatrix} or \code{dsCMatrix}), or \code{spam}
#' @param cholQp the Cholesky factor of class \code{dtCMatrix} of the permuted \code{Q} (if known already). If both \code{Q} and \code{cholQp} are specified, \code{Q} is ignored
#' @param return_perm_chol if 1, the Cholesky factor of the permuted \code{Q} is returned
#' @param P the permutation matrix of class \code{dgCMatrix} (if known already)
#' @param gc do garbage collection throughout (may increase computational time but useful for small memory machines)
#' @return if return_perm_chol == 0, the sparse inverse subset of Q is returned, where the non-zero elements correspond to those in the Cholesky factor of its permutation.
#' If !(return_perm_chol  == 0), a list with three elements is returned: \code{S} (the sparse inverse subset), Lp (the Cholesky factor of the permuted matrix) and P (the
#' permutation matrix)
#' @keywords Cholesky factor, sparse inverse subset
#' @note This package is a wrapper for C functions implemented by Timothy Davis in SuiteSparse. The author of this package has done no work on the sparse inverse routines themselves and any acknowledgment should include one to SuiteSparse (see below for reference). The author of this package was made aware of this methodology by Botond Cseke.
#' @export
#' @examples
#' require(Matrix)
#' Q = sparseMatrix(i = c(1, 1, 2, 2),
#'                  j = c(1, 2, 1, 2),
#'                  x = c(0.1, 0.2, 0.2, 1))
#' X <- cholPermute(Q)
#' S_partial = Takahashi_Davis(Q, cholQp = X$Qpermchol, P = X$P)
#' @references Takahashi, K., Fagan, J., Chin, M.-S., 1973. Formation of a sparse bus impedance matrix and its application to short circuit study. 8th PICA Conf. Proc. June 4--6, Minneapolis, Minn.
#'
#' Davis, T. A., 2014. sparseinv: Sparse Inverse Subset. URL https://au.mathworks.com/matlabcentral/fileexchange/33966-sparseinv--sparse-inverse-subset
#' Davis, T. A., 2006. Direct Methods for Sparse Linear Systems. SIAM, Philadelphia, PA.
Takahashi_Davis <- function(Q = NULL,
                            cholQp = NULL,
                            return_perm_chol = 0,
                            P = 0, gc = 0) {

    .check_args_Takahashi(Q = Q,
                          cholQp = cholQp,
                          return_perm_chol = return_perm_chol,
                          P = P,
                          gc = gc)

    n <- nrow(Q)

    if (is.null(cholQp)) {
        X <- cholPermute(Q = Q)
        L <- X$Qpermchol
        P <- X$P
    } else {
        L <- cholQp
        P <- P
    }

    if (return_perm_chol == 0) rm(cholQp)

    d <- diag(L)
    L <- tril(L %*% Diagonal(x = 1/d), -1)
    d <- d^2
    D <- Diagonal(x = d)

    #ii <- L@i + 1 # in {1,...,n}
    dp <- diff(L@p)
    jj <- rep(seq_along(dp), dp) # in {1,...,n}, non-decreasing

    if(gc) gc()
    Zpattern <- sparseMatrix(c(L@i + 1, jj, 1:n),
                             c(jj, L@i + 1, 1:n))
    rm(dp,jj)

    if(gc) gc()
    Z <- .sparseinv_wrapper(L, d, L, Zpattern)
    if (return_perm_chol == 0) {
        return(P %*% Z %*% t(P))
    } else {
        return(list(S= P %*% Z %*% t(P), Lp = cholQp, P=P))
        # Only possible for small problems
    }

}

#' @title Sparse Cholesky factorisation with fill-in reducing permutations
#' @description This function is similar to chol(A, pivot=T) when A is a sparse matrix. The fill-in reduction permutation is the approximate minimum degree permutation of
#' Davis' SuiteSparse package configured to be slightly more aggressive than that in the Matrix package.
#'
#' @param Q precision matrix of class \code{matrix}, \code{Matrix} (column-compressed, i.e., \code{dgCMatrix} or \code{dsCMatrix}), or \code{spam}
#' @return A list with two elements, Qpermchol (the permuted Cholesky factor) and P (the permutation matrix) of class Matrix. Note that \code{spam} matrices are not returned to comply with the Takahashi_Davis function which requires objects of class \code{Matrix}.
#' @keywords Cholesky factor
#' @export
#' @examples
#' require(Matrix)
#' cholPermute(sparseMatrix(i = c(1,1,2,2),
#'                          j = c(1, 2, 1, 2),
#'                          x = c(0.1, 0.2, 0.2, 1)))
#' @references Havard Rue and Leonhard Held (2005). Gaussian Markov Random Fields: Theory and Applications. Chapman & Hall/CRC Press
cholPermute <- function(Q)  {
  if(!is(Q,"Matrix") & !is(Q,"spam") & !is(Q,"matrix"))
      stop("Q needs to be a matrix, spam matrix or Matrix")

  if(!isSymmetric(Q)) stop("Q needs to be symmetric")
  n <- nrow(Q) # dimension

  ## Cast to Matrix
  if(is(Q,"matrix")) {
      Q <- as(Q,"dsCMatrix")
  } else if(is(Q, "spam"))   {
      Q <- as.dgCMatrix.spam(Q)
  }

  P <- .amd_Davis(Q)
  Qp <- Q[P,P]
  Qpermchol  <- t(Matrix::chol(Qp))
  P <- sparseMatrix(i = P, j = 1:n, x = 1)
  return(list(Qpermchol = Qpermchol,
              P = P))
}

#' @title Solve the equation Qx = y
#'
#' @description This function is similar to \code{solve(Q,y)} but with the added benefit that it allows for permuted matrices. This function does the job in order to minimise
#' user error when attempting to re-permute the matrices prior or after solving. The user also has an option for the permuted Cholesky factorisation of Q to be carried out
#' internally.
#'
#' @param Q matrix (if of class \code{Matrix} needs to be column-compressed, i.e., \code{dgCMatrix} or \code{dsCMatrix})), the Cholesky factor of which needs to be found
#' @param y matrix with the same number of rows as Q
#' @param perm if FLASE no permutation is carried out, if TRUE permuted Cholesky factors are used
#' @param cholQ the lower Cholesky factor of Q (if known already)
#' @param cholQp the lower Cholesky factor of a permuted Q (if known already)
#' @param P the permutation matrix (if known already)
#' @return x solution to Qx = y
#' @keywords Cholesky factor, linear solve
#' @export
#' @examples
#' require(Matrix)
#' Q = sparseMatrix(i = c(1, 1, 2, 2),
#'                  j = c(1, 2, 1, 2),
#'                  x = c(0.1, 0.2, 0.2, 1))
#' y = matrix(c(1, 2), 2, 1)
#' cholsolve(Q, y)
#' @references Havard Rue and Leonhard Held (2005). Gaussian Markov Random Fields: Theory and Applications. Chapman & Hall/CRC Press
cholsolve <- function(Q = NULL, y = NULL, perm = FALSE,
                      cholQ = NULL, cholQp = NULL, P = NULL)  {


    .check_args_cholsolve(Q = Q, y = y, perm = perm,
                          cholQ = cholQ, cholQp = cholQp, P = P)

    ## Solve Qx = y
    if (!perm) {
        if (is.null(cholQ)) {
            L <- t(chol(Q))
        }  else {
            L <- cholQ
        }
        v <- solve(L, y)
        x <- solve(t(L), v)
    }

    if (perm) {
        if (is.null(cholQp)) {
            QP <- cholPermute(Q)
            Lp <- QP$Qpermchol
            P <- QP$P
        } else {
            Lp <- cholQp
        }

        v <- solve(Lp, t(P) %*% y)
        w <- solve(t(Lp), v)
        x <- P %*% w
    }
    return(x)
}

#' @title Solve the equation X = AQ^{-1}t(A) under permutations
#'
#' @description This function is a wrapper of \code{solve()} for finding \code{X = AQ^{-1}t(A)} when the permuted Cholesky factor of Q is known.
#' #'
#' @param Q matrix (if of class \code{Matrix} needs to be column-compressed, i.e., \code{dgCMatrix} or \code{dsCMatrix})), the Cholesky factor of which needs to be found
#' @param A sparse or dense matrix
#' @param Lp the lower Cholesky factor of a permuted Q
#' @param P the permutation matrix
#' @return x solution to \code{X = AQ^{-1}t(A)}
#' @keywords Cholesky factor, linear solve
#' @export
#' @examples
#' require(Matrix)
#' Q = sparseMatrix(i = c(1, 1, 2, 2),
#'                  j = c(1, 2, 1, 2),
#'                  x = c(0.1, 0.2, 0.2, 1))
#' X <- cholPermute(Q)
#' y <- matrix(c(1,2), 2, 1)
#' A <- y %*% t(y)
#' cholsolveAQinvAT(Q,A,X$Qpermchol,X$P)
#' @references Havard Rue and Leonhard Held (2005). Gaussian Markov Random Fields: Theory and Applications. Chapman & Hall/CRC Press
cholsolveAQinvAT <- function(Q = NULL, A = NULL, Lp = NULL, P = NULL) {

    .check_args_cholsolveAQinvAT(Q = Q, A = A, Lp = Lp, P = P)

    if(is.null(Lp) & is.null(P)) {
        QP <- cholPermute(Q)
        Lp <- QP$Qpermchol
        P <- QP$P
    }

    #Solve X = AQ^{-1}t(A)
    W <- t(solve(Lp, t(A %*% P)))
    tcrossprod(W)
}

#' @title Return the symbolic representation of a Matrix
#'
#' @description This function takes an object of class Matrix and returns the same Matrix with all elements replaced with 1
#' #'
#' @param A object of class Matrix
#' @return object of class Matrix
#' @export
#' @examples
#' require(Matrix)
#' Q = sparseMatrix(i = c(1, 1, 2, 2),
#'                  j = c(1, 2, 1, 2),
#'                  x = c(0.1, 0.2, 0.2, 1))
#' Qsymb <- symb(Q)
#' Qsymb
symb <- function(A) {
  A@x <- rep(1, length(A@x))
  A
}

#' @title Densify with explicit zeroes
#'
#' @description This function takes two sparse matrices and returns the first matrix padded with explicit zeros so that it is at least dense (probably denser) than the second matrix. This function only works with matrices of class Matrix
#' #'
#' @param A object of class Matrix
#' @param B object of class Matrix
#' @return object of class Matrix
#' @export
#' @examples
#' require(Matrix)
#' Q1 <- sparseMatrix(i = c(1, 2, 2), j = c(1, 1, 2), x = c(0.1, 0.2, 1))
#' Q2 <- sparseMatrix(i = c(1, 1, 2, 2),j = c(1, 2, 1, 2), x = c(0.1, 0.3, 0.2, 1))
#' Q1dens <- densify(Q1, Q2)
#' Q1
#' Q1dens
densify <- function(A, B) {
  ## Makes A at least as dense as B
  As <- symb(A)
  Bs <- symb(B)
  delta <- as(As - Bs,"dgTMatrix")
  idx <- which(delta@x == -1)
  addon <- sparseMatrix(delta@i + 1,
                        delta@j + 1,
                        x = 0)
  A + addon
}

######## NOT EXPORTED #################

.sparseinv_wrapper <- function(L,d,U,Zpattern) {

    n <- nrow(L)
    Lp <- L@p
    Li <- L@i
    Lx <- L@x

    Up <- U@p
    Uj <- U@i
    Ux <- U@x

    Zpatp <- Zpattern@p
    Zpati <- Zpattern@i
    znz = Zpatp [n+1]


    if(0) {
        X <- .C("sparseinv",as.integer(n),as.integer(Lp),as.integer(Li),as.double(Lx),as.double(d),as.integer(Up),as.integer(Uj),as.double(Ux),as.integer(Zpatp),as.integer(Zpati),result = double(znz))
        X <- X$result
    } else {
        X <- sparseinv2(n = as.integer(n),Lp = as.integer(Lp), Li = as.integer(Li),Lx = as.double(Lx),d = as.double(d),Up = as.integer(Up),Uj = as.integer(Uj),Ux = as.double(Ux),Zp = as.integer(Zpatp),Zi = as.integer(Zpati))
    }

    rm(U,L,Zpattern,Ux,Uj,Up,Lp,Li,Lx)
    Z <- sparseMatrix(p = Zpatp, i =Zpati, x = X,index1=F)

    return(Z)
}

.amd_Davis <- function(Q) {
    n <- nrow(Q)
    Ap <- Q@p
    Ai <- Q@i

    X <- .C("AMD_order_wrapper",as.integer(n),as.integer(Ap),as.integer(Ai),
            P = integer(n), Control=double(5),Info=double(20))
    return(X$P + 1)
}


.amd_test <- function() {
    n=24
    Ap = c( 0, 9, 15, 21, 27, 33, 39, 48, 57, 61, 70, 76, 82, 88, 94, 100,
            106, 110, 119, 128, 137, 143, 152, 156, 160 )

    Ai = c(0, 5, 6, 12, 13, 17, 18, 19, 21,
           1, 8, 9, 13, 14, 17,
           2, 6, 11, 20, 21, 22,
           3, 7, 10, 15, 18, 19,
           4, 7, 9, 14, 15, 16,
           0, 5, 6, 12, 13, 17,
           0, 2, 5, 6, 11, 12, 19, 21, 23,
           3, 4, 7, 9, 14, 15, 16, 17, 18,
           1, 8, 9, 14,
           1, 4, 7, 8, 9, 13, 14, 17, 18,
           3, 10, 18, 19, 20, 21,
           2, 6, 11, 12, 21, 23,
           0, 5, 6, 11, 12, 23,
           0, 1, 5, 9, 13, 17,
           1, 4, 7, 8, 9, 14,
           3, 4, 7, 15, 16, 18,
           4, 7, 15, 16,
           0, 1, 5, 7, 9, 13, 17, 18, 19,
           0, 3, 7, 9, 10, 15, 17, 18, 19,
           0, 3, 6, 10, 17, 18, 19, 20, 21,
           2, 10, 19, 20, 21, 22,
           0, 2, 6, 10, 11, 19, 20, 21, 22,
           2, 20, 21, 22,
           6, 11, 12, 23 )
    #Q <- as(sparseMatrix(i=Ai,p=Ap,index1=F,x=1),"dgTMatrix")
    X <- .C("AMD_order_wrapper", as.integer(n), as.integer(Ap), as.integer(Ai),
            P = integer(n), Control = double(5), Info = double(20))
}

.check_args_Takahashi <- function(Q, cholQp, P, return_perm_chol, gc) {
    if(!is.null(Q))
        if(!(is(Q,"spam") | is(Q,"Matrix") | is(Q,"matrix")))
            stop("Q needs to be empty, of class matrix, Matrix, or spam")

    if(!is.null(cholQp))
        if(!(is(Q,"spam") | is(Q,"Matrix") | is(Q,"matrix")))
            stop("Q needs to be empty, of class matrix, Matrix, or spam")

    if(!is.null(P))
        if(!(is(Q,"spam") | is(Q,"Matrix") | is(Q,"matrix")))
            stop("Q needs to be empty, of class matrix, Matrix, or spam")

    if(is.null(Q) & is.null(cholQp))
        stop("At least one of Q and cholQp needs to be specified")

    if(!(return_perm_chol == 1 | return_perm_chol == 0)) {
        stop("return_perm_chol needs to be 1 or 0 (TRUE or FALSE)")
    }

    if(!(gc == 1 | gc == 0)) {
        stop("gc needs to be 1 or 0 (TRUE or FALSE)")
    }

}

.check_args_cholsolve <- function(Q = NULL, y = NULL, perm = FALSE,
                      cholQ = NULL, cholQp = NULL, P = NULL) {

    if(is.null(Q) & is.null(cholQ) & is.null(cholQp))
        stop("One of Q, cholQ or cholQp has to be specified")

    if(!is.null(cholQ) & !is.null(cholQp))
        stop("cholQ or cholQp need to be specified, not both")

    if(!is.null(cholQp) & is.null(P))
        stop("When specifying cholQp, P also needs to be set")

    if(!(perm == 1 | perm == 0))
        stop("perm needs to be 1 or 0 (TRUE or FALSE)")

    if(is.null(y))
        stop("y needs to be set")

    if(!perm & (is.null(Q) & is.null(cholQ)))
        stop("If perm == FALSE, then Q or cholQ need to be set")

    if(perm & (is.null(Q) & is.null(cholQp)))
        stop("If perm == TRUE, then Q or cholQp need to be set")

}

.check_args_cholsolveAQinvAT <- function(Q = NULL, A = NULL, Lp = NULL, P = NULL)  {
    if(is.null(A))
        stop("A needs to be specified")
    if(is.null(Q) & is.null(Lp) & is.null(P))
        stop("Q or (Lp and P) need to be specified")
    if(!is.null(Lp) & is.null(P))
        stop("Both Lp and P need to be specified")
    if(is.null(Lp) & !is.null(P))
        stop("Both Lp and P need to be specified")
}
