[![Build Status](https://travis-ci.org/andrewzm/sparseinv.svg)](https://travis-ci.org/andrewzm/sparseinv)

# sparseinv

`sparseinv` is an R package that wraps SuiteSparse routines that compute the Takahashi equations. These equations compute the elements of the inverse of a sparse matrix at locations where the (permuted) Cholesky factor is non-zero. The resulting matrix is known as a sparse inverse subset. Some helper functions (like the permuted Cholesky factorisation)   are also implemented. Support for `spam` matrices is currently limited and will be implemented in the future. 