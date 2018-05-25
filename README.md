[![Build Status](https://travis-ci.org/andrewzm/sparseinv.svg)](https://travis-ci.org/andrewzm/sparseinv)
[![codecov.io](http://codecov.io/github/andrewzm/sparseinv/coverage.svg?branch=master)](http://codecov.io/github/andrewzm/sparseinv?branch=master)

# sparseinv

`sparseinv` is an R package that wraps sparseinv routines coded by Timothy Davis that compute the Takahashi equations. These equations compute the elements of the inverse of a sparse matrix at locations where the (permuted) Cholesky factor is non-zero. The resulting matrix is known as a sparse inverse subset. Some helper functions (like the permuted Cholesky factorisation)   are also implemented. Support for `spam` matrices is currently limited and will be implemented in the future. 

Davis T (2014). sparseinv: Sparse Inverse Subset. URL https://au.mathworks.com/matlabcentral/fileexchange/33966-sparseinv--sparse-inverse-subset.
