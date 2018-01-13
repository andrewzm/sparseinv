pkgname <- "sparseinv"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "sparseinv-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('sparseinv')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Takahashi_Davis")
### * Takahashi_Davis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Takahashi_Davis
### Title: Takahashi equations
### Aliases: Takahashi_Davis
### Keywords: Cholesky factor, inverse sparse subset

### ** Examples

require(Matrix)
Q = sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
X <- cholPermute(Q)
S_partial = Takahashi_Davis(Q,cholQp = X$Qpermchol,P=X$P)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Takahashi_Davis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cholPermute")
### * cholPermute

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cholPermute
### Title: Sparse Cholesky Factorisation with fill-in reducing permutations
### Aliases: cholPermute
### Keywords: Cholesky factor

### ** Examples

require(Matrix)
cholPermute(sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1)))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cholPermute", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cholsolve")
### * cholsolve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cholsolve
### Title: Solve the equation Qx = y
### Aliases: cholsolve
### Keywords: Cholesky factor, linear solve

### ** Examples

require(Matrix)
Q = sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
y = matrix(c(1,2),2,1)
cholsolve(Q,y)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cholsolve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cholsolveAQinvAT")
### * cholsolveAQinvAT

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cholsolveAQinvAT
### Title: Solve the equation X = AQ^-1t(A) under permutations
### Aliases: cholsolveAQinvAT
### Keywords: Cholesky factor, linear solve

### ** Examples

require(Matrix)
Q <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
X <- cholPermute(Q)
y <- matrix(c(1,2),2,1)
A <- y %*% t(y)
cholsolveAQinvAT(Q,A,X$Qpermchol,X$P)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cholsolveAQinvAT", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("densify")
### * densify

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: densify
### Title: Densify with explicit zeroes
### Aliases: densify

### ** Examples

require(Matrix)
Q1 <- sparseMatrix(i=c(1,2,2),j=c(1,1,2),x=c(0.1,0.2,1))
Q2 <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.3,0.2,1))
Q1dens <- densify(Q1,Q2)
Q1
Q1dens



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("densify", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("symb")
### * symb

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: symb
### Title: Return the symbolic representation of a Matrix
### Aliases: symb

### ** Examples

require(Matrix)
Q <- sparseMatrix(i=c(1,2,2),j=c(1,1,2),x=c(0.1,0.2,1))
Qsymb <- symb(Q)
Qsymb



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("symb", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
