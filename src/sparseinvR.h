#ifndef _SPARSEINV_H_
#define _SPARSEINV_H_

#include <stddef.h>
#include <R.h>
#include <Rinternals.h>

SEXP _sparseinv_sparseinv2(SEXP nSEXP, SEXP LpSEXP, SEXP LiSEXP, SEXP LxSEXP,
                           SEXP dSEXP, SEXP UpSEXP, SEXP UjSEXP, SEXP UxSEXP,
                           SEXP ZpSEXP, SEXP ZiSEXP);
#endif
