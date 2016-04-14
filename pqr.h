/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _PQR_H_
#define _PQR_H_

#include "poly.h"

/*!
* \brief Add two polynomials modulo N.
*/
poly_t *pqr_add(poly_t *d, poly_t *a, poly_t *b, poly_t *N);

/*!
* \brief Subtract two polynomials modulo N.
*/
poly_t *pqr_sub(poly_t *d, poly_t *a, poly_t *b, poly_t *N);

/*!
* \brief Multiply two polynomials modulo N.
*        Note that d must not be the same as a or b.
*/
poly_t *pqr_mul_fast(poly_t *d, poly_t *a, poly_t *b, poly_t *N);

/*!
* \brief Multiply two polynomials modulo N.
*/
poly_t *pqr_mul(poly_t *d, poly_t *a, poly_t *b, poly_t *N);

/*!
* \brief Invert two polynomials modulo N.
*/
poly_t *pqr_inv(poly_t *d, poly_t *p, poly_t *N);

/*!
* \brief Exponentiate two polynomials modulo N.
*        Note that d must not be the same as p.
*/
poly_t *pqr_exp_fast(poly_t *d, poly_t *p, bn_t *e, poly_t *N);

/*!
* \brief Exponentiate two polynomials modulo N.
*/
poly_t *pqr_exp(poly_t *d, poly_t *p, bn_t *e, poly_t *N);

#endif
