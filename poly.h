/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _POLY_H_
#define _POLY_H_

#include <stdio.h>

#include "bn.h"

/*! Polynomial. */
/*
* Define a polynomial P of degree n as
* P(X) = a_0 + a_1*X + a_2*X^2 + ... + a_n*X^n
* with n + 1 coefficients a_0, a_1, ..., a_n.
*/
typedef struct _poly
{
	/*! The polynomial's degree. */
	int degree;
	/*! Coefficients (mon). */
	bn_t **coeffs;
	/*! Modulus. */
	bn_t *N;
} poly_t;

/*!
* \brief Read polynomial from file stream.
*/
poly_t *poly_read(FILE *fp, poly_t *dst);

/*!
* \brief Write polynomial to file stream.
*/
poly_t *poly_write(FILE *fp, poly_t *p);

/*!
* \brief Print polynomial to file stream.
*/
void poly_print(FILE *fp, const s8 *pre, poly_t *p, const s8 *post);

/*!
* \brief Allocate polynomial (with zero coefficients).
*/
poly_t *poly_alloc(int degree, bn_t *N, int alloc_coeffs);

/*!
* \brief Free polynomial.
*/
void poly_free(poly_t *p, int free_coeffs);

/*!
* \brief Adjust polynomial degree.
*        This will alloc/free coefficients of alloc_free_coeffs is true.
*/
poly_t *poly_adjust(poly_t *p, int degree, int alloc_free_coeffs);

/*!
* \brief Copy polynomial.
*        This will zero out coefficients if the source degree is smaller and adjust is false.
*        This will truncate coefficients if the destination degree is smaller and adjust is false.
*        If adjust is true, then it will set the destination degree to that of the source and copy every coefficient.
*/
poly_t *poly_copy(poly_t *d, poly_t *s, int adjust);

/*!
* \brief Concatenate two polynomials (as lists of coefficients).
*/
poly_t *poly_concat(poly_t *d, poly_t *a, poly_t *b);

/*!
* \brief Zero out coefficients.
*/
poly_t *poly_zero(poly_t *p);

/*!
* \brief Set linear coefficient to one.
*/
poly_t *poly_one(poly_t *p);

/*!
* \brief Checks whether the polynomial is linear one.
*/
int poly_is_one(poly_t *p);

/*!
* \brief Get mathematical degree of polynomial (highest non-zero power).
*/
int poly_deg(poly_t *p);

/*!
* \brief Set coefficient.
*/
int poly_set_coeff(poly_t *p, int i, bn_t *coeff);

/*!
* \brief Free coefficient.
*/
int poly_free_coeff(poly_t *p, int i);

/*!
* \brief Fill coefficients from format.
*/
poly_t *poly_from_fmt(poly_t *p, const char *fmt, ...);

/*!
* \brief Convert coefficients to montgomery form.
*/
poly_t *poly_to_mon(poly_t *p);

/*!
* \brief Convert coefficients from montgomery form.
*/
poly_t *poly_from_mon(poly_t *p);

/*!
* \brief Evaluate polynomial at given x.
*/
bn_t *poly_eval(poly_t *p, bn_t *dst, bn_t *x);

/*!
* \brief Add two polynomials.
*        Note that d must not be the same as a or b.
*/
poly_t *poly_add_fast(poly_t *d, poly_t *a, poly_t *b);

/*!
* \brief Add two polynomials.
*/
poly_t *poly_add(poly_t *d, poly_t *a, poly_t *b);

/*!
* \brief Subtract two polynomials.
*        Note that d must not be the same as a or b.
*/
poly_t *poly_sub_fast(poly_t *d, poly_t *a, poly_t *b);

/*!
* \brief Subtract two polynomials.
*/
poly_t *poly_sub(poly_t *d, poly_t *a, poly_t *b);

/*!
* \brief Multiply two polynomials.
*        Note that d must not be the same as a or b.
*/
poly_t *poly_mul_fast(poly_t *d, poly_t *a, poly_t *b);

/*!
* \brief Multiply two polynomials.
*/
poly_t *poly_mul(poly_t *d, poly_t *a, poly_t *b);

/*!
* \brief Multiply each polynomial coefficient with a constant.
*/
poly_t *poly_mulc(poly_t *d, poly_t *a, bn_t *b);

/*!
* \brief Divide two polynomials yielding quotient and remainder.
*/
poly_t *poly_div_fast(poly_t *q, poly_t *r, poly_t *a, poly_t *b);

/*!
* \brief Compute remainder of two polynomials.
*/
poly_t *poly_rem_fast(poly_t *r, poly_t *a, poly_t *b);

#endif
