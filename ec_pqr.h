/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _EC_PQR_H_
#define _EC_PQR_H_

#include "pqr.h"

/*! Elliptic curve (over polynomial quotient ring) point. */
typedef struct _ec_pqr_point
{
	/*! x coord. */
	poly_t *x;
	/*! y coord. */
	poly_t *y;
} ec_pqr_point_t;

/*! Elliptic curve (over polynomial quotient ring) group parameters (defining equation: y^2 = x^3 + ax + b). */
typedef struct _ec_pqr_group
{
	/*! Modulus. */
	poly_t *p;
	/*! Parameter a. (mon!) */
	poly_t *a;
	/*! Parameter b. (mon!) */
	poly_t *b;
} ec_pqr_group_t;

/*!
* \brief Allocate point.
*/
ec_pqr_point_t *ec_pqr_point_alloc(ec_pqr_group_t *ecg);

/*!
* \brief Free point.
*/
void ec_pqr_point_free(ec_pqr_point_t *p);

/*!
* \brief Copy point.
*/
ec_pqr_point_t *ec_pqr_point_copy(ec_pqr_point_t *d, ec_pqr_point_t *s);

/*!
* \brief Zero out point.
*/
ec_pqr_point_t *ec_pqr_point_zero(ec_pqr_point_t *p);

/*!
* \brief Test for zero point.
*/
int ec_pqr_point_is_zero(ec_pqr_point_t *p);

/*!
* \brief Convert point to montgomery form.
*/
void ec_pqr_point_to_mon(ec_pqr_point_t *p);

/*!
* \brief Convert point from montgomery form.
*/
void ec_pqr_point_from_mon(ec_pqr_point_t *p);

/*!
* \brief Double point.
*/
ec_pqr_point_t *ec_pqr_point_double(ec_pqr_point_t *r, ec_pqr_point_t *p, ec_pqr_group_t *ecg);

/*!
* \brief Add two points.
*/
ec_pqr_point_t *ec_pqr_point_add(ec_pqr_point_t *r, ec_pqr_point_t *p, ec_pqr_point_t *q, ec_pqr_group_t *ecg);

/*!
* \brief Multiply point with bignum (d = a * b).
*/
ec_pqr_point_t *ec_pqr_point_mul(ec_pqr_point_t *d, bn_t *a, ec_pqr_point_t *b, ec_pqr_group_t *ecg);

#endif // _EC_PQR_H_
