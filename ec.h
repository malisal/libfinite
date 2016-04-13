/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _EC_H_
#define _EC_H_

#include "bn.h"

/*! Elliptic curve point. */
typedef struct _ec_point
{
	/*! x coord. */
	bn_t *x;
	/*! y coord. */
	bn_t *y;
} ec_point_t;

/*! Elliptic curve group parameters (defining equation: y^2 = x^3 + ax + b). */
typedef struct _ec_group
{
	/*! Modulus. */
	bn_t *p;
	/*! Parameter a. (mon!) */
	bn_t *a;
	/*! Parameter b. (mon!) */
	bn_t *b;
} ec_group_t;

/*!
* \brief Allocate point.
*/
ec_point_t *ec_point_alloc(uint32_t n);

/*!
* \brief Free point.
*/
void ec_point_free(ec_point_t *p);

/*!
* \brief Copy point.
*/
ec_point_t *ec_point_copy(ec_point_t *d, ec_point_t *s);

/*!
* \brief Zero out point.
*/
ec_point_t *ec_point_zero(ec_point_t *p);

/*!
* \brief Test for zero point.
*/
int ec_point_is_zero(ec_point_t *p);

/*!
* \brief Convert point to montgomery form.
*/
void ec_point_to_mon(ec_point_t *p, ec_group_t *ecg);

/*!
* \brief Convert point from montgomery form.
*/
void ec_point_from_mon(ec_point_t *p, ec_group_t *ecg);

/*!
* \brief Double point.
*/
ec_point_t *ec_point_double(ec_point_t *r, ec_point_t *p, ec_group_t *ecg);

/*!
* \brief Add two points.
*/
ec_point_t *ec_point_add(ec_point_t *r, ec_point_t *p, ec_point_t *q, ec_group_t *ecg);

/*!
* \brief Multiply point with bignum (d = a * b).
*/
ec_point_t *ec_point_mul(ec_point_t *d, bn_t *a, ec_point_t *b, ec_group_t *ecg);

#endif // _EC_H_

