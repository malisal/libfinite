/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _PC_H_
#define _PC_H_

#include "bn.h"

/*! Pell conic point. */
typedef struct _pc_point
{
	/*! x coord. */
	bn_t *x;
	/*! y coord. */
	bn_t *y;
} pc_point_t;

/*! Pell conic group parameters (defining equation: x^2 - D*y^2 = 1). */
typedef struct _pc_group
{
	/*! Modulus. */
	bn_t *p;
	/*! Parameter D. (mon!) */
	bn_t *D;
} pc_group_t;

/*!
* \brief Allocate point.
*/
pc_point_t *pc_point_alloc(u32 n);

/*!
* \brief Free point.
*/
void pc_point_free(pc_point_t *p);

/*!
* \brief Copy point.
*/
pc_point_t *pc_point_copy(pc_point_t *d, pc_point_t *s);

/*!
* \brief Zero out point.
*/
pc_point_t *pc_point_zero(pc_point_t *p);

/*!
* \brief Test for zero point.
*/
int pc_point_is_zero(pc_point_t *p);

/*!
* \brief Convert point to montgomery form.
*/
void pc_point_to_mon(pc_point_t *p, pc_group_t *pcg);

/*!
* \brief Convert point from montgomery form.
*/
void pc_point_from_mon(pc_point_t *p, pc_group_t *pcg);

/*!
* \brief Add two points.
*/
pc_point_t *pc_point_add(pc_point_t *r, pc_point_t *p, pc_point_t *q, pc_group_t *pcg);

/*!
* \brief Multiply point with bignum (d = a * b).
*/
pc_point_t *pc_point_mul(pc_point_t *d, bn_t *a, pc_point_t *b, pc_group_t *pcg);

#endif
