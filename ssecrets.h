/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _SSECRETS_H_
#define _SSECRETS_H_

#include "bn.h"
#include "poly.h"

/*!
* \brief Create Shamir's secret sharing polynomial.
*/
poly_t *ssecrets_create_poly(bn_t *s, u32 k, bn_t *N);

/*!
* \brief Create share.
*/
bn_t *ssecrets_create_share(poly_t *p, bn_t *x);

/*!
* \brief Calculate secret from shares.
*/
bn_t *ssecrets_calc_secret(bn_t **x, bn_t **s, u32 cnt, bn_t *N);

#endif
