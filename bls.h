/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _BLS_H_
#define _BLS_H_

#include "ec_pqr.h"

typedef struct _bls_ctxt_t
{
	/*! Base field prime. */
	bn_t *p;
	/*! Elliptic curve group. */
	ec_pqr_group_t *ecg;
	/*! Generator point G1 (in E/GF(p)). */
	ec_pqr_point_t *G1;
	/*! Generator point G2 (in E/GF(p^12)). */
	ec_pqr_point_t *G2;
	/*! Group order r (for both points). */
	bn_t *r;
	/*! Public point P. */
	ec_pqr_point_t *P;
	/*! Private k. */
	bn_t *k;
} bls_ctxt_t;

void bls_sign(bls_ctxt_t *ctxt, ec_pqr_point_t *sig, bn_t *m);
int bls_verify(bls_ctxt_t *ctxt, ec_pqr_point_t *sig, bn_t *m);

#endif
