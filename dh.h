/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _DH_H_
#define _DH_H_

#include "bn.h"

/*! Diffie-Hellman key exchange context. */
typedef struct _dh_ctxt
{
	/*! Prime p. */
	bn_t *p;
	/*! Primitive root modulo p \in [2, p - 2]. */
	bn_t *g;
	/*! Private random number \in [1, p - 2]. */
	bn_t *c;
	/*! Public constant to exchange. */
	bn_t *C;
} dh_ctxt_t;

/*!
* \brief Initialize Diffie–Hellman key exchange context.
* \param p Prime p.
* \param g Primitive root modulo p \in [2, p - 2] (not readonly).
* \return NULL on error.
*/
dh_ctxt_t *dh_init(bn_t *p, bn_t *g);

/*!
* \brief Free context.
*/
void dh_free(dh_ctxt_t *ctxt);

/*!
* \brief Calculate shared secret.
* \param K Shared secret dest.
* \param ctxt Diffie–Hellman key exchange context.
* \param D Exchanged constant.
* \return Shared secret K.
*/
bn_t *dh_step(bn_t *K, dh_ctxt_t *ctxt, bn_t *D);

#endif
