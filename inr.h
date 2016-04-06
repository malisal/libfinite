/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _INR_H_
#define _INR_H_

#include "bn.h"

/*! Integer Nyberg-Rueppel context. */
typedef struct _inr_ctxt_t
{
   /*! Common prime p (= u*q + 1 with small u). */
   bn_t *p;
   /*! Common prime q. */
   bn_t *q;
   /*! Common primitive root g of q. (mon!) */
   bn_t *g;
   /*! Private x \in [1, q - 1] (mod q). (mon!) */
   bn_t *x;
   /*! Public y (mod p). (mon!) */
   bn_t *y;
} inr_ctxt_t;

/*! Integer Nyber-Rueppel signature. */
typedef struct _inr_sig_t
{
   /*! Signature r (mod p) */
   bn_t *r;
   /*! Signature s (mod q) */
   bn_t *s;
} inr_sig_t;

/*!
* \brief Sign message using INR.
* \param ctxt Integer Nyberg-Rueppel context.
* \param sig Integer Nyber-Rueppel signature.
* \param m Message to sign.
*/
void inr_sign(inr_ctxt_t *ctxt, inr_sig_t *sig, bn_t *m);

/*!
* \brief Verify message sing INR.
* \param ctxt Integer Nyberg-Rueppel context.
* \param sig Integer Nyber-Rueppel signature.
* \param m Message to verify.
*/
int inr_verify(inr_ctxt_t *ctxt, inr_sig_t *sig, bn_t *m);

#endif
