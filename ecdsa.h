/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _ECDSA_H_
#define _ECDSA_H_

#include "ec.h"

/*! ECDSA context. */
typedef struct _ecdsa_ctxt_t
{
   /*! Elliptic curve group. */
   ec_group_t *ecg;
   /*! Base field modulus. */
   bn_t *N;
   /*! Generator point G. (mon!) */
   ec_point_t *G;
   /*! Public point Q. (mon!) */
   ec_point_t *Q;
   /*! Private k. */
   bn_t *k;
} ecdsa_ctxt_t;

/*! ECDSA signature. */
typedef struct _ecdsa_sig_t
{
   /*! Signature R. */
   bn_t *R;
   /*! Signature S. */
   bn_t *S;
} ecdsa_sig_t;

/*!
* \brief Sign message using ECDSA.
*/
void ecdsa_sign(ecdsa_ctxt_t *ctxt, ecdsa_sig_t *sig, bn_t *H);

/*!
* \brief Verify message using ECDSA.
*/
int ecdsa_verify(ecdsa_ctxt_t *ctxt, ecdsa_sig_t *sig, bn_t *H);

/*!
* \brief Recover private key from flawed signatures (same R values).
*/
void ecdsa_recover_priv(bn_t *k, ecdsa_ctxt_t *ctxt, bn_t *H1, bn_t *H2, ecdsa_sig_t *sig1, ecdsa_sig_t *sig2);

#endif
