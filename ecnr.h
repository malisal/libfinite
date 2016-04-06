/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _ECNR_H_
#define _ECNR_H_

#include "ec.h"

/*! Elliptic curve Nyberg-Rueppel context. */
typedef struct _ecnr_ctxt_t
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
} ecnr_ctxt_t;

/*! Elliptic curve Nyberg-Rueppel signature. */
typedef struct _ecnr_sig_t
{
   /*! Signature R. */
   bn_t *R;
   /*! Signature S. */
   bn_t *S;
} ecnr_sig_t;

/*!
* \brief Sign message using ECNR.
*/
void ecnr_sign(ecnr_ctxt_t *ctxt, ecnr_sig_t *sig, bn_t *H);

/*!
* \brief Verify message using ECNR.
*/
int ecnr_verify(ecnr_ctxt_t *ctxt, ecnr_sig_t *sig, bn_t *H);

#endif
