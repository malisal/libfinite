/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _PCNR_H_
#define _PCNR_H_

#include "pc.h"

/*! Pell conic Nyberg-Rueppel context. */
typedef struct _pcnr_ctxt_t
{
   /*! Pell conic group. */
   pc_group_t *pcg;
   /*! Base field modulus. */
   bn_t *N;
   /*! Generator point G. (mon!) */
   pc_point_t *G;
   /*! Public point Q. (mon!) */
   pc_point_t *Q;
   /*! Private k. */
   bn_t *k;
} pcnr_ctxt_t;

/*! Pell conic Nyberg-Rueppel signature. */
typedef struct _pcnr_sig_t
{
   /*! Signature R. */
   bn_t *R;
   /*! Signature S. */
   bn_t *S;
} pcnr_sig_t;

/*!
* \brief Sign message using PCNR.
*/
void pcnr_sign(pcnr_ctxt_t *ctxt, pcnr_sig_t *sig, bn_t *H);

/*!
* \brief Verify message using PCNR.
*/
int pcnr_verify(pcnr_ctxt_t *ctxt, pcnr_sig_t *sig, bn_t *H);

#endif
