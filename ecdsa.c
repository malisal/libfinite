/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <stdlib.h>

#include "ecdsa.h"

void ecdsa_sign(ecdsa_ctxt_t *ctxt, ecdsa_sig_t *sig, bn_t *H)
{
   bn_t *e = bn_alloc(ctxt->N->n), 
      *kk = bn_alloc(ctxt->N->n), 
      *m = bn_alloc(ctxt->N->n), 
      *minv = bn_alloc(ctxt->N->n);
   ec_point_t *mG = ec_point_alloc(ctxt->ecg->p->n);

   // Create random(!) m.
   bn_reduce(bn_rand(m), ctxt->N);

   // R = (mG).x
   ec_point_mul(mG, m, ctxt->G, ctxt->ecg);
   ec_point_from_mon(mG, ctxt->ecg);
   bn_copy(sig->R, mG->x);

   // S = m^(-1)*(e + Rk) mod N
   bn_reduce(bn_copy(e, H), ctxt->N);
   bn_reduce(bn_copy(kk, ctxt->k), ctxt->N);
   bn_to_mon(m, ctxt->N);
   bn_to_mon(e, ctxt->N);
   bn_to_mon(sig->R, ctxt->N);
   bn_to_mon(kk, ctxt->N);

   bn_mon_mul(sig->S, sig->R, kk, ctxt->N);
   bn_add(kk, sig->S, e, ctxt->N);
   bn_mon_inv(minv, m, ctxt->N);
   bn_mon_mul(sig->S, minv, kk, ctxt->N);

   bn_from_mon(sig->R, ctxt->N);
   bn_from_mon(sig->S, ctxt->N);

   //Free temporaries.
   ec_point_free(mG);
   bn_free(minv);
   bn_free(m);
   bn_free(kk);
   bn_free(e);
}

int ecdsa_verify(ecdsa_ctxt_t *ctxt, ecdsa_sig_t *sig, bn_t *H)
{
   int res = 0;
   bn_t *Sinv = bn_alloc(ctxt->N->n), 
      *e = bn_alloc(ctxt->N->n), 
      *w1 = bn_alloc(ctxt->N->n), 
      *w2 = bn_alloc(ctxt->N->n), 
      *rr = bn_alloc(ctxt->N->n);
   ec_point_t *r1 = ec_point_alloc(ctxt->ecg->p->n), 
      *r2 = ec_point_alloc(ctxt->ecg->p->n);

   bn_reduce(bn_copy(e, H), ctxt->N);
   bn_to_mon(sig->R, ctxt->N);
   bn_to_mon(sig->S, ctxt->N);
   bn_to_mon(e, ctxt->N);

   bn_mon_inv(Sinv, sig->S, ctxt->N);

   bn_mon_mul(w1, e, Sinv, ctxt->N);
   bn_mon_mul(w2, sig->R, Sinv, ctxt->N);

   bn_from_mon(w1, ctxt->N);
   bn_from_mon(w2, ctxt->N);

   ec_point_mul(r1, w1, ctxt->G, ctxt->ecg);
   ec_point_mul(r2, w2, ctxt->Q, ctxt->ecg);

   ec_point_add(r1, r1, r2, ctxt->ecg);

   ec_point_from_mon(r1, ctxt->ecg);

   bn_copy(rr, r1->x);
   bn_reduce(rr, ctxt->N);

   bn_from_mon(sig->R, ctxt->N);
   bn_from_mon(sig->S, ctxt->N);

   res = (bn_cmp(rr, sig->R) == BN_CMP_E);

   //Free temporaries.
   ec_point_free(r2);
   ec_point_free(r1);
   bn_free(rr);
   bn_free(w2);
   bn_free(w1);
   bn_free(e);
   bn_free(Sinv);

   return res;
}

void ecdsa_recover_priv(bn_t *k, ecdsa_ctxt_t *ctxt, bn_t *H1, bn_t *H2, ecdsa_sig_t *sig1, ecdsa_sig_t *sig2)
{
   bn_t *z1 = bn_alloc(ctxt->N->n), 
      *z2 = bn_alloc(ctxt->N->n),
      *t1 = bn_alloc(ctxt->N->n),
      *t2 = bn_alloc(ctxt->N->n);

   bn_reduce(bn_copy(z1, H1), ctxt->N);
   bn_reduce(bn_copy(z2, H2), ctxt->N);

   //t_1 = (S_1 - S_2)^{-1}
   bn_to_mon(bn_sub(t2, sig1->S, sig2->S, ctxt->N), ctxt->N);
   bn_mon_inv(t1, t2, ctxt->N);

   //t_1 = (Z_1 - Z_2)/(S_1 - S_2)
   bn_to_mon(bn_sub(t2, z1, z2, ctxt->N), ctxt->N);
   bn_mon_mul(t1, t2, t1, ctxt->N);

   //k = (S*t_1 - Z_1)/R_1
   bn_mon_mul(t1, bn_to_mon(sig1->S, ctxt->N), t1, ctxt->N);
   bn_from_mon(sig1->S, ctxt->N);
   bn_from_mon(t1, ctxt->N);
   bn_sub(t1, t1, z1, ctxt->N);
   bn_to_mon(t1, ctxt->N);

   bn_to_mon(sig1->R, ctxt->N);
   bn_mon_inv(t2, sig1->R, ctxt->N);
   bn_from_mon(sig1->R, ctxt->N);

   bn_mon_mul(t1, t1, t2, ctxt->N);
   bn_from_mon(t1, ctxt->N);

   bn_copy(k, t1);

   bn_free(t2);
   bn_free(t1);
   bn_free(z2);
   bn_free(z1);
}
