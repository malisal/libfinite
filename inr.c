/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include "inr.h"

void inr_sign(inr_ctxt_t *ctxt, inr_sig_t *sig, bn_t *m)
{
   bn_t *k = bn_alloc(ctxt->p->n),
      *r2 = bn_alloc(ctxt->p->n),
      *t = bn_alloc(ctxt->p->n);

   // Try to create a valid signature.
   while(1)
   {
      bn_reduce(bn_rand(k), ctxt->q); // TODO: should be k \in [1,q-1].
      bn_to_mon(bn_reduce(bn_copy(sig->r, m), ctxt->p), ctxt->p); // r=m (mod p)
      bn_mon_pow(t, ctxt->g, ctxt->p, k);                         // t=g^k
      bn_mon_mul(sig->r, sig->r, t, ctxt->p);                     // r=m*g^k
      bn_from_mon(sig->r, ctxt->p);
      bn_reduce(bn_copy(r2, sig->r), ctxt->q);                    // r2 = r mod q
      bn_to_mon(r2, ctxt->q);
      bn_mon_mul(r2, r2, ctxt->x, ctxt->q);                       // r2 = r2*x
      bn_from_mon(r2, ctxt->q);
      bn_zero(sig->s);
      bn_sub(sig->s, sig->s, k, ctxt->q);                         // s = -k
      bn_sub(sig->s, sig->s, r2, ctxt->q);                        // s = -k - x*r2
      if(inr_verify(ctxt, sig, m))
         break;
   }

   bn_free(t);
   bn_free(r2);
   bn_free(k);
}

int inr_verify(inr_ctxt_t *ctxt, inr_sig_t *sig, bn_t *m)
{
   bn_t *r2 = bn_alloc(ctxt->p->n), 
      *m2 = bn_alloc(ctxt->p->n),
      *t = bn_alloc(ctxt->p->n);

   bn_to_mon(bn_copy(m2, sig->r), ctxt->p);
   bn_reduce(bn_copy(r2, sig->r), ctxt->q);  // r2 = r mod q
   bn_to_mon(sig->r, ctxt->p);

   bn_mon_pow(t, ctxt->g, ctxt->p, sig->s);  // t=g^s
   bn_mon_mul(m2, m2, t, ctxt->p);           // m2=r*g^s
   bn_mon_pow(t, ctxt->y, ctxt->p, r2);      // t=y^r2
   bn_mon_mul(m2, m2, t, ctxt->p);           // m2=r*g^s*y^r2

   bn_from_mon(sig->r, ctxt->p);
   bn_from_mon(m2, ctxt->p);

   int res = bn_cmp(bn_reduce(bn_copy(t, m), ctxt->p), m2) == BN_CMP_E;

   bn_free(t);
   bn_free(m2);
   bn_free(r2);

   return res;
}
