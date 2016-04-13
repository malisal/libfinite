/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include "ecnr.h"

void ecnr_sign(ecnr_ctxt_t *ctxt, ecnr_sig_t *sig, bn_t *H)
{
	bn_t *e = bn_alloc(ctxt->N->n),
		*kk = bn_alloc(ctxt->N->n),
		*m = bn_alloc(ctxt->N->n);
	ec_point_t *mG = ec_point_alloc(ctxt->ecg->p->n);

	// Create random(!) m.
	bn_reduce(bn_rand(m), ctxt->N);

	// R = (mG).x + e
	bn_reduce(bn_copy(e, H), ctxt->N);
	ec_point_mul(mG, m, ctxt->G, ctxt->ecg);
	ec_point_from_mon(mG, ctxt->ecg);
	bn_add(sig->R, mG->x, e, ctxt->N);

	// S = (m - kR) mod N
	bn_reduce(bn_copy(kk, ctxt->k), ctxt->N);
	bn_to_mon(kk, ctxt->N);
	bn_to_mon(sig->R, ctxt->N);
	bn_mon_mul(e, kk, sig->R, ctxt->N);
	bn_from_mon(sig->R, ctxt->N);
	bn_from_mon(e, ctxt->N);
	bn_sub(sig->S, m, e, ctxt->N);

	//Free temporaries.
	ec_point_free(mG);
	bn_free(m);
	bn_free(kk);
	bn_free(e);
}

int ecnr_verify(ecnr_ctxt_t *ctxt, ecnr_sig_t *sig, bn_t *H)
{
	int res = 0;
	bn_t *e = bn_alloc(ctxt->N->n),
		*z = bn_alloc(ctxt->N->n);
	ec_point_t *P1 = ec_point_alloc(ctxt->ecg->p->n),
		*P2 = ec_point_alloc(ctxt->ecg->p->n);

	//P1 = S*G + R*Q
	ec_point_mul(P1, sig->S, ctxt->G, ctxt->ecg);
	ec_point_mul(P2, sig->R, ctxt->Q, ctxt->ecg);
	ec_point_add(P1, P1, P2, ctxt->ecg);
	ec_point_from_mon(P1, ctxt->ecg);

	//z = R - P.x (mod N)
	bn_sub(z, sig->R, P1->x, ctxt->N);

	bn_reduce(bn_copy(e, H), ctxt->N);
	res = (bn_cmp(e, z) == BN_CMP_E);

	bn_free(z);
	bn_free(e);
	ec_point_free(P2);
	ec_point_free(P1);

	return res;
}
