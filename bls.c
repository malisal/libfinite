/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include "bls.h"
#include "pairing.h"

/*
* BLS:
*   pairing e: G1 x G2 -> GT (all of order r)
*   private parameters:
*     x \in [0, r - 1] (random)
*   public parameters:
*     G (a generator point of G2)
*     P = [x]G
*   signing a message m:
*     h = H(m) (e.g. H(m) = [m]Q where Q is a generator point of G1)
*     S = [x]h (the signature)
*   verification of a signature S:
*     e(S,G) == e(H(m),P)
*/

void bls_sign(bls_ctxt_t *ctxt, ec_pqr_point_t *sig, bn_t *m)
{
	//S = [k]H(m)
	ec_pqr_point_mul(sig, m, ctxt->G1, ctxt->ecg);
	ec_pqr_point_mul(sig, ctxt->k, sig, ctxt->ecg);
}

int bls_verify(bls_ctxt_t *ctxt, ec_pqr_point_t *sig, bn_t *m)
{
	int res = 0;

	poly_t *a = poly_alloc(ctxt->ecg->p->degree - 1, ctxt->p, 1);
	poly_t *b = poly_alloc(ctxt->ecg->p->degree - 1, ctxt->p, 1);
	ec_pqr_point_t *P = ec_pqr_point_alloc(ctxt->ecg);

	//e(H(m), P)
	ec_pqr_point_mul(P, m, ctxt->G1, ctxt->ecg);
	pairing_weil(a, P, ctxt->P, ctxt->r, ctxt->ecg);
	//e(S, G)
	pairing_weil(b, sig, ctxt->G2, ctxt->r, ctxt->ecg);

	res = (poly_cmp(a, b) == POLY_CMP_E ? 1 : 0);

	ec_pqr_point_free(P);
	poly_free(b, 1);
	poly_free(a, 1);

	return res;
}
