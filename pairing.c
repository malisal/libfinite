/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <assert.h>
#include <stdlib.h>

#include "pairing.h"

poly_t *_pairing_line(poly_t *d, ec_pqr_point_t *P, ec_pqr_point_t *R, ec_pqr_point_t *Q, ec_pqr_group_t *ecg)
{
	if (ec_pqr_point_is_zero(P) || ec_pqr_point_is_zero(R))
	{
		if (ec_pqr_point_cmp(P, R) == POINT_CMP_E)
			return poly_one(d);
		if (ec_pqr_point_is_zero(P))
			return pqr_sub(d, Q->x, R->x, ecg->p);
		if (ec_pqr_point_is_zero(R))
			return pqr_sub(d, Q->x, P->x, ecg->p);
	}
	else if (ec_pqr_point_cmp(P, R) != POINT_CMP_E)
	{
		if (poly_cmp(P->x, R->x) == POLY_CMP_E)
			return pqr_sub(d, Q->x, P->x, ecg->p);

		poly_t *l = poly_alloc(d->degree, d->N, 1);

		pqr_sub(l, R->x, P->x, ecg->p); //l = R.x - P.x
		pqr_inv(d, l, ecg->p);          //d = 1 / (R.x - P.x)
		pqr_sub(l, R->y, P->y, ecg->p); //l = R.y - P.y
		pqr_mul(l, l, d, ecg->p);       //l = (R.y - P.y) / (R.x - P.x)
		pqr_sub(d, Q->x, P->x, ecg->p); //d = Q.x - P.x
		pqr_mul(l, l, d, ecg->p);       //l = l * (Q.x - P.x)
		pqr_sub(d, Q->y, P->y, ecg->p); //d = Q.y - P.y
		pqr_sub(d, d, l, ecg->p);       //d = Q.y - P.y - l * (Q.x - P.x)

		poly_free(l, 1);

		return d;
	}

	poly_t *t = poly_alloc(d->degree, d->N, 1);
	poly_t *s = poly_alloc(d->degree, d->N, 1);

	pqr_mul(t, P->x, P->x, ecg->p); //t = P.x*P.x
	pqr_add(s, t, t, ecg->p);       //s = 2*P.x*P.x
	pqr_add(s, s, t, ecg->p);       //s = 3*P.x*P.x
	pqr_add(s, s, ecg->a, ecg->p);  //s = 3*px*px + a
	pqr_add(t, P->y, P->y, ecg->p); //t = 2*py

	if (poly_is_zero(t))
		pqr_sub(d, Q->x, P->x, ecg->p);
	else
	{
		pqr_inv(t, t, ecg->p);          //t = 1 / (2*py)
		pqr_mul(s, s, t, ecg->p);       //s = (3*px*px + a) / (2*py)
		pqr_sub(t, Q->x, P->x, ecg->p); //t = Q.x - P.x
		pqr_mul(s, s, t, ecg->p);       //s = s * (Q.x - P.x)
		pqr_sub(d, Q->y, P->y, ecg->p); //d = Q.y - P.y
		pqr_sub(d, d, s, ecg->p);       //d = Q.y - P.y - s * (Q.x - P.x)
	}

	poly_free(s, 1);
	poly_free(t, 1);

	return d;
}

poly_t *_pairing_miller(poly_t *d, ec_pqr_point_t *P, ec_pqr_point_t *Q, bn_t *n, ec_pqr_group_t *ecg)
{
	assert(!ec_pqr_point_is_zero(Q));
	assert(!bn_is_zero(n));

	int i;

	ec_pqr_point_t *R = ec_pqr_point_copy(ec_pqr_point_alloc(ecg), P);
	ec_pqr_point_t *S = ec_pqr_point_alloc(ecg);

	poly_t *l = poly_alloc(d->degree, d->N, 1);
	poly_t *v = poly_alloc(d->degree, d->N, 1);

	poly_one(d);

	for (i = bn_maxbit(n) - 1; i >= 0; i--)
	{
		ec_pqr_point_double(S, R, ecg); //S = [2]R
		_pairing_line(l, R, R, Q, ecg); //l_{R,R}
		ec_pqr_point_neg(R, S, ecg);
		_pairing_line(v, S, R, Q, ecg);
		pqr_inv(v, v, ecg->p);          //1 / v_{[2]R}
		pqr_mul(d, d, d, ecg->p);
		pqr_mul(d, d, l, ecg->p);
		pqr_mul(d, d, v, ecg->p);       //d = d^2 * l_{R,R} / v_{[2]R}

		ec_pqr_point_copy(R, S);

		if (!bn_getbit(n, i))
			continue;

		ec_pqr_point_add(S, R, P, ecg); //S = R + P
		_pairing_line(l, R, P, Q, ecg); //l_{R,P}
		ec_pqr_point_neg(R, S, ecg);
		_pairing_line(v, S, R, Q, ecg);
		pqr_inv(v, v, ecg->p);          //1 / v_{R+P}
		pqr_mul(d, d, l, ecg->p);
		pqr_mul(d, d, v, ecg->p);       //d = d * l_{R,P} / v_{R+P}

		ec_pqr_point_copy(R, S);

	}

	poly_free(v, 1);
	poly_free(l, 1);
	ec_pqr_point_free(S);
	ec_pqr_point_free(R);

	return d;
}

poly_t *pairing_weil(poly_t *d, ec_pqr_point_t *P, ec_pqr_point_t *Q, bn_t *n, ec_pqr_group_t *ecg)
{
	poly_t *t = poly_alloc(d->degree, d->N, 1);

	//d = f_{n,P}(D_Q) / f_{n,Q}(D_P)
	_pairing_miller(d, P, Q, n, ecg);
	_pairing_miller(t, Q, P, n, ecg);
	pqr_inv(t, t, ecg->p);
	pqr_mul(d, d, t, ecg->p);

	if (bn_getbit(n, 0))
	{
		//d = -d
		poly_zero(t);
		pqr_sub(d, t, d, ecg->p);
	}

	poly_free(t, 1);

	return d;
}

/*poly_t *pairing_tate(poly_t *d, ec_pqr_point_t *P, ec_pqr_point_t *Q, bn_t *n, bn_t *k, ec_pqr_group_t *ecg)
{
	bn_t *e = bn_alloc(ecg->p->N->n);

	_pairing_miller(d, P, Q, n, ecg);
	//(q^k - 1) / n
	bn_exp(e, ecg->p->N, k);
	bn_sub(e, e, one);
	bn_div(e, e, n);

	pqr_exp(d, d, e, ecg->p);

	bn_free(e);

	return d;
}*/
