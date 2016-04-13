/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <assert.h>
#include <stdlib.h>
#include "pqr.h"

poly_t *pqr_add_fast(poly_t *d, poly_t *a, poly_t *b, poly_t *N)
{
	int i;

	assert(d->degree == a->degree && d->degree == b->degree);

	for (i = 0; i < d->degree; i++)
		bn_add(d->coeffs[i], a->coeffs[i], b->coeffs[i], N->N);

	return d;
}

poly_t *pqr_sub_fast(poly_t *d, poly_t *a, poly_t *b, poly_t *N)
{
	int i;

	assert(d->degree == a->degree && d->degree == b->degree);

	for (i = 0; i < d->degree; i++)
		bn_sub(d->coeffs[i], a->coeffs[i], b->coeffs[i], N->N);

	return d;
}

poly_t *pqr_mul_fast(poly_t *d, poly_t *a, poly_t *b, poly_t *N)
{
	poly_mul_fast(d, a, b);
	return poly_rem_fast(d, d, N);
}

poly_t *pqr_inv_fast(poly_t *d, poly_t *p, poly_t *N)
{
	assert(d->degree == p->degree);

	poly_t *a = poly_alloc(p->degree, d->N, 1);
	poly_copy(a, p, 0);
	poly_t *b = poly_alloc(N->degree, d->N, 1);
	poly_copy(b, N, 0);

	poly_t *ea = poly_alloc(d->degree, d->N, 1);
	poly_t *eb = poly_alloc(d->degree, d->N, 1);

	poly_t *t1 = poly_alloc(d->degree, d->N, 1);
	poly_t *t2 = poly_alloc(d->degree, d->N, 1);

	poly_one(ea);  // ea = 1
	poly_zero(eb); // eb = 0

	while (1)
	{
		poly_div_fast(t2, a, a, b);   // q,a = quo_rem(a,b)
		poly_mul_fast(t1, t2, eb);    // t1 = q*eb
		poly_copy(t2, ea, 1);
		poly_sub_fast(ea, t2, t1);    // ea = ea - t1

		if (a->degree == 0)
			break;

		// Swap a and b.
		poly_t *x = b;
		b = a;
		a = x;

		// Swap ea and eb.
		x = eb;
		eb = ea;
		ea = x;
	}

	poly_copy(d, ea, 0);

	// Make it monic.
	if (!poly_is_one(a))
	{
		bn_mon_inv(a->coeffs[0], a->coeffs[0], a->N);
		poly_mulc(d, d, a->coeffs[0]);
	}

	poly_free(t2, 1);
	poly_free(t1, 1);
	poly_free(eb, 1);
	poly_free(ea, 1);
	poly_free(b, 1);
	poly_free(a, 1);

	return d;
}

poly_t *pqr_exp_fast(poly_t *d, poly_t *p, bn_t *e, poly_t *N)
{
	int i;

	assert(d->degree == p->degree);

	poly_one(d);

	poly_t *r = poly_alloc(p->degree, d->N, 1);

	for (i = bn_maxbit(e); i >= 0; i--)
	{
		pqr_mul_fast(r, d, d, N);
		if (bn_getbit(e, i))
			pqr_mul_fast(d, r, p, N);
		else
			poly_copy(d, r, 0);
	}

	poly_free(r, 1);

	return d;
}
