/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#if !defined(NAKED)
#include <stdlib.h>
#include <assert.h>
#endif

#include "dh.h"

dh_ctxt_t *dh_init(bn_t *p, bn_t *g)
{
	dh_ctxt_t *res;
	bn_t *t;

	if (p == NULL || g == NULL)
		return NULL;

	assert(p->n == g->n);

	if ((res = (dh_ctxt_t *)mem_alloc(sizeof(dh_ctxt_t))) == NULL)
		return NULL;

	res->g = g;
	res->p = p;

	t = bn_copy(bn_alloc(p->n), p);
	bn_sub_ui(t, t, 2, p);

	// Check g \in [2, p - 2].
	if (bn_cmp_ui(g, 2) < 0 || bn_cmp(g, t) > 0)
		goto outerr;

	// Generate c \in [1, p - 2].
	res->c = bn_alloc(p->n);
	bn_rand_range(res->c, 1, p, 2);

	// C = g^c mod p
	res->C = bn_alloc(p->n);
	bn_pow_mod(res->C, res->g, res->c, p);
	goto outok;

outerr:;
	mem_free(res);
	res = NULL;

outok:;
	bn_free(t);

	return res;
}

void dh_free(dh_ctxt_t *ctxt)
{
	if (ctxt == NULL)
		return;

	bn_free(ctxt->C);
	bn_free(ctxt->c);
	mem_free(ctxt);
}

bn_t *dh_step(bn_t *K, dh_ctxt_t *ctxt, bn_t *D)
{
	if (K == NULL || ctxt == NULL || D == NULL)
		return NULL;

	assert(K->n == D->n && D->n == ctxt->p->n);

	// K = D^c mod p
	return bn_pow_mod(K, D, ctxt->c, ctxt->p);
}
