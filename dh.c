/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <stdlib.h>
#include <assert.h>

#include "dh.h"

dh_ctxt_t *dh_init(bn_t *p, bn_t *g)
{
	dh_ctxt_t *res;
	bn_t *t;

	if(p == NULL || g == NULL)
		return NULL;

	assert(p->n == g->n);
	
	if((res = (dh_ctxt_t *) mem_alloc(sizeof(dh_ctxt_t))) == NULL)
		return NULL;

	t = bn_alloc(p->n);

	res->p = p;

	//Check g \in [2, p - 2].
	bn_set_ui(t, 2);;
	if(bn_cmp(g, t) < 0)
		goto outerr;
	bn_sub(t, p, t, p);
	if(bn_cmp(g, t) > 0)
		goto outerr;

	res->g = bn_to_mon(g, p);

	// Generate c \in [1, p - 2].
	res->c = bn_alloc(p->n);
	while(1)
	{
		bn_rand(res->c);

		bn_zero(t);
		bn_set_ui(t, 1);
		if(bn_cmp(res->c, t) < 0) // Check c >= 1.
			continue;

		bn_set_ui(t, 2);
		bn_sub(t, p, t, p);
		if(bn_cmp(res->c, t) > 0) // Check c <= p - 2.
			continue;

		break;
	}

	// C = g^c mod p
	res->C = bn_from_mon(bn_mon_pow(bn_alloc(p->n), res->g, p, res->c), p);
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
	if(ctxt == NULL)
		return;

	bn_free(ctxt->C);
	bn_free(ctxt->c);
	mem_free(ctxt);
}

bn_t *dh_step(bn_t *K, dh_ctxt_t *ctxt, bn_t *D)
{
	if(K == NULL || ctxt == NULL || D == NULL)
		return NULL;

	assert(K->n == D->n && D->n == ctxt->p->n);

	// K = D^c mod p
	return bn_from_mon(bn_mon_pow(K, bn_to_mon(D, ctxt->p), ctxt->p, ctxt->c), ctxt->p);
}

