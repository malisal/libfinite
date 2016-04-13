/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <stdlib.h>

#include "pc.h"

static void _pc_add(bn_t *d, bn_t *a, bn_t *b, pc_group_t *pcg)
{
	bn_add(d, a, b, pcg->p);
}

static void _pc_sub(bn_t *d, bn_t *a, bn_t *b, pc_group_t *pcg)
{
	bn_sub(d, a, b, pcg->p);
}

static void _pc_mul(bn_t *d, bn_t *a, bn_t *b, pc_group_t *pcg)
{
	bn_mon_mul(d, a, b, pcg->p);
}

static void _pc_square(bn_t *d, bn_t *a, pc_group_t *pcg)
{
	_pc_mul(d, a, a, pcg);
}

static void _pc_inv(bn_t *d, bn_t *a, pc_group_t *pcg)
{
	bn_t *t = bn_copy(bn_alloc(a->n), a);
	bn_mon_inv(d, t, pcg->p);
	bn_free(t);
}

pc_point_t *pc_point_alloc(u32 n)
{
	pc_point_t *res;

	if ((res = (pc_point_t *)mem_alloc(sizeof(pc_point_t))) == NULL)
		return NULL;

	res->x = bn_alloc(n);
	res->y = bn_alloc(n);

	return res;
}

void pc_point_free(pc_point_t *p)
{
	bn_free(p->x);
	bn_free(p->y);
	mem_free(p);
}

pc_point_t *pc_point_copy(pc_point_t *d, pc_point_t *s)
{
	bn_copy(d->x, s->x);
	bn_copy(d->y, s->y);

	return d;
}

pc_point_t *pc_point_zero(pc_point_t *p)
{
	bn_zero(p->x);
	bn_zero(p->y);

	bn_set_ui(p->x, 1);

	return p;
}

int pc_point_is_zero(pc_point_t *p)
{
	return (bn_is_zero(p->x) && bn_is_zero(p->y));
}

void pc_point_to_mon(pc_point_t *p, pc_group_t *pcg)
{
	bn_to_mon(p->x, pcg->p);
	bn_to_mon(p->y, pcg->p);
}

void pc_point_from_mon(pc_point_t *p, pc_group_t *pcg)
{
	bn_from_mon(p->x, pcg->p);
	bn_from_mon(p->y, pcg->p);
}

pc_point_t *pc_point_add(pc_point_t *r, pc_point_t *p, pc_point_t *q, pc_group_t *pcg)
{
	bn_t *x1 = p->x, *y1 = p->y, *x2 = q->x, *y2 = q->y;
	bn_t *s = bn_alloc(pcg->p->n), *t = bn_alloc(pcg->p->n);
	pc_point_t *rr = pc_point_alloc(pcg->p->n);

	// r.x = x1*x2 + D*y1*y2
	_pc_mul(s, x1, x2, pcg);
	_pc_mul(t, y1, y2, pcg);
	_pc_mul(t, pcg->D, t, pcg);
	_pc_add(rr->x, s, t, pcg);

	// r.y = x1*y2 + x2*y1
	_pc_mul(s, x1, y2, pcg);
	_pc_mul(t, x2, y1, pcg);
	_pc_add(rr->y, s, t, pcg);

	pc_point_copy(r, rr);

	pc_point_free(rr);
	bn_free(t);
	bn_free(s);

	return r;
}


pc_point_t *pc_point_mul(pc_point_t *d, bn_t *a, pc_point_t *b, pc_group_t *pcg)
{
	int i;
	pc_point_t *bt = pc_point_copy(pc_point_alloc(b->x->n), b);

	pc_point_zero(d);

	for (i = 0; i <= bn_maxbit(a); i++)
	{
		if (bn_getbit(a, i))
			pc_point_add(d, d, bt, pcg);
		pc_point_add(bt, bt, bt, pcg);
	}

	pc_point_free(bt);

	return d;
}
