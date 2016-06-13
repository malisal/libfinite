/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <stdlib.h>
#include "ec_pqr.h"

static void _ec_pqr_add(poly_t *d, poly_t *a, poly_t *b, ec_pqr_group_t *ecg)
{
	pqr_add(d, a, b, ecg->p);
}

static void _ec_pqr_sub(poly_t *d, poly_t *a, poly_t *b, ec_pqr_group_t *ecg)
{
	pqr_sub(d, a, b, ecg->p);
}

static void _ec_pqr_mul(poly_t *d, poly_t *a, poly_t *b, ec_pqr_group_t *ecg)
{
	pqr_mul(d, a, b, ecg->p);
}

static void _ec_pqr_square(poly_t *d, poly_t *a, ec_pqr_group_t *ecg)
{
	_ec_pqr_mul(d, a, a, ecg);
}

static void _ec_pqr_inv(poly_t *d, poly_t *a, ec_pqr_group_t *ecg)
{
	pqr_inv(d, a, ecg->p);
}

ec_pqr_point_t *ec_pqr_point_alloc(ec_pqr_group_t *ecg)
{
	ec_pqr_point_t *res;

	if ((res = (ec_pqr_point_t *)mem_alloc(sizeof(ec_pqr_point_t))) == NULL)
		return NULL;

	//(degree - 1) since x, y are in the polynomial quotient ring.
	res->x = poly_alloc(ecg->p->degree - 1, ecg->p->N, 1);
	res->y = poly_alloc(ecg->p->degree - 1, ecg->p->N, 1);

	return res;
}

void ec_pqr_point_free(ec_pqr_point_t *p)
{
	poly_free(p->x, 1);
	poly_free(p->y, 1);
	mem_free(p);
}

ec_pqr_point_t *ec_pqr_point_copy(ec_pqr_point_t *d, ec_pqr_point_t *s)
{
	//TODO: maybe we don't need to force adjust here.
	poly_copy(d->x, s->x, 1);
	poly_copy(d->y, s->y, 1);

	return d;
}

int ec_pqr_point_cmp(ec_pqr_point_t *p, ec_pqr_point_t *q)
{
	if (poly_cmp(p->x, q->x) == POLY_CMP_E && poly_cmp(p->y, q->y) == POLY_CMP_E)
		return POINT_CMP_E;
	return POINT_CMP_NE;
}

ec_pqr_point_t *ec_pqr_point_zero(ec_pqr_point_t *p)
{
	poly_zero(p->x);
	poly_zero(p->y);

	return p;
}

int ec_pqr_point_is_zero(ec_pqr_point_t *p)
{
	return (poly_is_zero(p->x) && poly_is_zero(p->y));
}

void ec_pqr_point_to_mon(ec_pqr_point_t *p)
{
	poly_to_mon(p->x);
	poly_to_mon(p->y);
}

void ec_pqr_point_from_mon(ec_pqr_point_t *p)
{
	poly_from_mon(p->x);
	poly_from_mon(p->y);
}

ec_pqr_point_t *ec_pqr_point_neg(ec_pqr_point_t *r, ec_pqr_point_t *p, ec_pqr_group_t *ecg)
{
	poly_copy(r->x, p->x, 1);
	poly_zero(r->y);
	pqr_sub(r->y, r->y, p->y, ecg->p);
	return r;
}

ec_pqr_point_t *ec_pqr_point_double(ec_pqr_point_t *r, ec_pqr_point_t *p, ec_pqr_group_t *ecg)
{
	// Handle trivial case.
	if (poly_is_zero(p->y))
	{
		ec_pqr_point_zero(r);
		return r;
	}

	poly_t *s = poly_alloc(ecg->p->degree - 1, ecg->p->N, 1), 
		*t = poly_alloc(ecg->p->degree - 1, ecg->p->N, 1);
	ec_pqr_point_t *pp = ec_pqr_point_copy(ec_pqr_point_alloc(ecg), p);
	poly_t *px = pp->x, *py = pp->y, *rx = r->x, *ry = r->y;

	_ec_pqr_square(t, px, ecg);       // t = px*px
	_ec_pqr_add(s, t, t, ecg);        // s = 2*px*px
	_ec_pqr_add(s, s, t, ecg);        // s = 3*px*px
	_ec_pqr_add(s, s, ecg->a, ecg);   // s = 3*px*px + a
	_ec_pqr_add(t, py, py, ecg);      // t = 2*py
	_ec_pqr_inv(t, t, ecg);           // t = 1/(2*py)
	_ec_pqr_mul(s, s, t, ecg);        // s = (3*px*px+a)/(2*py)

	_ec_pqr_square(rx, s, ecg);       // rx = s*s
	_ec_pqr_add(t, px, px, ecg);      // t = 2*px
	_ec_pqr_sub(rx, rx, t, ecg);      // rx = s*s - 2*px

	_ec_pqr_sub(t, px, rx, ecg);      // t = -(rx-px)
	_ec_pqr_mul(ry, s, t, ecg);       // ry = -s*(rx-px)
	_ec_pqr_sub(ry, ry, py, ecg);     // ry = -s*(rx-px) - py

	poly_free(s, 1);
	poly_free(t, 1);
	ec_pqr_point_free(pp);

	return r;
}

ec_pqr_point_t *ec_pqr_point_add(ec_pqr_point_t *r, ec_pqr_point_t *p, ec_pqr_point_t *q, ec_pqr_group_t *ecg)
{
	// Handle trivial cases.
	if (ec_pqr_point_is_zero(p))
	{
		ec_pqr_point_copy(r, q);
		return r;
	}

	if (ec_pqr_point_is_zero(q))
	{
		ec_pqr_point_copy(r, p);
		return r;
	}

	poly_t *s = poly_alloc(ecg->p->degree - 1, ecg->p->N, 1),
		*t = poly_alloc(ecg->p->degree - 1, ecg->p->N, 1),
		*u = poly_alloc(ecg->p->degree - 1, ecg->p->N, 1);
	ec_pqr_point_t *pp = ec_pqr_point_copy(ec_pqr_point_alloc(ecg), p),
		*qq = ec_pqr_point_copy(ec_pqr_point_alloc(ecg), q);
	poly_t *px = pp->x, *py = pp->y, *qx = qq->x, *qy = qq->y, *rx = r->x, *ry = r->y;

	// Handle limit cases.
	_ec_pqr_sub(u, qx, px, ecg);
	if (poly_is_zero(u))
	{
		_ec_pqr_sub(u, qy, py, ecg);
		if (poly_is_zero(u))
			ec_pqr_point_double(r, pp, ecg);
		else
			ec_pqr_point_zero(r);
		return r;
	}

	_ec_pqr_inv(t, u, ecg);        // t = 1/(qx-px)
	_ec_pqr_sub(u, qy, py, ecg);   // u = qy-py
	_ec_pqr_mul(s, t, u, ecg);     // s = (qy-py)/(qx-px)

	_ec_pqr_square(rx, s, ecg);    // rx = s*s
	_ec_pqr_add(t, px, qx, ecg);   // t = px+qx
	_ec_pqr_sub(rx, rx, t, ecg);   // rx = s*s - (px+qx)

	_ec_pqr_sub(t, px, rx, ecg);   // t = -(rx-px)
	_ec_pqr_mul(ry, s, t, ecg);    // ry = -s*(rx-px)
	_ec_pqr_sub(ry, ry, py, ecg);  // ry = -s*(rx-px) - py

	poly_free(s, 1);
	poly_free(t, 1);
	poly_free(u, 1);
	ec_pqr_point_free(pp);
	ec_pqr_point_free(qq);

	return r;
}

ec_pqr_point_t *ec_pqr_point_mul(ec_pqr_point_t *d, bn_t *a, ec_pqr_point_t *b, ec_pqr_group_t *ecg)
{
	int i;
	ec_pqr_point_t *bt = ec_pqr_point_copy(ec_pqr_point_alloc(ecg), b);

	ec_pqr_point_zero(d);

	for (i = 0; i <= bn_maxbit(a); i++)
	{
		if (bn_getbit(a, i))
			ec_pqr_point_add(d, d, bt, ecg);
		ec_pqr_point_double(bt, bt, ecg);
	}

	ec_pqr_point_free(bt);

	return d;
}
