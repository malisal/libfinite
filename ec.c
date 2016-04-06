/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <stdlib.h>
#include "ec.h"

static void _ec_add(bn_t *d, bn_t *a, bn_t *b, ec_group_t *ecg)
{
   bn_add(d, a, b, ecg->p);
}

static void _ec_sub(bn_t *d, bn_t *a, bn_t *b, ec_group_t *ecg)
{
   bn_sub(d, a, b, ecg->p);
}

static void _ec_mul(bn_t *d, bn_t *a, bn_t *b, ec_group_t *ecg)
{
   bn_mon_mul(d, a, b, ecg->p);
}

static void _ec_square(bn_t *d, bn_t *a, ec_group_t *ecg)
{
   _ec_mul(d, a, a, ecg);
}

static void _ec_inv(bn_t *d, bn_t *a, ec_group_t *ecg)
{
   bn_mon_inv(d, a, ecg->p);
}

ec_point_t *ec_point_alloc(uint32_t n)
{
   ec_point_t *res;

   if((res = (ec_point_t *)mem_alloc(sizeof(ec_point_t))) == NULL)
      return NULL;

   res->x = bn_alloc(n);
   res->y = bn_alloc(n);

   return res;
}

void ec_point_free(ec_point_t *p)
{
   bn_free(p->x);
   bn_free(p->y);
   mem_free(p);
}

ec_point_t *ec_point_copy(ec_point_t *d, ec_point_t *s)
{
   bn_copy(d->x, s->x);
   bn_copy(d->y, s->y);

   return d;
}

ec_point_t *ec_point_zero(ec_point_t *p)
{
   bn_zero(p->x);
   bn_zero(p->y);

   return p;
}

int ec_point_is_zero(ec_point_t *p)
{
   return (bn_is_zero(p->x) && bn_is_zero(p->y));
}

void ec_point_to_mon(ec_point_t *p, ec_group_t *ecg)
{
   bn_to_mon(p->x, ecg->p);
   bn_to_mon(p->y, ecg->p);
}

void ec_point_from_mon(ec_point_t *p, ec_group_t *ecg)
{
   bn_from_mon(p->x, ecg->p);
   bn_from_mon(p->y, ecg->p);
}

ec_point_t *ec_point_double(ec_point_t *r, ec_point_t *p, ec_group_t *ecg)
{
   // Handle trivial case.
   if(bn_is_zero(p->y))
   {
      ec_point_zero(r);
      return r;
   }

   bn_t *s = bn_alloc(ecg->p->n), *t = bn_alloc(ecg->p->n);
   ec_point_t *pp = ec_point_copy(ec_point_alloc(ecg->p->n), p);
   bn_t *px = pp->x, *py = pp->y, *rx = r->x, *ry = r->y;

   _ec_square(t, px, ecg);       // t = px*px
   _ec_add(s, t, t, ecg);        // s = 2*px*px
   _ec_add(s, s, t, ecg);        // s = 3*px*px
   _ec_add(s, s, ecg->a, ecg);   // s = 3*px*px + a
   _ec_add(t, py, py, ecg);      // t = 2*py
   _ec_inv(t, t, ecg);           // t = 1/(2*py)
   _ec_mul(s, s, t, ecg);        // s = (3*px*px+a)/(2*py)

   _ec_square(rx, s, ecg);       // rx = s*s
   _ec_add(t, px, px, ecg);      // t = 2*px
   _ec_sub(rx, rx, t, ecg);      // rx = s*s - 2*px

   _ec_sub(t, px, rx, ecg);      // t = -(rx-px)
   _ec_mul(ry, s, t, ecg);       // ry = -s*(rx-px)
   _ec_sub(ry, ry, py, ecg);     // ry = -s*(rx-px) - py

   bn_free(s);
   bn_free(t);
   ec_point_free(pp);

   return r;
}

ec_point_t *ec_point_add(ec_point_t *r, ec_point_t *p, ec_point_t *q, ec_group_t *ecg)
{
   // Handle trivial cases.
   if(ec_point_is_zero(p))
   {
      ec_point_copy(r, q);
      return r;
   }

   if(ec_point_is_zero(q))
   {
      ec_point_copy(r, p);
      return r;
   }

   bn_t *s = bn_alloc(ecg->p->n), *t = bn_alloc(ecg->p->n), *u = bn_alloc(ecg->p->n);
   ec_point_t *pp = ec_point_copy(ec_point_alloc(ecg->p->n), p), *qq = ec_point_copy(ec_point_alloc(ecg->p->n), q);
   bn_t *px = pp->x, *py = pp->y, *qx = qq->x, *qy = qq->y, *rx = r->x, *ry = r->y;

   // Handle limit cases.
   _ec_sub(u, qx, px, ecg);
   if(bn_is_zero(u))
   {
      _ec_sub(u, qy, py, ecg);
      if(bn_is_zero(u))
         ec_point_double(r, pp, ecg);
      else
         ec_point_zero(r);
      return r;
   }

   _ec_inv(t, u, ecg);        // t = 1/(qx-px)
   _ec_sub(u, qy, py, ecg);   // u = qy-py
   _ec_mul(s, t, u, ecg);     // s = (qy-py)/(qx-px)

   _ec_square(rx, s, ecg);    // rx = s*s
   _ec_add(t, px, qx, ecg);   // t = px+qx
   _ec_sub(rx, rx, t, ecg);   // rx = s*s - (px+qx)

   _ec_sub(t, px, rx, ecg);   // t = -(rx-px)
   _ec_mul(ry, s, t, ecg);    // ry = -s*(rx-px)
   _ec_sub(ry, ry, py, ecg);  // ry = -s*(rx-px) - py

   bn_free(s);
   bn_free(t);
   bn_free(u);
   ec_point_free(pp);
   ec_point_free(qq);

   return r;
}

ec_point_t *ec_point_mul(ec_point_t *d, bn_t *a, ec_point_t *b, ec_group_t *ecg)
{
   bn_t *at = bn_copy(bn_alloc(a->n), a);
   ec_point_t *bt = ec_point_copy(ec_point_alloc(b->x->n), b);

   ec_point_zero(d);

   while(bn_cmp_ui(at, 0))
   {
      if(bn_lsb(at))
         ec_point_add(d, d, bt, ecg);

      ec_point_double(bt, bt, ecg);

      bn_rshift(at, 1);
   }

   bn_free(at);
   ec_point_free(bt);

   return d;
}

