/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include <assert.h>
#include <stdlib.h>
#include <stdarg.h>

#include "poly.h"

/* Wrapper for safe calling fast poly op. */
#define _FAST_WRAPPER(op, d, a, b) \
   do { \
      if (d == a && d != b) \
      { \
         poly_t *t = poly_alloc(d->degree, d->N, 1); \
         poly_copy(t, a, 0); \
         poly_##op##_fast(d, t, b); \
         poly_free(t, 1); \
         return d; \
      } \
      else if (d == b && d != a) \
      { \
         poly_t *t = poly_alloc(d->degree, d->N, 1); \
         poly_copy(t, b, 0); \
         poly_##op##_fast(d, a, t); \
         poly_free(t, 1); \
         return d; \
      } \
      else if (d == a && d == b) \
      { \
         poly_t *t = poly_alloc(d->degree, d->N, 1); \
         poly_copy(t, a, 0); \
         poly_##op##_fast(d, t, t); \
         poly_free(t, 1); \
         return d; \
      } \
      return poly_##op##_fast(d, a, b); \
   } while (0)

poly_t *poly_read(FILE *fp, poly_t *dst)
{
	int i;

	for (i = 0; i <= dst->degree; i++)
		bn_read(fp, dst->coeffs[i]);

	return dst;
}

poly_t *poly_write(FILE *fp, poly_t *p)
{
	int i;

	for (i = 0; i <= p->degree; i++)
		bn_write(fp, p->coeffs[i]);

	return p;
}

void poly_print(FILE *fp, const s8 *pre, poly_t *p, const s8 *post)
{
	int i, not_first = 0;

	if (p == NULL)
	{
		fprintf(fp, "%s(NULL)%s", pre, post);
		return;
	}

	fputs((char *)pre, fp);
	for (i = 0; i <= p->degree; i++)
	{
		if (bn_is_zero(p->coeffs[i]))
			continue;

		if (not_first)
			fputs(" + ", fp);
		else
			not_first = 1;

		bn_print(fp, "", bn_from_mon(p->coeffs[i], p->N), "");
		bn_to_mon(p->coeffs[i], p->N);

		if (i == 1)
			fputs("*X", fp);
		else if (i > 1)
			fprintf(fp, "*X^%d", i);
	}
	fputs((char *)post, fp);
}

poly_t *poly_alloc(int degree, bn_t *N, int alloc_coeffs)
{
	poly_t *res;

	if ((res = (poly_t *)mem_alloc(sizeof(poly_t))) == NULL)
		return NULL;

	// Degree + 1 coeffs.
	if ((res->coeffs = (bn_t **)mem_alloc(sizeof(bn_t *) * (degree + 1))) == NULL)
	{
		mem_free(res);
		return NULL;
	}

	res->degree = degree;
	res->N = N;

	if (alloc_coeffs)
	{
		int i;
		for (i = 0; i < degree + 1; i++)
			res->coeffs[i] = bn_alloc(N->n);
	}

	return res;
}

void poly_free(poly_t *p, int free_coeffs)
{
	if (free_coeffs)
	{
		int i;
		for (i = 0; i < p->degree + 1; i++)
			bn_free(p->coeffs[i]);
	}

	mem_free(p->coeffs);
	mem_free(p);
}

poly_t *poly_adjust(poly_t *p, int degree, int alloc_free_coeffs)
{
	if (p->degree == degree)
		return p;

	if (p->degree < degree)
	{
		// Grow polynomial.
		p->coeffs = (bn_t **)realloc(p->coeffs, sizeof(bn_t *) * (degree + 1));

		if (alloc_free_coeffs)
		{
			// Allocate new coefficients.
			int i;
			for (i = p->degree + 1; i < degree + 1; i++)
				p->coeffs[i] = bn_zero(bn_alloc(p->N->n));
		}

		p->degree = degree;
	}
	else if (p->degree > degree)
	{
		if (alloc_free_coeffs)
		{
			// Free old coefficients
			int i;
			for (i = degree + 1; i < p->degree + 1; i++)
				bn_free(p->coeffs[i]);
		}

		// Shrink polynomial.
		p->coeffs = (bn_t **)realloc(p->coeffs, sizeof(bn_t *) * (degree + 1));
		p->degree = degree;
	}

	return p;
}

poly_t *poly_copy(poly_t *d, poly_t *s, int adjust)
{
	int i;

	if (adjust)
	{
		// Adjust degree and copy all of the coefficients.
		if (d->degree != s->degree)
			poly_adjust(d, s->degree, 1);
		for (i = 0; i < d->degree + 1; i++)
			bn_copy(d->coeffs[i], s->coeffs[i]);

		return d;
	}

	int to = d->degree < s->degree ? d->degree + 1 : s->degree + 1;

	// Copy or truncate.
	for (i = 0; i < to; i++)
		bn_copy(d->coeffs[i], s->coeffs[i]);

	// Zero out higher degrees.
	for (; i < d->degree + 1; i++)
		bn_zero(d->coeffs[i]);

	return d;
}

poly_t *poly_concat(poly_t *d, poly_t *a, poly_t *b)
{
	int i;

	poly_adjust(d, a->degree + b->degree + 1, 1);
	for (i = 0; i < a->degree + 1; i++)
		bn_copy(d->coeffs[i], a->coeffs[i]);
	for (; i < a->degree + 1 + b->degree + 1; i++)
		bn_copy(d->coeffs[i], b->coeffs[i - a->degree - 1]);

	return d;
}

poly_t *poly_zero(poly_t *p)
{
	int i;

	for (i = 0; i < p->degree + 1; i++)
		bn_zero(p->coeffs[i]);

	return p;
}

int poly_is_zero(poly_t *p)
{
	int i;

	assert(p->degree >= 0);

	// Check if all the coeffs are 0.
	for (i = 0; i < p->degree + 1; i++)
		if (!bn_is_zero(p->coeffs[i]))
			return 0;

	return 1;
}

poly_t *poly_one(poly_t *p)
{
	assert(p->degree >= 0 && p->coeffs[0] != NULL);

	poly_zero(p);
	bn_set_ui(p->coeffs[0], 1);
	bn_to_mon(p->coeffs[0], p->N);

	return p;
}

int poly_is_one(poly_t *p)
{
	int i, is_one = 0;

	assert(p->degree >= 0 && p->coeffs[0] != NULL);

	//Check if the linear coeff is 1.
	bn_from_mon(p->coeffs[0], p->N);
	is_one = bn_cmp_ui(p->coeffs[0], 1) == BN_CMP_E;
	bn_to_mon(p->coeffs[0], p->N);
	if (!is_one)
		return 0;

	// Return 1 if we only have a linear term.
	if (p->degree == 0)
		return 1;

	// Check if the other coeffs are 0.
	for (i = 1; i < p->degree + 1; i++)
		if (!bn_is_zero(p->coeffs[i]))
			return 0;

	return 1;
}

int poly_deg(poly_t *p)
{
	int i;

	for (i = p->degree; i >= 0 && bn_is_zero(p->coeffs[i]); i--);

	return i;
}

int poly_set_coeff(poly_t *p, int i, bn_t *coeff)
{
	if (i > p->degree)
		return 0;

	// Convert coefficients for faster eval.
	p->coeffs[i] = bn_to_mon(coeff, p->N);
	return 1;
}

int poly_free_coeff(poly_t *p, int i)
{
	if (i > p->degree)
		return 0;

	bn_free(p->coeffs[i]);
	p->coeffs[i] = NULL;
	return 1;
}

poly_t *poly_from_fmt(poly_t *p, const char *fmt, ...)
{
	va_list ap;
	int i;

	va_start(ap, fmt);

	for (i = 0; fmt[i] && i <= p->degree; i++)
	{
		switch (fmt[i])
		{
		case 'i':
			bn_set_ui(p->coeffs[i], va_arg(ap, int));
			break;
		case 'I':
			bn_to_mon(bn_set_ui(p->coeffs[i], va_arg(ap, int)), p->N);
			break;
		case 'q':
			bn_set_ui(p->coeffs[i], va_arg(ap, ull_t));
			break;
		case 'Q':
			bn_to_mon(bn_set_ui(p->coeffs[i], va_arg(ap, ull_t)), p->N);
			break;
		case 's':
			bn_from_str(p->coeffs[i], va_arg(ap, char *));
			break;
		case 'S':
			bn_to_mon(bn_from_str(p->coeffs[i], va_arg(ap, char *)), p->N);
			break;
		case '0':
			bn_zero(p->coeffs[i]);
			break;
		}
	}

	va_end(ap);

	return p;
}

poly_t *poly_to_mon(poly_t *p)
{
	int i;

	for (i = 0; i <= p->degree; i++)
		bn_to_mon(p->coeffs[i], p->N);

	return p;
}

poly_t *poly_from_mon(poly_t *p)
{
	int i;

	for (i = 0; i <= p->degree; i++)
		bn_from_mon(p->coeffs[i], p->N);

	return p;
}

bn_t *poly_eval(poly_t *p, bn_t *dst, bn_t *x)
{
	int i;
	bn_t *e = bn_alloc(4);
	bn_t *t = bn_alloc(p->N->n);
	bn_t *tx = bn_to_mon(bn_copy(bn_alloc(x->n), x), p->N);

	// TODO: Use Horner's method for faster eval.
	for (i = 0; i <= p->degree; i++)
	{
		bn_set_ui(e, i);
		// t = x^e
		bn_mon_pow(t, tx, p->N, e);
		// dst += t * a_i
		bn_add(dst, dst, bn_mon_mul(t, t, p->coeffs[i], p->N), p->N);
	}

	bn_free(e);
	bn_free(t);
	bn_free(tx);

	return dst;
}

poly_t *poly_add_fast(poly_t *d, poly_t *a, poly_t *b)
{
	int i, md;

	assert(d != a && d != b);

	// Grow d if the degree is too small.
	md = MAX(a->degree, b->degree);
	if (d->degree < md)
		poly_adjust(d, md, 1);

	if (a->degree < b->degree)
	{
		for (i = 0; i < a->degree + 1; i++)
			bn_add(d->coeffs[i], a->coeffs[i], b->coeffs[i], d->N);
		for (; i < b->degree + 1; i++)
			bn_copy(d->coeffs[i], b->coeffs[i]);
	}
	else
	{
		for (i = 0; i < b->degree + 1; i++)
			bn_add(d->coeffs[i], b->coeffs[i], a->coeffs[i], d->N);
		for (; i < a->degree + 1; i++)
			bn_copy(d->coeffs[i], a->coeffs[i]);
	}

	return d;
}

poly_t *poly_add(poly_t *d, poly_t *a, poly_t *b)
{
	_FAST_WRAPPER(add, d, a, b);
}

poly_t *poly_sub_fast(poly_t *d, poly_t *a, poly_t *b)
{
	int i, md;

	assert(d != a && d != b);

	//Grow d if the degree is too small.
	md = MAX(a->degree, b->degree);
	if (d->degree < md)
		poly_adjust(d, md, 1);

	if (a->degree < b->degree)
	{
		for (i = 0; i < a->degree + 1; i++)
			bn_sub(d->coeffs[i], a->coeffs[i], b->coeffs[i], d->N);
		for (; i < b->degree + 1; i++)
		{
			bn_zero(d->coeffs[i]);
			bn_sub(d->coeffs[i], d->coeffs[i], b->coeffs[i], d->N);
		}
	}
	else
	{
		for (i = 0; i < b->degree + 1; i++)
			bn_sub(d->coeffs[i], a->coeffs[i], b->coeffs[i], d->N);
		for (; i < a->degree + 1; i++)
			bn_copy(d->coeffs[i], a->coeffs[i]);
	}

	return d;
}

poly_t *poly_sub(poly_t *d, poly_t *a, poly_t *b)
{
	_FAST_WRAPPER(sub, d, a, b);
}

poly_t *poly_mul_fast(poly_t *d, poly_t *a, poly_t *b)
{
	int i, j;

	assert(d != a && d != b);

	// TODO: could just always keep this in the poly_t struct for speed.
	bn_t *t = bn_alloc(d->N->n);

	// Adjust degree of d and zero out coefficients.
	poly_adjust(d, a->degree + b->degree, 1);
	poly_zero(d);

	for (i = 0; i < a->degree + 1; i++)
		for (j = 0; j < b->degree + 1; j++)
		{
			bn_mon_mul(t, a->coeffs[i], b->coeffs[j], d->N);
			bn_add(d->coeffs[i + j], d->coeffs[i + j], t, d->N);
		}

	bn_free(t);

	return d;
}

poly_t *poly_mul(poly_t *d, poly_t *a, poly_t *b)
{
	_FAST_WRAPPER(mul, d, a, b);
}

poly_t *poly_mulc(poly_t *d, poly_t *a, bn_t *b)
{
	int i;

	if (d->degree != a->degree)
		poly_adjust(d, a->degree, 1);

	for (i = 0; i < d->degree + 1; i++)
		bn_mon_mul(d->coeffs[i], a->coeffs[i], b, d->N);

	return d;
}

poly_t *poly_div(poly_t *q, poly_t *r, poly_t *a, poly_t *b)
{
	int i, dega = poly_deg(a), degb = poly_deg(b);

	poly_t *t = poly_alloc(-1, q->N, 0);
	poly_t *tt = poly_alloc(-1, q->N, 0);
	poly_t *ta = poly_alloc(dega, q->N, 1);
	poly_t *tq = poly_alloc(-1, q->N, 0);

	bn_t *bc = bn_alloc(q->N->n);
	bn_t *c = bn_alloc(q->N->n);
	bn_t *ctb = bn_alloc(q->N->n);

	// Truncate a.
	poly_copy(ta, a, 0);

	// Calculate inverse of lead coefficient.
	bn_mon_inv(bc, b->coeffs[degb], q->N);

	while (ta->degree + 1 >= degb + 1)
	{
		poly_adjust(t, degb, 1);
		for (i = 0; i < degb + 1; i++)
			bn_copy(t->coeffs[i], ta->coeffs[ta->degree - degb + i]);

		bn_mon_mul(c, t->coeffs[t->degree], bc, q->N);
		poly_adjust(tq, tq->degree + 1, 1);
		for (i = tq->degree; i > 0; i--)
			bn_copy(tq->coeffs[i], tq->coeffs[i - 1]);
		bn_copy(tq->coeffs[0], c);

		for (i = 0; i < degb + 1; i++)
		{
			bn_mon_mul(ctb, c, b->coeffs[i], q->N);
			bn_sub(t->coeffs[i], t->coeffs[i], ctb, q->N);
		}

		poly_adjust(tt, ta->degree - degb - 1, 1);
		poly_copy(tt, ta, 0); //Truncate ta.
		poly_adjust(t, t->degree - 1, 1);
		poly_concat(ta, tt, t);
	}

	poly_copy(q, tq, 1);
	poly_copy(r, ta, 1);

	poly_free(tq, 1);
	poly_free(ta, 1);
	poly_free(tt, 1);
	poly_free(t, 1);

	bn_free(ctb);
	bn_free(c);
	bn_free(bc);

	return q;
}

poly_t *poly_rem(poly_t *r, poly_t *a, poly_t *b)
{
	int i, dega = poly_deg(a), degb = poly_deg(b);

	poly_t *t = poly_alloc(-1, r->N, 0);
	poly_t *tt = poly_alloc(-1, r->N, 0);
	poly_t *ta = poly_alloc(dega, r->N, 1);

	bn_t *bc = bn_alloc(r->N->n);
	bn_t *c = bn_alloc(r->N->n);
	bn_t *ctb = bn_alloc(r->N->n);

	// Truncate a.
	poly_copy(ta, a, 0);

	// Calculate inverse of lead coefficient.
	bn_mon_inv(bc, b->coeffs[degb], r->N);

	while (ta->degree + 1 >= degb + 1)
	{
		poly_adjust(t, degb, 1);
		for (i = 0; i < degb + 1; i++)
			bn_copy(t->coeffs[i], ta->coeffs[ta->degree - degb + i]);

		bn_mon_mul(c, t->coeffs[t->degree], bc, r->N);

		for (i = 0; i < degb + 1; i++)
		{
			bn_mon_mul(ctb, c, b->coeffs[i], r->N);
			bn_sub(t->coeffs[i], t->coeffs[i], ctb, r->N);
		}

		poly_adjust(tt, ta->degree - degb - 1, 1);
		poly_copy(tt, ta, 0); //Truncate ta.
		poly_adjust(t, t->degree - 1, 1);
		poly_concat(ta, tt, t);
	}

	poly_copy(r, ta, 1);

	poly_free(ta, 1);
	poly_free(tt, 1);
	poly_free(t, 1);

	bn_free(ctb);
	bn_free(c);
	bn_free(bc);

	return r;
}
