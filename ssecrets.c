/*
* Copyright 2016 naehrwert
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#include "bn.h"
#include "poly.h"

static poly_t *_ssecrets_random_poly(u32 k, bn_t *N)
{
   u32 i;
   poly_t *res = poly_alloc(k - 1, N, 0);

   for(i = 1; i <= k - 1; i++)
   {
      bn_t *c = bn_reduce(bn_rand(bn_alloc(N->n)), N);
      poly_set_coeff(res, i, c);
   }

   return res;
}

poly_t *ssecrets_create_poly(bn_t *s, u32 k, bn_t *N)
{
   poly_t *res;

   // Create random polynomial with k-1 coefficients.
   res = _ssecrets_random_poly(k, N);
   // Set shared secret as coefficient 0.
   if(s != NULL)
      poly_set_coeff(res, 0, bn_copy(bn_alloc(N->n), s));
   else
      poly_set_coeff(res, 0, bn_alloc(N->n));

   return res;
}

bn_t *ssecrets_create_share(poly_t *p, bn_t *x)
{
   // Evaluate polynomial at x.
   return poly_eval(p, bn_alloc(p->N->n), x);
}

bn_t *ssecrets_calc_secret(bn_t **x, bn_t **s, u32 cnt, bn_t *N)
{
   u32 i, j;
   BOOL first;

   bn_t *res = bn_alloc(N->n), *t = bn_alloc(N->n), *a = bn_alloc(N->n), *b = bn_alloc(N->n);

   // Compute secret by using Lagrange polynomial interpolation algorithm for x = 0.
   // s = \sum s_i \prod_{j \ne i} (-x_j)(x_i - x_j)^{-1} \mod N
   for(i = 0; i < cnt; i++)
   {
      bn_zero(t);
      first = TRUE;
      for(j = 0; j < cnt; j++)
      {
         if(j != i)
         {
            // b = (x_i - x_j)^(-1)
            bn_mon_inv(b, bn_to_mon(bn_sub(b, x[i], x[j], N), N), N);
            // a = -x_j*b
            bn_mon_mul(a, bn_to_mon(bn_sub(a, bn_zero(a), x[j], N), N), b, N);
            // t *= a
            if(first == TRUE)
            {
               bn_to_mon(bn_copy(t, a), N);
               first = FALSE;
            }
            else
               bn_mon_mul(t, t, a, N);
         }
      }
      // res += s_i * t
      bn_add(res, res, bn_from_mon(bn_mon_mul(t, s[i], t, N), N), N);
   }

   bn_free(t);
   bn_free(a);
   bn_free(b);

   return res;
}
