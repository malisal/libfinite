#include <stdio.h>

#include "bn.h"
#include "poly.h"
#include "pqr.h"

int main(int argc, char **argv)
{
	bn_t *N = bn_from_str(bn_alloc(9), "01000000000000000D");

	poly_t *I = poly_from_fmt(poly_alloc(3, N, 1), "IIII", 1, 1, 0, 1); //1 + x + x^3
	poly_print(stdout, "I = ", I, "\n");

	poly_t *q = poly_from_fmt(poly_alloc(2, N, 1), "III", 1, 2, 1); //1 + 2x + x^2
	poly_print(stdout, "q = ", q, "\n");

	poly_t *qi = poly_alloc(2, N, 1);
	pqr_inv_fast(qi, q, I);
	poly_print(stdout, "qi = ", qi, "\n");

	poly_t *p = poly_alloc(2, N, 1);
	pqr_mul_fast(p, q, qi, I);
	poly_print(stdout, "q*qi = ", p, "\n");

	return 0;
}
