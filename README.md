## libfinite

```
   Copyright 2016 naehrwert
   Copyright 2016 Luka Malisa <luka.malisha@gmail.com>
   Licensed under the terms of the GNU GPL, version 2
   http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
```

### Introduction

`libfinite` is a small and fast bignum library. 

```bash

```

##### `bn.c`
Bignum library for arithmetic in Z/nZ (the positive half, since we don't care about signs) which can also be used to describe operations over some finite field GF(p) (F_p). Note that in some of the cases below, the use of a finite field can be substituted for a more general Z/nZ.

##### `dh.c`
This implements Diffie-Hellman key exchange.

##### `ec.c`
Given an elliptic curve in Weierstrass form (y^2 = x^3 + ax + b) over a finite field, this gives a representation of the group of points on the elliptic curve, i.e. point addition, point doubling, multiplication by a number.

##### `ecdsa.c`
This implements the Elliptic Curve Digital Signature Algorithm (ECDSA).

##### `ecnr.c`
This implements Elliptic Curve Nyberg Rueppel (ECNR), the Nyberg-Rueppel signature scheme over the elliptic curve group.

##### `inr.c`
This implements the Nyberg-Rueppel signature scheme over a finite field.

##### `pc.c`
Given a Pell conic (x^2 - Dy^2 = 1) over a finite field, this gives a representation of the group of points on the Pell conic, i.e. point addition, multiplication by a number.
Assuming the finite field is F_p, the parameter D should be a non-square in the field since then the group is known to be isomorphic to the multiplicative subgroup of F_p^2. Otherwise it will simply be isomorphic to the multiplicative subgroup of F_p.

##### `pcnr.c`
This implements Pell Conic Nyberg Rueppel (PCNR), the Nyberg-Rueppel signature scheme over the Pell conic group.

##### `poly.c`
This implements operations for the ring of polynomials R over some finite field F_p (R = F_p[X]), i.e. addition, subtraction, multiplication, division and remainder.

##### `pqr.c`
Given a polynomial ring R and a polynomial I in R, this implements operations for the polynomial quotient ring R/(I), i.e. addition, subtraction, multiplication, inversion, exponentiation. This can be used to construct a representation of certain finite fields, e.g. F_q where q = p^n, p prime, n some positive integer, i.e. take I to be an irreducible polynomial of degree n over the finite field F_p.

##### `ssecrets.c`
This implements Shamir's Secret Sharing.

### Examples

TODO
