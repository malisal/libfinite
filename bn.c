/*
* Copyright 2016 Luka Malisa <luka.malisha@gmail.com>
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#if !defined(NAKED)
   #include <stdio.h>
   #include <stdlib.h>
   #include <string.h>
   #include <assert.h>
   #include <sys/types.h>
   #include <sys/stat.h>
   #include <fcntl.h>

   #if defined(_WIN32) || defined(_MSC_VER)
      #include <windows.h>
   #endif
#endif

#include "bn.h"

// Fast Montgomery initialization (taken from PolarSSL)
static void _bn_mon_init(bn_t *n)
{
   ul_t x, m0 = n->l[0];

   x  = m0;
   x += ((m0 + 2) & 4) << 1;

   for(int i = BN_LIMB_BITS; i >= 8; i /= 2)
      x *= (2 - (m0 * x));

   n->mp = ~x + 1;
}

static u8 _bn_str_to_u8(const s8 *str)
{
   u32 i = 2;
   u8 t, res = 0;
   u8 c;

   while(i--)
   {
      c = *str++;
      if(c >= '0' && c <= '9')
         t = c - '0';
      else if(c >= 'a' && c <= 'f')
         t = c - 'a' + 10;
      else if(c >= 'A' && c <= 'F')
         t = c - 'A' + 10;
      else
         t = 0;
      res |= t << (i * 4);
   }

   return res;
}

static int _bn_add(bn_t *d, bn_t *a, bn_t *b)
{
   ull_t C = 0;

//   assert(a->n == b->n);

   for(int i = 0; i < MIN(a->n_limbs, b->n_limbs); i++)
   {
      C += (ull_t)a->l[i] + b->l[i];
      d->l[i] = C;

      C >>= BN_LIMB_BITS;
   }

   return C;
}

static int _bn_add_ui(bn_t *d, bn_t *a, ul_t b)
{
   ull_t C = b;

   for(int i = 0; i < a->n_limbs; i++)
   {
      C += (ull_t)a->l[i];
      d->l[i] = C;

      C >>= BN_LIMB_BITS;
   }

   return C;
}

static int _bn_sub(bn_t *d, bn_t *a, bn_t *b)
{
   ull_t C = 1;

//   assert(a->n == b->n);

   for(int i = 0; i < MIN(a->n_limbs, b->n_limbs); i++)
   {
      C += (ull_t)a->l[i] + BN_MAX_DIGIT - b->l[i];
      d->l[i] = C;

      C >>= BN_LIMB_BITS;
   }

   return 1 - C;
}

static int _bn_sub_ui(bn_t *d, bn_t *a, ul_t b)
{
   ull_t C = 1 + BN_MAX_DIGIT - b;

   for(int i = 0; i < a->n_limbs; i++)
   {
      C += (ull_t)a->l[i];
      d->l[i] = C;

      C >>= BN_LIMB_BITS;
   }

   return 1 - C;
}

// Runtime endianess detection. A bit slower, but no macro magic needed
static int _bn_big_endian()
{
   u32 t = 0x11223344;
   u8 *p = (u8 *)&t;

   if(p[0] == 0x11)
      return 1;
   else
      return 0;
}

static bn_t *_bn_lshift_limbs(bn_t *a, int n)
{
   int i;

   for(i = a->n_limbs; i >= n; i--)
      a->l[i] = a->l[i-n];

   while(i >= 0)
      a->l[i--] = 0;

   return a;
}

static bn_t *_bn_lshift(bn_t *a, int b)
{
   ul_t mask = -1;
   ul_t prev_c = 0;
   ul_t c = 0;

   if(b != BN_LIMB_BITS)
      mask = ~(mask >> b);

   for(int x = 0; x < a->n_limbs; x++)
   {
      c = a->l[x] & mask;

      // Shift left by amount
      a->l[x] <<= b;

      // Add the carry part from the previous limb
      a->l[x] |= prev_c >> (BN_LIMB_BITS - b);

      prev_c = c;
   }

   return a;
}

static bn_t *_bn_rshift_limbs(bn_t *a, int n)
{
   // Makes the code run faster
   memmove(a->l, &a->l[n], (a->n_limbs - n) * BN_LIMB_BYTES);

/*
   int x;

   for(x = 0; x < a->n_limbs - n; x++)
      a->l[x] = a->l[x+n];

   while(x <= n)
      a->l[x++] = 0;
*/
   return a;
}

static bn_t *_bn_rshift(bn_t *a, int b)
{
   ull_t mask;
   ul_t prev_c = 0;
   ul_t c = 0;

   // Create the mask
   mask = ((ull_t)1 << b) - 1;

   for(int x = a->n_limbs - 1; x >= 0; x--)
   {
      c = a->l[x] & mask;

      // Shift right by amount
      a->l[x] >>= b;

      // Add the carry part from the previous limb
      a->l[x] |= prev_c << (BN_LIMB_BITS - b);

      prev_c = c;
   }

   return a;
}

// Helper function which multiplies and adds in a single run.
static bn_t *_bn_mad_ui(bn_t *d, bn_t *a, ul_t b)
{
   // D = D + A * b
   int x;
   ull_t S = 0;

   // Can a hold the result?
   assert(d->n >= (a->n + 2));

   for(x = 0; x < a->n_limbs; x++)
   {
      S += (ull_t)a->l[x] * (ull_t)b + (ull_t)d->l[x];
      d->l[x] = S;

      S >>= BN_LIMB_BITS;
   }

   // Add in the remaining carry
   while(S)
   {
      S += d->l[x];
      d->l[x] = S;

      S >>= BN_LIMB_BITS;
      x++;
   }

   return d;
}

int bn_maxbit(bn_t *a)
{
   for(int x = a->n_limbs * BN_LIMB_BITS; x >= 0; x--)
      if(bn_getbit(a, x) != 0)
         return x;

   return 0;
}

int bn_getbit(bn_t *a, int x)
{
   return (a->l[x / BN_LIMB_BITS] >> (x % BN_LIMB_BITS)) & 1;
}

void bn_setbit(bn_t *a, int x)
{
   a->l[x / BN_LIMB_BITS] |= 1 << (x % BN_LIMB_BITS);
}

bn_t *bn_from_bin(bn_t *a, s8 *s, int len)
{
   int x, y, z, w;

   ul_t limb;
   u8 *p_limb = (u8 *)&limb;

   for(x = len - 1, y = 0; x >= 0; y++)
   {
      limb = 0;

      for(z = 0, w = 0; z < BN_LIMB_BYTES; z += 1, w += 1)
      {
         if(x < 0)
            break;

         p_limb[w] = s[x];

         x -= 1;
      }

      if(_bn_big_endian())
         limb = SWAP(limb);

      a->l[y] = limb;
   }

   return a;
}

u8 *bn_to_bin(u8 *s, bn_t *a)
{
   ul_t *p = (ul_t *)s;

   int be = _bn_big_endian();

   for(int x = a->n_limbs - 1, y = 0; x >= 0; x--, y++)
   {
      if(be)
         p[y] = a->l[x];
      else
         p[y] = SWAP(a->l[x]);
   }

   return s;
}

bn_t *bn_from_str(bn_t *a, const s8 *s)
{
   int len = strlen(s);
   int x, y, z, w;

   ul_t limb;
   u8 *p_limb = (u8 *)&limb;

   // The x2 is for two ascii chars representing a single byte
   int step = (BN_LIMB_BYTES * 2);

   if(len % 2)
      return NULL;

   for(x = len - 2, y = 0; x >= 0; y++)
   {
      limb = 0;

      for(z = 0, w = 0; z < step && x >= 0; z += 2, w += 1)
      {
         p_limb[w] = _bn_str_to_u8(&s[x]);

         x -= 2;
      }

      if(_bn_big_endian())
         limb = SWAP(limb);

      a->l[y] = limb;
   }

   return a;
}

bn_t *bn_zero(bn_t *a)
{
   memset((char *)a->l, 0, a->n_limbs * BN_LIMB_BYTES);
   return a;
}

bn_t *bn_alloc(int size)
{
   int s;

   bn_t *ret = (bn_t *)mem_alloc(sizeof(bn_t));
   memset((char *)ret, 0x00, sizeof(bn_t));

   ret->n = size;
   ret->n_limbs = BYTES_TO_LIMBS(size);

   // Always allocate 4 limbs more than we need, so that potential bn_mon_mul is faster
   s = sizeof(ul_t) * (ret->n_limbs + 4);
   ret->l = (ul_t *)mem_alloc(s);
   memset((char *)ret->l, 0x00, s);

   return ret;
}

bn_t *bn_alloc_limbs(int limbs)
{
   return bn_alloc(LIMBS_TO_BYTES(limbs));
}

bn_t *bn_copy(bn_t *a, bn_t *b)
{
   int s = MIN(a->n_limbs, b->n_limbs);

   bn_zero(a);
   memcpy((s8 *)a->l, (s8 *)b->l, sizeof(ul_t) * s);

   return a;
}

void bn_free(bn_t *a)
{
   // Zero it out, just for good measure
   bn_zero(a);

   mem_free(a->l);
   mem_free(a);
}

inline bn_t *bn_set_ui(bn_t *a, u64 val)
{
   for(int x = 0; x < sizeof(val) / sizeof(ul_t); x++)
   {
      a->l[x] = val & (ul_t)-1;
      val >>= BN_LIMB_BITS;
   }

   return a;
}

bn_t *bn_add(bn_t *d, bn_t *a, bn_t *b, bn_t *n)
{
   // D = A + B % N

   // Prevent overflow
   if(_bn_add(d, a, b))
      _bn_sub(d, d, n);

   bn_reduce(d, n);

   return d;
}

bn_t *bn_add_ui(bn_t *d, bn_t *a, unsigned int b, bn_t *n)
{
   // D = A + B % N

   // Prevent overflow
   if(_bn_add_ui(d, a, b))
      _bn_sub(d, d, n);

   bn_reduce(d, n);

   return d;
}

bn_t *bn_sub(bn_t *d, bn_t *a, bn_t *b, bn_t *n)
{
   // D = A - B % N

   // Prevent underflow
   if(_bn_sub(d, a, b))
      _bn_add(d, d, n);

   bn_reduce(d, n);

   return d;
}

bn_t *bn_sub_ui(bn_t *d, bn_t *a, unsigned int b, bn_t *n)
{
   // D = A - B % N

   // Prevent underflow
   if(_bn_sub_ui(d, a, b))
      _bn_add(d, d, n);

   bn_reduce(d, n);

   return d;
}

int bn_cmp(bn_t *a, bn_t *b)
{
   //   assert(a->n_limbs == b->n_limbs);

   // First check the high limbs, if any
   if(a->n_limbs > b->n_limbs)
   {
      for(int i = b->n_limbs; i < a->n_limbs; i++)
         if(a->l[i])
            return BN_CMP_G;
   }
   else if(b->n_limbs > a->n_limbs)
   {
      for(int i = a->n_limbs; i < b->n_limbs; i++)
         if(b->l[i])
            return BN_CMP_G;
   }

   // ...then check the main limbs
   for(int x = MIN(a->n_limbs, b->n_limbs); x >= 0; x--)
   {
      if(a->l[x] < b->l[x])
         return BN_CMP_L;
      else if(a->l[x] > b->l[x])
         return BN_CMP_G;
   }

   return BN_CMP_E;
}

int bn_cmp_ui(bn_t *a, ul_t b)
{
   int ret = 0;

   if(a->l[0] < b)
      ret = BN_CMP_L;
   else if(a->l[0] > b)
      ret = BN_CMP_G;
   else if(a->l[0] == b)
      ret = BN_CMP_E;

   // Let's walk over all other digits of a (if any)
   for(int x = 1; x < a->n_limbs; x++)
   {
      if(a->l[x] != 0)
      {
         // So a has some non-zero digits. It's clearly bigger than b
         ret = BN_CMP_G;
         break;
      }
   }

   return ret;
}

int bn_is_zero(bn_t *a)
{
   for(int x = 0; x < a->n_limbs; x++)
      if(a->l[x] != 0)
         return 0;

   return 1;
}

bn_t *bn_reduce(bn_t *a, bn_t *n)
{
   while(bn_cmp(a, n) >= 0)
      _bn_sub(a, a, n);

   return a;
}

bn_t *bn_reduce_slow(bn_t *a, bn_t *n)
{
   bn_t *q = bn_alloc_limbs(n->n_limbs);
   bn_t *r = bn_alloc_limbs(n->n_limbs);

   bn_divrem(q, r, a, n);
   bn_copy(a, r);

   bn_free(q);
   bn_free(r);

   return a;
}

bn_t *bn_lshift(bn_t *a, int b)
{
   // Single largest shift we can do is one limb
   while(b > BN_LIMB_BITS)
   {
      _bn_lshift_limbs(a, 1);
      b -= BN_LIMB_BITS;
   }

   return _bn_lshift(a, b);
}

bn_t *bn_rshift(bn_t *a, int b)
{
   // Single largest shift we can do is one limb
   while(b > BN_LIMB_BITS)
   {
      _bn_rshift_limbs(a, 1);
      b -= BN_LIMB_BITS;
   }

   return _bn_rshift(a, b);
}

int bn_lsb(bn_t *a)
{
   return a->l[0] & 1;
}

int bn_msb(bn_t *a)
{
   return a->l[a->n - 1] >> (BN_LIMB_BITS - 1);
}

#if defined(BN_PRINT_FUNCS)
void bn_print(FILE *fp, const s8 *pre, bn_t *a, const s8 *post)
{
   int i, init = 1;

   fputs((char *)pre, fp);

   //Skip zero limbs.
   for(i = a->n_limbs - 1; i >= 0; i--)
      if(a->l[i] != 0)
         break;

   for(; i >= 0; i--)
   {
      if(init)
      {
         init = 0;
         fprintf(fp, BN_PRINT_FORMAT_I, a->l[i]);
      }
      else
         fprintf(fp, BN_PRINT_FORMAT, a->l[i]);
   }

   fputs((char *)post, fp);
}

bn_t *bn_read(FILE *fp, bn_t *dst)
{
   s8 *data = mem_alloc(dst->n);

   fread(data, sizeof(u8), dst->n, fp);
   bn_from_bin(dst, data, dst->n);
   mem_free(data);

   return dst;
}

bn_t *bn_write(FILE *fp, bn_t *num)
{
   u8 *data = mem_alloc(num->n);

   bn_to_bin(data, num);
   fwrite(data, 1, num->n, fp);
   mem_free(data);

   return num;
}
#endif

bn_t *bn_mul(bn_t *d, bn_t *a, bn_t *b)
{
   // D = A * b

   // Can a hold the result?
   assert(d->n >= (a->n * 2 + 1));

   bn_zero(d);

   for(int i = 0; i < a->n_limbs; i++)
   {
      ull_t S = 0;

      for(int j = 0; j < a->n_limbs; j++)
      {
         S += (ull_t)a->l[i] * b->l[j];
         d->l[i+j] += S;

         S >>= BN_LIMB_BITS;
      }

      d->l[i + a->n_limbs] = S;
   }

   return d;
}

bn_t *bn_mul_ui(bn_t *d, bn_t *a, ul_t b)
{
   // D = A * b
   ull_t S = 0;

   // Can a hold the result?
   assert(d->n >= (a->n + 1));

   bn_zero(d);

   for(int x = 0; x < a->n_limbs; x++)
   {
      S += (ull_t)a->l[x] * b;
      d->l[x] = S;

      S >>= BN_LIMB_BITS;
   }

   d->l[a->n_limbs] = S;

   return d;
}

bn_t *bn_divrem(bn_t *q, bn_t *r, bn_t *a, bn_t *b)
{
   bn_zero(q);
   bn_zero(r);

   for(int i = bn_maxbit(a); i >= 0 ;i--)
   {
      bn_lshift(r, 1);
      bn_lshift(q, 1);

      if(bn_getbit(a, i))
         bn_setbit(r, 0);

      if(bn_cmp(r, b) >= 0)
      {
         _bn_sub(r, r, b);
         bn_setbit(q, 0);
      }
   }

   return q;
}

bn_t *bn_rand(bn_t *a)
{
   int size = a->n;
   u8 *tmp = (u8 *)mem_alloc(size);

   #if defined(_WIN32) || defined(_MSC_VER)
      HCRYPTPROV hProvider;

      CryptAcquireContext(&hProvider, 0, 0, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT | CRYPT_SILENT);
      CryptGenRandom(hProvider, size, tmp);
   #else
      int fd = open("/dev/urandom", O_RDONLY, 0);
      read(fd, tmp, size);
      close(fd);
   #endif

   bn_zero(a);
   bn_from_bin(a, (s8 *)tmp, size);

   mem_free(tmp);

   return a;
}

// Generate random a \in [x, b - y].
bn_t *bn_rand_range(bn_t *a, int x, bn_t *b, int y)
{
   bn_t *t = bn_alloc(b->n);
   _bn_sub_ui(t, b, y);

   while(1)
   {
      bn_rand(a);

      if(bn_cmp_ui(a, x) <= 0)   // Check a < x
         continue;

      if(bn_cmp(a, t) > 0)       // Check a > b - y
         continue;

      break;
   }

   bn_free(t);

   return a;
}

bn_t *bn_to_mon(bn_t *a, bn_t *n)
{
   bn_t *at = bn_copy(bn_alloc_limbs(a->n_limbs + 1), a);
   bn_t *nt = bn_copy(bn_alloc_limbs(n->n_limbs + 1), n);

   // We can't loop bn_add here since bn_add calls bn_reduce which in turn calls bn_to_mon.
   // POOF, infinite recursion.
   for(int x = 0; x < BN_LIMB_BITS * a->n_limbs; x++)
   {
      _bn_add(at, at, at);
      bn_reduce(at, nt);
   }

   bn_copy(a, at);

   bn_free(at);
   bn_free(nt);

   return a;
}

bn_t *bn_from_mon(bn_t *a, bn_t *n)
{
   bn_t *t = bn_alloc(a->n);
   bn_set_ui(t, 1);

   bn_mon_mul(a, a, t, n);

   bn_free(t);

   return a;
}

bn_t *bn_mon_mul(bn_t *d, bn_t *a, bn_t *b, bn_t *n)
{
   ul_t q;
   ull_t r = (ull_t)1 << BN_LIMB_BITS;

   // The calculation of the mp value needs to be done only once.
   if(n->mp == 0)
      _bn_mon_init(n);

   // This num needs to be 4 digits bigger so we prevent overflows.
   // The mul_ui increases digit count by 1, and add's possibly increase
   // the count by one (each).
   bn_t *t = bn_alloc_limbs(n->n_limbs + 4);

   for(int x = 0; x < a->n_limbs; x++)
   {
      q = (t->l[0] + a->l[x] * b->l[0]) * n->mp;
      // q % r == q & (r-1), for r power of two
      q = q & (r - 1);

      _bn_mad_ui(t, n, q);
      _bn_mad_ui(t, b, a->l[x]);

      // Shift right by one limb
      _bn_rshift_limbs(t, 1);
   }

   if(bn_cmp(t, n) >= 0)
      _bn_sub(t, t, n);

   bn_copy(d, t);
   bn_free(t);

   return d;
}

// Montgomery reduction
bn_t *bn_mon_reduce(bn_t *a, bn_t *n)
{
   ull_t r = (ull_t)1 << BN_LIMB_BITS;
   ul_t mu;

   // The calculation of the mp value needs to be done only once.
   if(n->mp == 0)
      _bn_mon_init(n);

   bn_t *at = bn_copy(bn_alloc_limbs(a->n_limbs * 2 + 1), a);
   bn_t *tmp = bn_alloc_limbs(a->n_limbs * 2 + 1);

   for(int x = 0; x < a->n_limbs; x++)
   {
      mu = at->l[x] * n->mp % r;

      bn_mul_ui(tmp, n, mu);
      _bn_lshift_limbs(tmp, x);

      _bn_add(at, at, tmp);
   }

   _bn_rshift_limbs(at, a->n_limbs);

   bn_copy(a, at);

   bn_free(at);
   bn_free(tmp);

   return a;
}

// Square and multiply exponentiation
bn_t *bn_mon_pow_sm(bn_t *d, bn_t *a, bn_t *e, bn_t *n)
{
   // D = A**E % N
   bn_t *s = bn_copy(bn_alloc(a->n), a);
   bn_t *t = bn_copy(bn_alloc(d->n), d);

   bn_set_ui(t, 1);
   bn_to_mon(t, n);

   for(int i = 0; i <= bn_maxbit(e); i++)
   {
      if(bn_getbit(e, i))
         bn_mon_mul(t, t, s, n);

      bn_mon_mul(s, s, s, n);
   }

   bn_copy(d, t);

   bn_free(s);
   bn_free(t);

   return d;
}

// Square and multiply montgomery power ladder
bn_t *bn_mon_pow_ml(bn_t *d, bn_t *a, bn_t *e, bn_t *n)
{
   // D = A**E % N
   bn_t *r0 = bn_set_ui(bn_alloc(a->n), 1);
   bn_t *r1 = bn_copy(bn_alloc(a->n), a);

   bn_to_mon(r0, n);

   for(int i = e->n * 8; i >= 0; i--)
   {
      if(bn_getbit(e, i))
      {
         bn_mon_mul(r0, r0, r1, n);
         bn_mon_mul(r1, r1, r1, n);
      }
      else
      {
         bn_mon_mul(r1, r0, r1, n);
         bn_mon_mul(r0, r0, r0, n);
      }
   }

   bn_copy(d, r0);

   bn_free(r0);
   bn_free(r1);

   return d;
}

// Sliding-window exponentiation (HAC 14.83)
bn_t *bn_mon_pow_sw(bn_t *d, bn_t *a, bn_t *e, bn_t *n)
{
   // D = A**E % N
   bn_t *s = bn_copy(bn_alloc(a->n), a);
   bn_t *t = bn_copy(bn_alloc(d->n), d);

   // Select which window size to use
   int blen = BN_LIMB_BITS * n->n_limbs;
   int wsize = (blen > 671) ? 6 : (blen > 239) ? 5 : (blen >  79) ? 4 : (blen > 23) ? 3 : 1;

   //
   // Initialize the cache
   //
   bn_t *cache[1 << wsize];
   bn_set_ui(t, 1);
   bn_to_mon(t, n);

   // Initialize the cache first
   for(int i = 0; i < (1 << wsize); i++)
      cache[i] = bn_zero(bn_alloc(a->n));

   // 1st and 2nd elements are always the same
   bn_copy(cache[1], a);
   bn_mon_mul(cache[2], a, a, n);

   for(int i = 1; i < 1 << (wsize - 1); i++)
      bn_mon_mul(cache[2*i+1], cache[2*i-1], cache[2], n);

   // And iterate...
   for(int i = bn_maxbit(e); i >= 0;)
   {
      // In the 0-bit case, just square
      if(!bn_getbit(e, i))
      {
         bn_mon_mul(t, t, t, n);
         i--;
      }
      else
      {
         int num = 0;
         int sub = 0;
         int idx = 0;

         for(int j = 0; j < wsize; j++)
         {
            if(i - j < 0)
               break;

            int bit = bn_getbit(e, i - j);
            idx = (idx << 1) | bit;

            if(bit)
            {
               num = idx;
               sub = j + 1;
            }
         }

         // Square first
         for(int j = 0; j < sub; j++)
            bn_mon_mul(t, t, t, n);

         // ...then multiply with cache
         bn_mon_mul(t, t, cache[num], n);

         i -= sub;
      }
   }

   bn_copy(d, t);

   bn_free(s);
   bn_free(t);

   for(int i = 0; i < (1 << wsize); i++)
      bn_free(cache[i]);

   return d;
}

bn_t *bn_inv(bn_t *d, bn_t *a, bn_t *n)
{
   // D = A**-1 % N
   bn_t *q = bn_alloc_limbs(n->n_limbs);
   bn_t *r = bn_alloc_limbs(n->n_limbs);
   bn_t *m = bn_alloc_limbs(n->n_limbs);
   bn_t *t = bn_alloc_limbs(2*n->n_limbs);
   bn_t *b = bn_copy(bn_alloc_limbs(n->n_limbs), n);
   bn_t *nn = bn_copy(bn_alloc_limbs(2*n->n_limbs), n);

   bn_t *x = bn_alloc_limbs(n->n_limbs);
   bn_t *u = bn_alloc_limbs(n->n_limbs);

   bn_set_ui(x, 0);
   bn_set_ui(u, 1);

   while(bn_cmp_ui(a, 0) != 0)
   {
      bn_divrem(q, r, b, a);

      bn_mul(t, u, q);
      bn_reduce_slow(t, nn);
      bn_sub(m, x, t, n);

      bn_copy(b, a);
      bn_copy(a, r);
      bn_copy(x, u);
      bn_copy(u, m);
   }

   bn_copy(d, x);

   bn_free(q);
   bn_free(r);
   bn_free(m);
   bn_free(t);
   bn_free(b);
   bn_free(nn);
   bn_free(x);
   bn_free(u);

   return d;
}

bn_t *bn_mon_inv(bn_t *d, bn_t *a, bn_t *n)
{
   // D = A**-1 % N

   // Note: Only for prime modulus
   bn_t *t = bn_copy(bn_alloc(n->n), n);
   bn_t *s = bn_alloc(n->n);

   bn_set_ui(s, 2);

   _bn_sub(t, t, s);
   bn_mon_pow_sw(d, a, t, n);

   bn_free(t);
   bn_free(s);

   return d;
}

bn_t *bn_pow_mod(bn_t *d, bn_t *a, bn_t *e, bn_t *n)
{
   // D = A**E mod N
   bn_t *t = bn_copy(bn_alloc(a->n), a);
   bn_to_mon(t, n);

   bn_from_mon(bn_mon_pow_sw(d, t, e, n), n);
//   bn_from_mon(bn_mon_pow_ml(d, t, e, n), n);

   bn_free(t);

   return d;
}
