/*
* Copyright 2016 Luka Malisa <luka.malisha@gmail.com>
* Licensed under the terms of the GNU GPL, version 2
* http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
*/

#ifndef _BN_H_
#define _BN_H_

#if !defined(NAKED)
   #include <stdio.h>
   #include <stdint.h>
   #include <inttypes.h>

   #if !defined(_WIN32) && !defined(_MSC_VER)
      #include <unistd.h>
   #endif
#endif

#include "config.h"

#if BN_LIMB_SIZE == 8
   #define l_t       int8_t
   #define ul_t      uint8_t
   #define ll_t      int16_t
   #define ull_t     uint16_t
   #define BN_PRINT_FORMAT_I "%"PRIX8
   #define BN_PRINT_FORMAT "%02"PRIX8
   #define SWAP SWAP8
   #define BN_LIMB_BYTES 1
#elif BN_LIMB_SIZE == 16
   #define l_t       int16_t
   #define ul_t      uint16_t
   #define ll_t      int32_t
   #define ull_t     uint32_t
   #define BN_PRINT_FORMAT_I "%"PRIX16
   #define BN_PRINT_FORMAT "%04"PRIX16
   #define SWAP SWAP16
   #define BN_LIMB_BYTES 2
#elif BN_LIMB_SIZE == 32
   #define l_t       int32_t
   #define ul_t      uint32_t
   #define ll_t      int64_t
   #define ull_t     uint64_t
   #define BN_PRINT_FORMAT_I "%"PRIX32
   #define BN_PRINT_FORMAT "%08"PRIX32
   #define SWAP SWAP32
   #define BN_LIMB_BYTES 4
#elif BN_LIMB_SIZE == 64
   #define l_t       int64_t
   #define ul_t      uint64_t
   #define ll_t      __int128
   #define ull_t     unsigned __int128
   #define BN_PRINT_FORMAT_I "%"PRIX64
   #define BN_PRINT_FORMAT "%016"PRIX64
   #define SWAP SWAP64
   #define BN_LIMB_BYTES 8
#else
   #error "Are you a wizard?"
#endif

#if !defined(BN_ASSERT)
   #define assert(...) 
#endif

typedef char      s8;
typedef uint8_t   u8;
typedef int16_t   s16;
typedef uint16_t  u16;
typedef int32_t   s32;
typedef uint32_t  u32;
typedef int64_t   s64;
typedef uint64_t  u64;

#define BOOL int
#define TRUE 1
#define FALSE 0

#define BN_LIMB_BITS (BN_LIMB_BYTES * 8)
#define BN_MAX_DIGIT ((ul_t)-1)

/*! Comparison results. */
/*! Less. */
#define BN_CMP_L (-1)
/*! Greater. */
#define BN_CMP_G (1)
/*! Equal. */
#define BN_CMP_E (0)

/*! Generic helper macros. */
#define MAX(x, y) ((x >= y) ? x : y)
#define MIN(x, y) ((x >= y) ? y : x)

/*! Byte-order swapping. */
#define SWAP64(x) ((SWAP32(((x) & 0xffffffff)) << 32) | SWAP32(((x >> 32) & 0xffffffff)))
#define SWAP32(x) ((SWAP16(((x) & 0xffff)) << 16) | SWAP16(((x >> 16) & 0xffff)))
#define SWAP16(x) (((x << 8) & 0xff00) | ((x >> 8) & 0xff))
#define SWAP8(x) x

/*! Convert bits to number of limbs. */
#define BYTES_TO_LIMBS(x) ((x * 8) / BN_LIMB_BITS + ((x * 8) % BN_LIMB_BITS ? 1 : 0))
#define LIMBS_TO_BYTES(x) (x * BN_LIMB_BYTES)

/*! Bignum struct. */
typedef struct _bn_t
{
   /*! Length in bytes. */
   int n;
   /*! Length in limbs. */
   int n_limbs;
   /*! A temporary value needed by Montgomery multiplication. */
   ull_t mp;
   /*! The limbs. */
   ul_t *l;
} bn_t;

/*!
* \brief Returns the position of the highest-placed non-zero bit.
*/
int bn_maxbit(bn_t *a);

/*!
* \brief Returns the x-th bit.
*/
int bn_getbit(bn_t *a, int x);

/*!
* \brief Read bignum from char array (binary data).
*/
bn_t *bn_from_bin(bn_t *a, s8 *s, int len);

/*!
* \brief Export bignum to a char array.
*/
u8 *bn_to_bin(u8 *s, bn_t *a);

/*!
* \brief Read bignum from a string.
*/
bn_t *bn_from_str(bn_t *a, const s8 *s);

/*!
* \brief Zero out the bignum.
*/
bn_t *bn_zero(bn_t *num);

/*!
* \brief Allocate a new bignum of chosen size (in bytes).
*/
bn_t *bn_alloc(int size);

/*!
* \brief Allocate a new bignum of chosen size (in limbs).
*/
bn_t *bn_alloc_limbs(int size);

/*!
* \brief Copy the bignum.
*/bn_t *bn_copy(bn_t *a, bn_t *b);

/*!
* \brief Free the bignum.
*/
void bn_free(bn_t *a);

/*!
* \brief Set the bignum to an integer value.
*/
bn_t *bn_set_ui(bn_t *a, u64 val);

/*!
* \brief Add two bignums.
*/
bn_t *bn_add(bn_t *a, bn_t *b, bn_t *c, bn_t *n);

/*!
* \brief Add a bignum and an integer.
*/
bn_t *bn_add_ui(bn_t *d, bn_t *a, unsigned int b, bn_t *n);

/*!
* \brief Subtract two bignums.
*/
bn_t *bn_sub(bn_t *a, bn_t *b, bn_t *c, bn_t *n);

/*!
* \brief Subtract a bignum and an integer.
*/
bn_t *bn_sub_ui(bn_t *d, bn_t *a, unsigned int b, bn_t *n);

/*!
* \brief Compare two bignums. Returns: BN_CMP_L, BN_CMP_G or BN_CMP_E.
*/
int bn_cmp(bn_t *a, bn_t *b);

/*!
* \brief Compare a bignum to an integer. Returns: BN_CMP_L, BN_CMP_G or BN_CMP_E.
*/
int bn_cmp_ui(bn_t *a, ul_t b);

/*!
* \brief Return 1 if the bignum is zero, 0 otherwise.
*/
int bn_is_zero(bn_t *a);

/*!
* \brief a = a % n
*/
bn_t *bn_reduce(bn_t *a, bn_t *n);

/*!
* \brief Shift left by a number of bits.
*/
bn_t *bn_lshift(bn_t *a, int b);

/*!
* \brief Shift right by a number of bits.
*/
bn_t *bn_rshift(bn_t *a, int b);

/*!
* \brief Return the most significant bit.
*/
int bn_msb(bn_t *a);

/*!
* \brief Return the least significant bit.
*/
int bn_lsb(bn_t *a);

#if defined(BN_PRINT_FUNCS)
/*!
* \brief Prints the bignum in hexadecimal format.
*/
void bn_print(FILE *fp, const s8 *pre, bn_t *a, const s8 *post);

/*!
* \brief Read bignum from file stream.
*/
bn_t *bn_read(FILE *fp, bn_t *dst);

/*!
* \brief Write bignum to file stream.
*/
bn_t *bn_write(FILE *fp, bn_t *num);
#endif

/*!
* \brief Textbook multiplication of two bignums.
*/
bn_t *bn_mul(bn_t *d, bn_t *a, bn_t *b);

/*!
* \brief Multiply bignum with an integer.
*/
bn_t *bn_mul_ui(bn_t *d, bn_t *a, ul_t b);

/*!
* \brief Divide two bignums, while tracking both the quotient (q) as well as
*        and the remainder (r).
*/
bn_t *bn_divrem(bn_t *q, bn_t *r, bn_t *a, bn_t *b);

/*!
* \brief Fill the bignum with random data.
*/
bn_t *bn_rand(bn_t *a);

/*!
* \brief Generate random a \in [x, b - y].
*/
bn_t *bn_rand_range(bn_t *a, int x, bn_t *b, int y);

/*!
* \brief Convert to Montgomery form.
*/
bn_t *bn_to_mon(bn_t *a, bn_t *n);

/*!
* \brief Convert from Montgomery form.
*/
bn_t *bn_from_mon(bn_t *a, bn_t *n);

/*!
* \brief D = A * B % N
*/
bn_t *bn_mon_mul(bn_t *d, bn_t *a, bn_t *b, bn_t *n);

/*!
* \brief A = A % N
*/
bn_t *bn_mon_reduce(bn_t *a, bn_t *n);

/*!
* \brief D = A**E % N
*/
bn_t *bn_mon_pow(bn_t *d, bn_t *a, bn_t *e, bn_t *n);

/*!
* \brief D = A**-1 % N
*/
bn_t *bn_mon_inv(bn_t *d, bn_t *a, bn_t *n);

/*!
* \brief This is a helper function which does *NOT* take Montgomery form 
*        numbers. It will make conversions internally.
*
*        D = A**B % N
*/
bn_t *bn_pow_mod(bn_t *d, bn_t *a, bn_t *b, bn_t *n);

#endif // _BN_H_
