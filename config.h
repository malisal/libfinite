#ifndef _CONFIG_H_
#define _CONFIG_H_

/*! Size (in bits) of each limb */
#if !defined(BN_LIMB_SIZE)
   #define BN_LIMB_SIZE 32
#endif

/*! Custom memory functions, e.g., for embedded code. */
#define mem_alloc(x) malloc(x)
#define mem_free(x) free(x)

#endif // _CONFIG_H_
