#include <stdio.h>
#include "bn.h"

int main()
{
   bn_t *N = bn_from_bin(bn_alloc(32), "\xa0\xe1\x53\x75\xe4\xe6\x94\xbb\x41\x8f\xdd\x18\xab\x4f\xae\x55\x42\xe9\x78\xd1\xdf\xef\x04\xbb\x90\x07\x2a\xa6\xc9\xc1\xf9\xc1", 32);
   bn_t *e = bn_from_bin(bn_alloc(32), "\x01\x00\x01", 3);
   bn_t *d = bn_from_bin(bn_alloc(32), "\x96\xb8\xe3\x6f\x45\x47\x3d\x3a\x7e\x3e\xe0\xfd\xd6\xa9\x6d\x02\x19\x97\xb2\xd0\x22\x9b\x69\xff\xc0\x52\x47\x60\x20\xb3\x89\x9d", 32);

   bn_t *r = bn_from_bin(bn_alloc(32), "\x01\x00\x01", 3);
   bn_t *rez = bn_alloc(32);

   bn_t *m1 = bn_from_bin(bn_alloc(32), "Hello World", 12);
   bn_t *m2 = bn_alloc(32);

   bn_t *c = bn_alloc(32);

  bn_inv(rez, d, N);
  bn_print(stdout, "REZ = ", rez, "\n");

   // RSA Encryption
   // c = m1 ^ e % N
   bn_pow_mod(c, m1, e, N);
   bn_print(stdout, "", d, "\n");

   // RSA Decryption
   // m2 = c ^ d % N
   bn_pow_mod(m2, c, d, N);
   bn_print(stdout, "", m2, "\n");

   s8 plain[32];
   bn_to_bin(plain, m2);

   printf("%s\n", &plain[20]);

   bn_free(N);
   bn_free(e);
   bn_free(d);
   bn_free(m1);
   bn_free(m2);
   bn_free(c);

   return 0;
}
