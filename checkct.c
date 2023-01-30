#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>

#include "csidh.h"
#include "crypto_classify.h"
#include "crypto_declassify.h"

void print_private_key(private_key const *priv)
{
    printf("[%" PRId8, priv->e[0]);
    for (int i = 1; i < primes_num; i++)
    {
        printf(", %" PRId8, priv->e[i]);
    }
    printf("]");
}

void fp_print(fp const *n)
{
    printf("0x");
    for (size_t i = 8*UINTBIG_LIMBS-1; i < 8*UINTBIG_LIMBS; --i)
        printf("%02hhx", i[(unsigned char *) n->x.c]);
}

int main()
{
  private_key priv_alice, priv_bob;
  public_key pub_alice, pub_bob;
  public_key shared_alice, shared_bob;

  csidh_private(&priv_alice);
  csidh_private(&priv_bob);

  //printf("Alice private key:    ");
  //print_private_key(&priv_alice);
  //printf("\nBob private key:      ");
  //print_private_key(&priv_bob);
  //printf("\n");

  // just in case csidh_private isn't using randombytes:
  crypto_classify(&priv_alice,sizeof priv_alice);
  crypto_classify(&priv_bob,sizeof priv_bob);

  action(&pub_alice, &base, &priv_alice);
  action(&pub_bob, &base, &priv_bob);

  //crypto_declassify(&pub_alice,sizeof pub_alice);
  //crypto_declassify(&pub_bob,sizeof pub_bob);
  //printf("Alice public key:       0x");
  //uint_print(&pub_alice.A.x);
  //printf("\nBob public key:         0x");
  //uint_print(&pub_bob.A.x);
  //printf("\n");

  // could declassify public keys at this point
  // but the action doesn't branch based on those

  action(&shared_alice, &pub_bob, &priv_alice);
  action(&shared_bob, &pub_alice, &priv_bob);

  // end of constant-time testing
  crypto_declassify(&shared_alice,sizeof shared_alice);
  crypto_declassify(&shared_bob,sizeof shared_bob);
  assert(!memcmp(&shared_alice, &shared_bob, sizeof(public_key)));

  // Print the shared secret in Montgomery form
  printf("Shared secret: ");
  fp_print(&shared_alice.A);
  printf("\n");

  return 0;
}
