#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "steps.h"
#include "elligator.h"
#include "csidh.h"
#include "primes.h"
#include "kat.h"

// Function for debugging
static void fp_print(fp const *n)
{
    printf("0x");
    for (size_t i = 8*UINTBIG_LIMBS-1; i < 8*UINTBIG_LIMBS; --i)
        printf("%02hhx", i[(unsigned char *) n->x.c]);
}

static bool is_key_valid(private_key const *priv)
{
  if (priv->e[0] == 0) {
    return false;
  }
  for (int i = 0; i < primes_batches; i++) {
    int batch_sum = 0;
    for (int j = primes_batchstart[i]; j < primes_batchstop[i]; j++) {
      int k = abs(priv->e[j]);
      batch_sum += k;
    }
    if (batch_sum > primes_batchbound[i]) {
      return false;
    }
  }
  return true;
}

static void test_iszero(void)
{
  printf("uintbig_iszero\n");
  fflush(stdout);

  uintbig u;
  unsigned char *x = (void *) &u;

  memset(x,0,sizeof u);
  assert(uintbig_iszero(&u));

  for (unsigned long long i = 0;i < 8*sizeof u;++i) {
    memset(x,0,sizeof u);
    x[i/8] = 1<<(i&7);
    assert(!uintbig_iszero(&u));
    for (unsigned long long j = 0;j < 8*sizeof u;++j) {
      memset(x,0,sizeof u);
      x[i/8] = 1<<(i&7);
      x[j/8] = 1<<(j&7);
      assert(!uintbig_iszero(&u));
    }
  }
}

static void test_sqrt(void)
{
  printf("fp_sqrt\n");
  fflush(stdout);

  for (long long loop = 0;loop < 1000;++loop) {
    fp x;
    fp xneg;
    fp x2;
    fp x2neg;
    fp_random(&x);
    fp_sq2(&x2,&x);
    fp_neg2(&xneg,&x);
    fp_neg2(&x2neg,&x2);
    assert(fp_sqrt(&x) != fp_sqrt(&xneg));
    assert(fp_sqrt(&x2));
    assert(!fp_sqrt(&x2neg));
  }
}

static void test_validate(void)
{
  printf("validate\n");
  fflush(stdout);

  private_key priv;
  public_key pub;

  pub.A = fp_2;
  assert(!validate(&pub));

  fp_neg1(&pub.A);
  assert(!validate(&pub));

  for (long long loop = 0;loop < 10;++loop) {
    fp_random(&pub.A);
    assert(!validate(&pub));
  }

  for (long long j = 0; j < 10; j++)
  {
      csidh_private(&priv);
      assert(csidh(&pub, &base, &priv));
      assert(validate(&pub));
  }

  uintbig_add3(&pub.A.x,&pub.A.x,&uintbig_p);
  assert(!validate(&pub));

  pub.A.x = uintbig_p;
  assert(!validate(&pub));
}

static void curve_setup(proj *A)
{
  private_key priv;
  public_key pub;
  csidh_private(&priv);
  assert(csidh(&pub,&base,&priv));
  fp_random(&A->z);
  fp_mul3(&A->x,&pub.A,&A->z);
}

static void kernel_setup(proj *K,long long k,const proj *A)
{
  uintbig cof = uintbig_1;
  uintbig_mul3_64(&cof, &cof, 4);
  uintbig_mul3_64(&cof, &cof, 3);
  for (long long i = 0;i < primes_num;++i)
    if (primes[i] != k)
      uintbig_mul3_64(&cof, &cof, primes[i]);
      
  for (;;) {
    proj P;
    fp_random(&P.x);
    fp_random(&P.z);
    xMUL_vartime(K,A,0,&P,&cof);
    if (memcmp(&K->z, &fp_0, sizeof(fp))) break;
  }
}

static void test_isog(void)
{
  printf("xISOG\n");
  fflush(stdout);

  proj A,K,P[10];
  proj A1,P1[10]; /* overwritten by xISOG */
  proj A2,P2[10]; /* overwritten by xISOG variants */

  // Curve setup
  private_key priv;
  public_key pub;
  csidh_private(&priv);
  assert(csidh(&pub,&base,&priv));
  A.x = pub.A;
  A.z = fp_1;
  xISOG_two(&A, 0, 0);

  for (long long t = 0;t <= 10;++t) {
    for (long long i = 1;i < primes_num;++i) {
      kernel_setup(&K,primes[i],&A);
      for (long long j = 0;j < t;++j) {
        fp_random(&P[j].x);
        fp_random(&P[j].z);
      }
      A1 = A;
      for (long long j = 0;j < t;++j)
        P1[j] = P[j];
      xISOG(&A1,P1,t,&K,primes[i]);

      for (long long matryoshka = 0;matryoshka < 10;++matryoshka) {
        long long ilower = (random() % i) + 1;
        long long iupper = i+(random()%(primes_num-i));
        assert (primes[ilower] <= primes[i]);
        assert (primes[i] <= primes[iupper]);
        A2 = A;
        for (long long j = 0;j < t;++j)
          P2[j] = P[j];
        xISOG_matryoshka(&A2,P2,t,&K,primes[i],primes[ilower],primes[iupper]);
        assert(proj_equal(&A2,&A1));
        for (long long j = 0;j < t;++j)
          assert(proj_equal(&P2[j],&P1[j]));
      }

      for (long long j = 0;j < t;++j) {
        A2 = A;
        P2[j] = P[j];
        xISOG_old(&A2,&P2[j],&K,primes[i]);
        assert(proj_equal(&A2,&A1));
        assert(proj_equal(&P2[j],&P1[j]));
      }
      A = A1;
    }
  }
}

static void test_elligator(void)
{
  printf("elligator\n");
  fflush(stdout);

  proj A;
  proj plusminus[10][2];
  fp X[10][2];

  for (long long curves = 0;curves < 10;++curves) {
    if (curves&1)
      curve_setup(&A);
    else {
      A.x = fp_0;
      fp_random(&A.z);
    }
    for (long long i = 0;i < 10;++i) {
      elligator(&plusminus[i][0],&plusminus[i][1],&A);
      fp_inv(&A.z); fp_mul2(&A.x,&A.z);

      for (long long pm = 0;pm < 2;++pm) {
        fp_inv(&plusminus[i][pm].z);
        fp_mul3(&X[i][pm],&plusminus[i][pm].x,&plusminus[i][pm].z);
        fp T;
        fp_add3(&T,&X[i][pm],&A.x);
        fp_mul2(&T,&X[i][pm]);
        fp_add2(&T,&fp_1);
        fp_mul2(&T,&X[i][pm]);
        assert(fp_sqrt(&T) == !pm);
      }
    }
    for (long long i = 0;i < 10;++i)
      for (long long j = 0;j < 10;++j) {
        assert(memcmp(&X[i][0],&X[j][1],sizeof(fp)));
        if (i != j) {
          assert(memcmp(&X[i][0],&X[j][0],sizeof(fp)));
          assert(memcmp(&X[i][1],&X[j][1],sizeof(fp)));
        }
      }
  }
}

/* compute x^3 + Ax^2 + x */
static void montgomery_rhs(fp *rhs, fp const *A, fp const *x)
{
    fp tmp;
    *rhs = *x;
    fp_sq1(rhs);
    fp_mul3(&tmp, A, x);
    fp_add2(rhs, &tmp);
    fp_add2(rhs, &fp_1);
    fp_mul2(rhs, x);
}

/* totally not constant-time. */
static void action_old(public_key *out, public_key const *in, private_key const *priv)
{
    proj A = {in->A, fp_1};

    xISOG_two_old(&A, priv->e[0], primes_batchbound[0]);

    uintbig k[2];
    uintbig_set(&k[0], 12); // Extra powers in p+1
    uintbig_set(&k[1], 12); // Extra powers in p+1

    uint8_t e[2][primes_num];

    for (int64_t i = 0; i < primes_num; ++i) {

        int8_t t = priv->e[i];

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            uintbig_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            uintbig_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            uintbig_mul3_64(&k[0], &k[0], primes[i]);
            uintbig_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_sqrt(&rhs);

        if (done[sign])
            continue;

        xMUL_vartime(&P, &A, 0, &P, &k[sign]);

        done[sign] = true;

        for (int64_t i = primes_num-1; i >= 1; --i) {  //changed loop direction

            if (e[sign][i]) {

                uintbig cof = uintbig_1;
                for (int64_t j = i - 1; j >= 0; --j)   //changed loop direction
                    if (e[sign][j])
                        uintbig_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL_vartime(&K, &A, 0, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG(&A, &P, 1, &K, primes[i]);

                    if (!--e[sign][i])
                        uintbig_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));
    to_montgomery_minus_curve(&A);
    fp_inv(&A.z);
    fp_mul2(&A.x, &A.z);
    A.z = fp_1;
    out->A = A.x;
}

static void test_nike(void)
{
  printf("nike\n");
  fflush(stdout);

  private_key priv_alice, priv_bob;
  public_key pub_alice, pub_bob;
  public_key shared_alice, shared_bob;
  public_key check;
  bool ok;

  for (long long bs = 0;bs <= 64;bs += 2) {
    for (long long gs = 0;;++gs) {
      if (!gs) if (bs) continue;
      if (!bs) if (gs) break;
      if (2*bs*gs > (primes[primes_num-1]-1)/2) break;
      if (gs > 4*bs) continue;
      if (bs > 4*gs) continue;

      printf("trying alice bs=%lld gs=%lld, bob bs=0 gs=0\n",bs,gs);
      fflush(stdout);

      steps_override(bs,gs);
    
      csidh_private(&priv_alice);
      ok = csidh(&pub_alice, &base, &priv_alice);
      assert(ok);
      action_old(&check,&base,&priv_alice);
      assert(!memcmp(&check,&pub_alice,sizeof pub_alice));
    
      steps_override(0,0);

      csidh_private(&priv_bob);
      ok = csidh(&pub_bob, &base, &priv_bob);
      assert(ok);
      action_old(&check,&base,&priv_bob);
      assert(!memcmp(&check,&pub_bob,sizeof pub_bob));
    
      ok = csidh(&shared_bob, &pub_alice, &priv_bob);
      assert(ok);
      action_old(&check,&pub_alice,&priv_bob);
      assert(!memcmp(&check,&shared_bob,sizeof shared_bob));

      steps_override(bs,gs);
    
      ok = csidh(&shared_alice, &pub_bob, &priv_alice);
      assert(ok);
      action_old(&check,&pub_bob,&priv_alice);
      assert(!memcmp(&check,&shared_alice,sizeof shared_alice));
  
      assert(!memcmp(&shared_alice, &shared_bob, sizeof(public_key)));
    }
  }
}

static void test_skgen(void)
{
  printf("skgen\n");
  fflush(stdout);
  for (int t = 0; t < 1000; t++) {
    private_key priv;
    csidh_private(&priv);
    assert(is_key_valid(&priv));
  }

  private_key invalid;
  for (int b = 0; b < primes_batches; b++) {
    for (int i = primes_batchstart[b]; i < primes_batchstop[b]; i++) {
      invalid.e[i] = primes_batchbound[i] + 1;
    }
  }
  assert(!is_key_valid(&invalid));
}

static void test_kat(void)
{
  printf("kat\n");
  fflush(stdout);

  for (int test = 0; test < known_answer_tests; test++) {
    private_key alice_priv;
    memcpy(alice_priv.e, alice_private_key_kats[test], sizeof alice_private_key_kats[test]);
    private_key bob_priv;
    memcpy(bob_priv.e, bob_private_key_kats[test], sizeof bob_private_key_kats[test]);

    assert(is_key_valid(&alice_priv));
    assert(is_key_valid(&bob_priv));

    public_key alice_pub;
    public_key bob_pub;
    public_key alice_shared_secret;
    public_key bob_shared_secret;

    public_key alice_pub_kat;
    public_key bob_pub_kat;
    public_key alice_shared_secret_kat;
    public_key bob_shared_secret_kat;

    alice_pub_kat.A = alice_public_key_kats[test];
    bob_pub_kat.A = bob_public_key_kats[test];
    alice_shared_secret_kat.A = alice_shared_secret_kats[test];
    bob_shared_secret_kat.A = bob_shared_secret_kats[test];

    assert(sizeof alice_shared_secret_kat == sizeof bob_shared_secret_kat);
    assert(!memcmp(&alice_shared_secret_kat, &bob_shared_secret_kat, sizeof alice_shared_secret_kat)); // Sanity check for KAT


    csidh(&alice_pub, &base, &alice_priv);
    csidh(&bob_pub, &base, &bob_priv);

    assert(sizeof alice_pub == sizeof alice_pub_kat);
    assert(sizeof bob_pub == sizeof bob_pub_kat);
    assert(!memcmp(&alice_pub, &alice_pub_kat, sizeof alice_pub));
    assert(!memcmp(&bob_pub, &bob_pub_kat, sizeof bob_pub));

    csidh(&alice_shared_secret, &bob_pub, &alice_priv);
    csidh(&bob_shared_secret, &alice_pub, &bob_priv);

    assert(sizeof alice_shared_secret == sizeof alice_shared_secret_kat);
    assert(sizeof bob_shared_secret == sizeof bob_shared_secret_kat);
    assert(sizeof alice_shared_secret == sizeof bob_shared_secret);
    assert(!memcmp(&alice_shared_secret, &alice_shared_secret_kat, sizeof alice_shared_secret));
    assert(!memcmp(&bob_shared_secret, &bob_shared_secret_kat, sizeof bob_shared_secret));
    assert(!memcmp(&alice_shared_secret, &bob_shared_secret, sizeof alice_shared_secret));


  }
}

int main()
{
  test_iszero();
  test_sqrt();
  test_elligator();
  test_skgen();
  test_validate();
  test_isog();
  test_kat();
  test_nike();
  return 0;
}
