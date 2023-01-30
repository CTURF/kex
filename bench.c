#include "cpucycles.h"
#include "csidh.h"
#include "fp.h"
#include "primes.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define KEYS 1000

// 1 to bench group action, 0 to bench validation
#define BENCH_GROUP_ACTION 0

private_key alice_private[KEYS];
private_key bob_private[KEYS];
public_key alice_public[KEYS];
public_key bob_public[KEYS];
public_key alice_shared_secret[KEYS];
public_key bob_shared_secret[KEYS];

int main(int argc, char **argv)
{
    long long iterations = 1;
    if (argc >= 2)
    {
        iterations = atoll(argv[1]);
    }

    #if BENCH_GROUP_ACTION
        long long alice_public_start_cycles[iterations][KEYS];
        long long alice_public_start_mulsq_count[iterations][KEYS];
        long long alice_public_end_cycles[iterations][KEYS];
        long long alice_public_end_mulsq_count[iterations][KEYS];

        long long bob_public_start_cycles[iterations][KEYS];
        long long bob_public_start_mulsq_count[iterations][KEYS];
        long long bob_public_end_cycles[iterations][KEYS];
        long long bob_public_end_mulsq_count[iterations][KEYS];

        long long alice_shared_secret_start_cycles[iterations][KEYS];
        long long alice_shared_secret_start_mulsq_count[iterations][KEYS];
        long long alice_shared_secret_end_cycles[iterations][KEYS];
        long long alice_shared_secret_end_mulsq_count[iterations][KEYS];

        long long bob_shared_secret_start_cycles[iterations][KEYS];
        long long bob_shared_secret_start_mulsq_count[iterations][KEYS];
        long long bob_shared_secret_end_cycles[iterations][KEYS];
        long long bob_shared_secret_end_mulsq_count[iterations][KEYS];

        for (int k = 0; k < KEYS; k++)
        {
            csidh_private(&alice_private[k]);
            csidh_private(&bob_private[k]);
        }

        assert(fp_mulsq_count == 0);
        for (int i = 0; i < iterations; i++)
        {
            for (int k = 0; k < KEYS; k++)
            {
                alice_public_start_mulsq_count[i][k] = fp_mulsq_count;
                alice_public_start_cycles[i][k] = cpucycles();
                action(&alice_public[k], &base, &alice_private[k]);
                alice_public_end_cycles[i][k] = cpucycles();
                alice_public_end_mulsq_count[i][k] = fp_mulsq_count;

                bob_public_start_mulsq_count[i][k] = fp_mulsq_count;
                bob_public_start_cycles[i][k] = cpucycles();
                action(&bob_public[k], &base, &bob_private[k]);
                bob_public_end_cycles[i][k] = cpucycles();
                bob_public_end_mulsq_count[i][k] = fp_mulsq_count;

                alice_shared_secret_start_mulsq_count[i][k] = fp_mulsq_count;
                alice_shared_secret_start_cycles[i][k] = cpucycles();
                action(&alice_shared_secret[k], &bob_public[k], &alice_private[k]);
                alice_shared_secret_end_mulsq_count[i][k] = fp_mulsq_count;
                alice_shared_secret_end_cycles[i][k] = cpucycles();

                bob_shared_secret_start_mulsq_count[i][k] = fp_mulsq_count;
                bob_shared_secret_start_cycles[i][k] = cpucycles();
                action(&bob_shared_secret[k], &alice_public[k], &bob_private[k]);
                bob_shared_secret_end_mulsq_count[i][k] = fp_mulsq_count;
                bob_shared_secret_end_cycles[i][k] = cpucycles();
            }
        }
        printf("GROUP_ACTION\n");
        printf("NUMBER_OF_KEYS %d\n", KEYS);
        printf("NUMBER_OF_ITERATIONS %lld\n", iterations);
        printf("Ran %lld key exchanges\n", iterations * KEYS);
        printf("Ran %lld group action evaluations\n", iterations * KEYS * 4);
        printf("TOTAL_MULS_AND_SQUARES %lld\n", fp_mulsq_count);
        printf("TOTAL_MULS %lld\n", fp_mulsq_count - fp_sq_count);
        printf("TOTAL_SQUARES %lld\n", fp_sq_count);
        printf("Format below: KEY_NUMBER ITERATION_NUMBER MULSQ_COUNT CYCLES\n");
        printf("--------------------------------------------------------------------------------\n");

        long long mulsqs, cycles;
        for (int k = 0; k < KEYS; k++)
        {
            for (int i = 0; i < iterations; i++)
            {
                mulsqs = alice_public_end_mulsq_count[i][k] - alice_public_start_mulsq_count[i][k];
                cycles = alice_public_end_cycles[i][k] - alice_public_start_cycles[i][k];
                printf("ALICE_PUBLIC %d %d %lld %lld\n", k, i, mulsqs, cycles);

                mulsqs = bob_public_end_mulsq_count[i][k] - bob_public_start_mulsq_count[i][k];
                cycles = bob_public_end_cycles[i][k] - bob_public_start_cycles[i][k];
                printf("BOB_PUBLIC %d %d %lld %lld\n", k, i, mulsqs, cycles);

                mulsqs = alice_shared_secret_end_mulsq_count[i][k] - alice_shared_secret_start_mulsq_count[i][k];
                cycles = alice_shared_secret_end_cycles[i][k] - alice_shared_secret_start_cycles[i][k];
                printf("ALICE_SHARED_SECRET %d %d %lld %lld\n", k, i, mulsqs, cycles);

                mulsqs = bob_shared_secret_end_mulsq_count[i][k] - bob_shared_secret_start_mulsq_count[i][k];
                cycles = bob_shared_secret_end_cycles[i][k] - bob_shared_secret_start_cycles[i][k];
                printf("BOB_SHARED_SECRET %d %d %lld %lld\n", k, i, mulsqs, cycles);
            }
        }
    #else
        long long start_cycles[iterations][KEYS];
        long long start_mulsq_count[iterations][KEYS];
        long long end_cycles[iterations][KEYS];
        long long end_mulsq_count[iterations][KEYS];

        assert(fp_mulsq_count == 0);
        for (int k = 0; k < KEYS; k++)
        {
            csidh_private(&alice_private[k]);
            action(&alice_public[k], &base, &alice_private[k]);
        }

        long long inital_mulsq_count = fp_mulsq_count;
        long long initial_sq_count = fp_sq_count;

        for (int i = 0; i < iterations; i++)
        {
            for (int k = 0; k < KEYS; k++)
            {
                start_mulsq_count[i][k] = fp_mulsq_count;
                start_cycles[i][k] = cpucycles();
                validate(&alice_public[k]);
                end_cycles[i][k] = cpucycles();
                end_mulsq_count[i][k] = fp_mulsq_count;
            }
        }

        long long total_mulsq_count = fp_mulsq_count - inital_mulsq_count;
        long long total_sq_count = fp_sq_count - initial_sq_count;
        long long total_mul_count = total_mulsq_count - total_sq_count;
        printf("VALIDATION\n");
        printf("NUMBER_OF_KEYS %d\n", KEYS);
        printf("NUMBER_OF_ITERATIONS %lld\n", iterations);
        printf("Ran %lld validations\n", iterations * KEYS);
        printf("TOTAL_MULS_AND_SQUARES %lld\n", total_mulsq_count);
        printf("TOTAL_MULS %lld\n", total_mul_count);
        printf("TOTAL_SQUARES %lld\n", total_sq_count);
        printf("Format below: KEY_NUMBER ITERATION_NUMBER MULSQ_COUNT CYCLES\n");
        printf("--------------------------------------------------------------------------------\n");

        long long mulsqs, cycles;
        for (int k = 0; k < KEYS; k++)
        {
            for (int i = 0; i < iterations; i++)
            {
                mulsqs = end_mulsq_count[i][k] - start_mulsq_count[i][k];
                cycles = end_cycles[i][k] - start_cycles[i][k];
                printf("VALIDATION %d %d %lld %lld\n", k, i, mulsqs, cycles);
            }
        }
    #endif

    return 0;
}
