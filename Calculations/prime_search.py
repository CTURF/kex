#!/usr/bin/env python
# coding: utf-8

import math
import sys
import itertools

from sympy import nextprime, primerange, prime, isprime, binomial, factorint

def cturf_prime_search(max_bits : int, n : int, k : int, search_end : int):
    """
    Computes the best prime for CTURF according to some reasonable heuristics.
    Prime will be of the form p = 3 * 4 * (p_1 * p_2 * ... * p_n) - 1

    Args:
        desired_bits: Maximum number of bits for the prime.
        n: Desired number of component primes.
        k: Number of last primes to vary.
           The first (n - k) primes in (p_1, ..., p_n) will be the first (n - k) primes (i.e. 2, 3, 5 ...).
           Larger values of k will search over a larger search space and thus be more likely to find a better value.
           However, it will take a long time to compute this function with large values of k.
           Larger values increase memory usage.
        search_end: Maximum component prime to consider. Larger values increase memory usage.
    """

    # Computed constants
    max_size = 2**max_bits
    f = n - k # Number of fixed primes

    fixed_primes = [prime(i) for i in range(1, f + 1)]
    assert len(fixed_primes) == f
    largest_fixed = max(fixed_primes)
    initial_product = 3 * 4 * math.prod(fixed_primes)

    # Computing the search space
    search_start = largest_fixed + 1

    candidate_primes = list(primerange(search_start, search_end + 1))

    # Start-up info
    print('Searching for a prime p of at most', max_bits, 'bits')
    print('p + 1 should have', n, 'distinct prime factors')
    print('Fixed the first', f, 'smallest prime factors')
    print('Fixed primes:', fixed_primes)
    print('Searching for best candidates for the last', k, 'prime factors')
    print('Largest fixed prime is', largest_fixed)
    print('Searching for primes between', search_start, 'and', search_end)
    print('There are', len(candidate_primes), 'primes in this range')
    print('Search space has', binomial(len(candidate_primes), k), 'elements')
    print('Searching for tuples...')

    print('Making tuples generator...')
    component_prime_combinations = itertools.combinations(candidate_primes, k)
    print('Made tuples generator')

    print('Processing tuples...')
    possible_tuples = []
    for tup in component_prime_combinations:
        p = initial_product * math.prod(tup) - 1
        if p < max_size:
            possible_tuples.append(tup)
    print(f'Processed {len(possible_tuples)} tuples')
    
    print('Sorting tuples...')
    possible_tuples.sort(key=lambda tup : sum(math.sqrt(x) for x in tup))
    print('Sorted tuples')

    print('Searching for best tuple...')
    for tup in possible_tuples:
        p = initial_product * math.prod(tup) - 1
        assert p % 8 == 7

        if isprime(p):
            print()
            print(p)
            print(math.log(p, 2), 'bits')
            print()
            print('p + 1 prime factors:')
            prime_factors = sorted(list(factorint(p + 1).keys()))
            print(prime_factors)
            print('Success:', isprime(p))

            # This is extra info, potentially useful for deciding how large k should be.
            i = 0
            last_prime = 1
            for f in prime_factors:
                if f == nextprime(last_prime):
                    last_prime = f
                    i += 1
                else:
                    print('The first', i, 'primes are consecutive')
                    print(n - i, 'primes are not consecutive')
                    break
            return p
    print('No primes found')
    return None

if __name__ == '__main__':
    desired_bits = 512
    if len(sys.argv) > 1:
        desired_bits = int(sys.argv[1])
    
    n = 74
    if len(sys.argv) > 2:
        n = int(sys.argv[2])
    # 512: n = 74
    # 1024: n = 131

    k = 7
    if len(sys.argv) > 3:
        k = int(sys.argv[3])
    
    search_end = nextprime(800)
    if len(sys.argv) > 4:
        search_end = int(sys.argv[4])

    print(cturf_prime_search(desired_bits, n, k, search_end))
