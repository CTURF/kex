# Non-constant time python implementation of CSURF
# This was translated from the Magma implementation of CSURF at https://github.com/TDecru/CSURF/blob/master/csurf_512.m
# 
# This is used to generate the known answer tests for CTURF.

import random

import sympy
import math

class CSURF:

  def __init__(self, primes, primes_batchsize, primes_batchbound) -> None:
    self.primes_batchstart = []
    pos = 0
    for s in primes_batchsize:
      self.primes_batchstart.append(pos)
      pos += s

    self.primes_batchstop = []
    pos = 0
    for s in primes_batchsize:
      pos += s
      self.primes_batchstop.append(pos)

    self.primes = primes
    self.primes_batchsize = primes_batchsize
    self.primes_batchbound = primes_batchbound
    self.primes_batches = len(primes_batchbound)

    for prime in self.primes:
      assert sympy.isprime(prime)

    assert sum(self.primes_batchsize) == len(primes)
    assert len(self.primes_batchsize) == len(primes_batchbound)
    assert primes[0] == 2


    ells = list(self.primes)
    assert 2 in ells
    assert 2 in primes
    ells.remove(2)
    self.ells = tuple(ells)
    assert 2 not in self.ells
    assert 2 in self.primes
    self.p = 2**3 * 3 * math.prod(self.ells) - 1
    self.F = sympy.GF(self.p)
    self.ONE = self.F(1)
    self.ZERO = self.F(0)
    self.IS_SQUARE_EXP = (self.p - 1) // 2
    assert self.p == 3 * 4 * math.prod(self.primes) - 1
    assert sympy.isprime(self.p)
    assert len(self.ells) + 1 == len(self.primes)

  def is_square(self, a):
    # Legendre symbol calculation with Euler's criterion
    return pow(int(a), self.IS_SQUARE_EXP, self.p) <= 1

  def mod_division(self, a, b):
    return a * pow(int(b), -1, self.p)

  def step_zero_Montgomery_min(self, X1, Z1, X2, Z2, X3, Z3, A):
    X1 = X1
    Z1 = Z1
    X2 = X2
    Z2 = Z2
    X3 = X3
    Z3 = Z3
    A = A
    return (X2**2+Z2**2)**2, 4*X2*Z2*(X2**2+A*X2*Z2-Z2**2), 4*(X2*X3+Z2*Z3)**2*Z1, 4*(X2*Z3-Z2*X3)**2*X1


  def step_one_Montgomery_min(self, X1, Z1, X3, Z3, X2, Z2, A):
    return 4*(X2*X3+Z2*Z3)**2*Z1, 4*(X2*Z3-Z2*X3)**2*X1, (X2**2+Z2**2)**2, 4*X2*Z2*(X2**2+A*X2*Z2-Z2**2) 


  def scalar_multiplication_Montgomery_min(self, n, X1, Z1, A):
    X2 = self.ONE
    Z2 = self.ZERO
    X3 = X1
    Z3 = Z1
    nbits = bin(n)[2:]
    if Z1 == 0:
      return 'Error'
    else:
      for b in nbits:
        if b == '0':
          X2, Z2, X3, Z3 = self.step_zero_Montgomery_min(X1, Z1, X2, Z2, X3, Z3, A)
        else:
          X2, Z2, X3, Z3 = self.step_one_Montgomery_min(X1, Z1, X2, Z2, X3, Z3, A)
    return X2, Z2


  def differential_addition_Montgomery_min(self, X1, Z1, X2, Z2, X3, Z3, A):
    if X1 == 0 or Z1 == 0 or (X2,Z2) == (0,0) or (X3,Z3) == (0,0):
      return self.ZERO, self.ZERO
    else:
      return 4*(X2*X3+Z2*Z3)**2*Z1, 4*(X2*Z3-Z2*X3)**2*X1


  def double_point_Montgomery_min(self, X2, Z2, A):
    if Z2 == 0 or X2**3+A*X2**2+X2 == 0:
      return self.ZERO, self.ZERO
    else:
      return (X2**2+Z2**2)**2, 4*X2*Z2*(X2**2+A*X2*Z2-Z2**2)


  def act_with_ell_on_Montgomery_min(self, A, ell, xTors, xPush):
    XQ = xTors
    ZQ = self.ONE
    pi = XQ
    sigma = XQ + self.mod_division(1, XQ)
    fXPush = xPush*(XQ*xPush+1)**2
    fZPush = (xPush-XQ)**2
    if ell == 3:
      return pi**2*(A-6*sigma), fXPush, fZPush
    else:
      XQ, ZQ = self.double_point_Montgomery_min(XQ, ZQ, A)
    xQ = self.mod_division(XQ, ZQ)
    pi *= xQ
    sigma += xQ + self.mod_division(1, xQ)
    fXPush *= (xPush*xQ+1)**2
    fZPush *= (xPush - xQ)**2
    XPrev = xTors
    ZPrev = 1
    for i in range(3, ((ell-1) // 2) + 1):
      XTemp = XQ
      ZTemp = ZQ
      XQ, ZQ = self.differential_addition_Montgomery_min(XPrev, ZPrev, xTors, 1, XQ, ZQ, A)
      xQ = self.mod_division(XQ, ZQ)
      pi *= xQ
      sigma += xQ + self.mod_division(1, xQ)
      fXPush *= (xPush*xQ+1)**2
      fZPush *= (xPush - xQ)**2
      XPrev = XTemp
      ZPrev = ZTemp
    return pi**2*(A - 6*sigma), fXPush, fZPush


  def act_with_two_on_Montgomery_min(self, A, exp):
    A *= sympy.sign(exp)
    sq_exp = (self.p+1) // 4
    delta = pow(A**2 + 4, sq_exp)
    assert exp != 0
    A = self.mod_division(2*(A-3*delta), (A+delta))
    for i in range(1, abs(exp)):
      eps = pow(A**2 - 4, sq_exp)
      A = 2*(3 - A*(A - eps))
    eps = pow(A**2 - 4, sq_exp)
    scalefac = self.mod_division(eps*(eps + A), 2) ** sq_exp
    return self.mod_division(sympy.sign(exp) * (-A-3*eps), 2*scalefac)


  def csurf_action(self, private_key, A):
    A = self.F(A)
    A = self.act_with_two_on_Montgomery_min(A,private_key[0])
    private_key = private_key[1:]
    while any(i != 0 for i in private_key):
      xP = self.F(random.randint(0, self.p))

      twist = not self.is_square(xP**3+A*xP**2-xP)
      if twist:
        A = -A
        xP = -xP
      indices_ells_correct_sign = []
      k = 1
      for i in range(len(self.ells)):
        if private_key[i] > 0 and not twist:
          indices_ells_correct_sign.append(i)
          k *= self.ells[i]
        elif private_key[i] < 0 and twist:
          indices_ells_correct_sign.append(i)
          k *= self.ells[i]
      XQ, ZQ = self.scalar_multiplication_Montgomery_min((self.p + 1) // k, xP, 1, A)
      for i in indices_ells_correct_sign:
        if self.ells[i] == 3:
          if ZQ != 0:
            xQ = self.mod_division(XQ, ZQ)
            XR9, ZR9 = self.scalar_multiplication_Montgomery_min(k // 3, xQ, self.ONE, A)
            if ZR9 != 0:
              XR, ZR = self.scalar_multiplication_Montgomery_min(3, self.mod_division(XR9, ZR9), self.ONE, A)
              if ZR == 0:
                A, XQ, ZQ = self.act_with_ell_on_Montgomery_min(A, 3, self.mod_division(XR9, ZR9), xQ)
                assert i == 0
                if twist:
                  private_key[i] += 1
                else:
                  private_key[i] -= 1
              else:
                A, XQ, ZQ = self.act_with_ell_on_Montgomery_min(A, 3, self.mod_division(XR, ZR), xQ)
                assert i == 0
                if twist:
                  private_key[i] += 1
                else:
                  private_key[i] -= 1
        else:
          if ZQ != 0:
            xQ = self.mod_division(XQ, ZQ)
            ell = self.ells[i]
            XR, ZR = self.scalar_multiplication_Montgomery_min(k // ell, xQ, 1, A)
            if ZR != 0:
              A, XQ, ZQ = self.act_with_ell_on_Montgomery_min(A, ell, self.mod_division(XR, ZR), xQ)
              assert i != 0
              if twist:
                private_key[i] += 1
              else:
                private_key[i] -= 1
      if twist:
        A = -A
    return int(A)

  def is_valid_key(self, key):
    if len(key) != len(self.primes):
      return False
    if key[0] == 0:
        return False
    for i in range(self.primes_batches):
      batch = key[self.primes_batchstart[i] : self.primes_batchstop[i]]
      total = sum(abs(e) for e in batch)
      if total > self.primes_batchbound[i]:
        return False
    return True


  def csurf_private_keygen(self):
    # This does not generate uniformly random keys!!!
    # It is good enough for generating keys to test correctness, but like the rest of this implementation,
    # should not be used in a real-world setting.
    key = []

    two_isos = random.randint(-self.primes_batchbound[0], self.primes_batchbound[0])
    while two_isos == 0:
      two_isos = random.randint(-self.primes_batchbound[0], self.primes_batchbound[0])
    key.append(two_isos)

    for i in range(1, self.primes_batches):
      batch_key = [0] * self.primes_batchsize[i]
      for j in range(self.primes_batchbound[i]):
        index = random.randrange(self.primes_batchsize[i])
        if batch_key[index] == 0:
          batch_key[index] = random.randint(-1, 1)
        else:
          batch_key[index] += sympy.sign(batch_key[index])
      assert len(batch_key) == self.primes_batchsize[i]

      assert sum(abs(b) for b in batch_key) <= self.primes_batchbound[i]
      key += batch_key

    assert len(key) == len(self.primes)
    assert self.is_valid_key(key)
    return key

  def csurf_key_exchange(self):
    alice_private = self.csurf_private_keygen()
    bob_private = self.csurf_private_keygen()
    alice_public = self.csurf_action(alice_private, 0)
    bob_public = self.csurf_action(bob_private, 0)
    alice_bob = self.csurf_action(alice_private, bob_public)
    bob_alice = self.csurf_action(bob_private, alice_public)
    assert alice_bob == bob_alice

    alice_public = alice_public % self.p
    bob_public = bob_public % self.p
    alice_bob = alice_bob % self.p
    bob_alice = bob_alice % self.p
    assert alice_bob == bob_alice

    return (alice_private, bob_private), (alice_public, bob_public), (alice_bob, bob_alice)

if __name__ == '__main__':
  # Do one key exchange
  # Useful for playing around with this file
  p = {}
  p['cturf-505'] = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 373, 379, 383)
  p['cturf-1023'] = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 719, 727, 733, 751, 761)


  batchsize_512 = (1,2,3,4,4,5,5,6,7,7,8,7,7,7,1) # 15 batches
  batchbound_512 = (14,9,13,16,16,17,17,17,17,17,17,17,17,11,1)
  csurf_512 = CSURF(p['cturf-505'], batchsize_512, batchbound_512)
  priv, pub, shared = csurf_512.csurf_key_exchange()
  print('Alice private key:   ', priv[0])
  print('Bob private key:     ', priv[1])
  print('Alice public key:    ', pub[0])
  print('Bob public key:      ', pub[1])
  print('Alice shared secret: ', shared[0])
  print('Bob shared secret:   ', shared[1])
