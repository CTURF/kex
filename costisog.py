#!/usr/bin/env python3

M = 1
S = 1

from costpoly import tree1,multiprod2,multiprod2_selfreciprocal,multieval_precompute,multieval_postcompute
import costs

x2 = S+S
x2DBL = M+M+M+M # but save 1M if affine
xDBL = x2+x2DBL
xADD = M+M+S+S+M+M

def isog(lmax,push,bsgs):
  assert lmax > 2
  assert lmax%2

  bs,gs = bsgs
  result = 0
  lbits = 6

  while lmax>>lbits: lbits += 2

  assert bs >= 0
  assert gs >= 0

  if bs == 0 or gs == 0:
    result += xDBL
    result += ((lmax-5)//2)*xADD
    result += ((lmax-3)//2)*2*M
    result += push*(2*M)
    result += ((lmax-3)//2)*4*M*push
  
  else: # velusqrt
    multiples = set()
    multiples.add(2)
    for j in range(3,2*bs+1,2): multiples.add(j)
    for j in range(4,lmax+1-4*bs*gs,2): multiples.add(j)
    multiples.add(2*bs)
    multiples.add(4*bs)
    for j in range(6*bs,2*bs*(2*gs+1),4*bs): multiples.add(j)
  
    result += len(multiples)*xADD # actually xADD or xDBL

    result += tree1[bs]
    result += gs*(3*S+4*M) # biquad_curve
    result += push*(3*S+M) # biquad_precompute_point
    result += push*gs*6*M # biquad_postcompute_point
    result += push*multiprod2[gs] # pushing point
    result += 2*multiprod2_selfreciprocal[gs] # pushing curve
    result += multieval_precompute[bs][gs*2+1] # reciprocal of root of product tree
    result += 2*(push+1)*multieval_postcompute[bs][gs*2+1] # scaled remainder tree
    result += 2*M*(push+1)*(bs-1) # accumulating multieval results
    result += 2*M*(1+2*push)*((lmax-1)//2-2*bs*gs) # stray points at the end
  
  result += push*(S+S+M+M) # final point evaluation
  result += 2*(S+M+M+(S+S+M)*(lbits//2-1)) # powpow8mod
  return result

def two_isog(prime_list, batching):
  assert prime_list is not None
  assert batching is not None
  result = 0
  loop_iterations = max(0, batching - 1)


  sqrt_cost = costs.sqrt(prime_list) - S
  result += S + sqrt_cost # Initial 2-isogeny that also moves to Montgomery+ curve
  result += loop_iterations * (2 * S + M + sqrt_cost)

  return result

def optimize(l,push):
  best,bestbsgs = isog(l,push,(0,0)),(0,0)

  # XXX: precompute more tree1 etc.; extend bs,gs limits
  for bs in range(2,33,2):
    for gs in range(1,2*bs+1):
      if 2*bs*gs > (l-1)//2: break
      if gs >= 32: break
      if bs > 3*gs: continue
      result = isog(l,push,(bs,gs))
      if result < best:
        best,bestbsgs = result,(bs,gs)

  return best,bestbsgs
