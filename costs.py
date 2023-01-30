#!/usr/bin/env python3

from memoized import memoized
import chain
import costisog

M = 1
S = 1

x2 = S+S
x2DBL = M+M+M+M # but save 1M if affine
xDBL = x2+x2DBL
xADD = M+M+S+S+M+M
clear2 = xDBL

def daccost(n):
  return xDBL+xADD+n*xADD

@memoized
def bigprime(primes):
  p = 3 * 4
  for l in primes: p *= l
  p -= 1
  return p

@memoized
def inv(primes):
  p = bigprime(primes)
  invchain = chain.chain2(p-2)
  invchaincost = chain.cost2(invchain)
  return invchaincost[0]*M+invchaincost[1]*S

@memoized
def div(primes):
  return inv(primes)+M

@memoized
def sqrt(primes):
  p = bigprime(primes)
  sqrtchain = chain.chain2((p+1)//4)
  sqrtchaincost = chain.cost2(sqrtchain)
  return sqrtchaincost[0]*M+(sqrtchaincost[1]+1)*S

@memoized
def elligator(primes):
  return S+M+M+M+S+M+M+sqrt(primes)

def dac_search(target,r0,r1,r2,chain,chainlen,best,bestlen):
  if chainlen >= bestlen: return best,bestlen
  if r2 > target: return best,bestlen
  if r2<<(bestlen-1-chainlen) < target: return best,bestlen
  if r2 == target: return chain,chainlen
  chain *= 2
  chainlen += 1
  best,bestlen = dac_search(target,r0,r2,r0+r2,chain+1,chainlen,best,bestlen)
  best,bestlen = dac_search(target,r1,r2,r1+r2,chain,chainlen,best,bestlen)
  return best,bestlen

def dac(target):
  if target == 2:
    return 0, 0
  best = None
  bestlen = 0
  while best == None:
    bestlen += 1
    best,bestlen = dac_search(target,1,2,3,0,0,best,bestlen)
  return best,bestlen

@memoized
def daclen(primes):
  return [dac(primes[j])[1] for j in range(len(primes))]

@memoized
def batchstart(batchsize):
  B = len(batchsize)
  return [sum(batchsize[:j]) for j in range(B)]

@memoized
def batchstop(batchsize):
  B = len(batchsize)
  return [sum(batchsize[:j+1]) for j in range(B)]

@memoized
def maxdaclen(primes,batchsize):
  B = len(batchsize)
  return [max(daclen(primes)[j]
              for j in range(batchstart(batchsize)[b],
                             batchstop(batchsize)[b]))
          for b in range(B)]

@memoized
def maxdac(primes,batchsize,b):
  B = len(batchsize)
  M = maxdaclen(primes,batchsize)
  return daccost(M[b])

@memoized
def eachdac(primes,batchsize,b):
  B = len(batchsize)
  D = daclen(primes)
  return sum(daccost(D[j])
             for j in range(batchstart(batchsize)[b],
                            batchstop(batchsize)[b]))

@memoized
def bsgs(primes,batchsize,b):
  return costisog.optimize(primes[batchstart(batchsize)[b]],1)[1]

@memoized
def isog(push,primes,batchsize,b):
  bs,gs = bsgs(primes,batchsize,b)
  return costisog.isog(primes[batchstop(batchsize)[b]-1],push,(bs,gs))

def mults(x,primes,batchsize,batchbound):
  B = len(batchsize)
  assert batchsize[0] == 1
  mults = 0
  mults += div(primes)
  mults += x['elligator']*elligator(primes)
  mults += x['clear2']*clear2
  for b in range(B):
    if b == 0:
      mults += costisog.two_isog(primes, batchbound[0])
    else:
      mults += x['maxdac',b]*maxdac(primes,batchsize,b)
      mults += x['eachdac',b]*eachdac(primes,batchsize,b)
      mults += x['isog',0,b]*isog(0,primes,batchsize,b)
      mults += x['isog',1,b]*isog(1,primes,batchsize,b)
      mults += x['isog',2,b]*isog(2,primes,batchsize,b)
  return mults

def strstats(x,prefix,format,primes,batchsize,batchbound):
  B = len(batchsize)
  result = prefix
  result += 'mults %s ' % (format%mults(x,primes,batchsize,batchbound))
  result += 'AB %s ' % (format%x['AB'])
  result += 'elligator %s ' % (format%x['elligator'])
  result += 'clear2 %s ' % (format%x['clear2'])
  result += 'clear2 %s ' % (format%x['clear2'])
  result += 'isog0 %s ' % ' '.join(format%x['isog',0,b] for b in range(1, B))
  result += 'isog1 %s ' % ' '.join(format%x['isog',1,b] for b in range(1, B))
  result += 'isog2 %s ' % ' '.join(format%x['isog',2,b] for b in range(1, B))
  result += 'maxdac %s ' % ' '.join(format%x['maxdac',b] for b in range(B))
  result += 'eachdac %s ' % ' '.join(format%x['eachdac',b] for b in range(B))
  return result
