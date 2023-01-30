#!/usr/bin/env python3

import costs
import random
import sys

def trial(primes,batchsize,batchbound):
  B = len(batchsize)
  assert B == len(batchbound)
  assert sum(batchsize) == len(primes)
  batchstart = [sum(batchsize[:j]) for j in range(B)]
  batchstop = [sum(batchsize[:j+1]) for j in range(B)]

  todo = list(batchbound)
  todo[0] = 0 # 2 isogenies have a fixed cost, no randomness involved

  x = {}
  x['AB'] = 0
  x['elligator'] = 0
  x['clear2'] = 0
  for b in range(B):
    x['eachdac',b] = 0
    x['maxdac',b] = 0
    x['isog',0,b] = 0
    x['isog',1,b] = 0
    x['isog',2,b] = 0

  def clearnonselected():
    # clear 4 from P:
    x['clear2'] += 2

    # clear primes outside batches from P:
    for b in range(B):
      if outsidebatch[b]:
        x['eachdac',b] += 1

    # clear non-selected primes in batches from P:
    for t in range(targetlen):
      x['maxdac',target[t]] += batchsize[target[t]]-1

  while sum(todo) > 0:
    x['AB'] += 1

    outsidebatch = [todo[b] == 0 for b in range(B)]
    target = [b for b in range(B) if todo[b] > 0]
    targetlen = len(target)
    if targetlen > 3:
      target.reverse()
      # 7 6 5 4 3 2 1 0
      target = target[1:2]+target[3:]+target[2:3]+target[0:1]
      # 6 4 3 2 1 0 5 7

    # generate P0:
    x['elligator'] += 1

    for t in range(targetlen):
      if t == 0:
        clearnonselected() # from P0

      # kernel point:
      for u in range(t+1,targetlen):
        x['maxdac',target[u]] += 1

      success = random.randrange(primes[batchstart[target[t]]]) != 0

      if success:
        if t == targetlen-1:
          push = 0
        elif t == 0:
          push = 1
        else:
          push = 2
        if t == targetlen-2 and targetlen > 2:
          push = 1

        x['isog',push,target[t]] += 1
        todo[target[t]] -= 1

      if t == 0:
        x['elligator'] += 1 # generate P1
        clearnonselected() # from P1

      if t == targetlen-2 and targetlen > 2:
        x['maxdac',target[t]] += 1 # P0
      elif t < targetlen-1:
        x['maxdac',target[t]] += 2 # P0, P1

  assert x['AB'] >= max(batchbound)
  assert x['elligator'] == 2*x['AB']
  assert x['clear2'] == 4*x['AB']
  for b in range(1, B):
    assert x['isog',0,b]+x['isog',1,b]+x['isog',2,b] == batchbound[b]

  return x
