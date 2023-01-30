#!/bin/bash

START_SEARCH=2
END_SEARCH=20

for (( B = $START_SEARCH; B <= $END_SEARCH; B++))
do
    /usr/bin/time --output cturf_512_batch_search.txt --append pypy3 greedy cturf-505 256 $B 0 16 >> cturf_512_batch_search.txt
done
