#!/bin/bash

START_SEARCH=2
END_SEARCH=20

for (( B = $START_SEARCH; B <= $END_SEARCH; B++))
do
    /usr/bin/time --output cturf_505_batch_search.txt --append pypy3 greedy cturf-505 252.0979 $B 0 16 >> cturf_505_batch_search.txt
done
