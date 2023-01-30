#!/bin/bash

START_SEARCH=2
END_SEARCH=25

for (( B = $START_SEARCH; B <= $END_SEARCH; B++))
do
    /usr/bin/time --output cturf_1024_batch_search.txt --append pypy3 greedy cturf-1023 256 $B 0 16 >> cturf_1024_batch_search.txt
done
