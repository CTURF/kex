# DO NOT EDIT! generated by ./autogen

SCC=gcc
CC=clang -O3 -Os -march=native -mtune=native -std=gnu99 -pedantic -Wall -Wextra -Wno-language-extension-token -fwrapv -DTIMECOP -DGETRANDOM

default: generic 505 512 1024

all: default timecop

generic: testrandom.out

505: costpoly505.out checkct505untuned bench505untuned test505.out \
checkct505mults bench505mults \
checkct505cycles bench505cycles \
ubench505 umults505

512: costpoly512.out checkct512untuned bench512untuned test512.out \
checkct512mults bench512mults \
checkct512cycles bench512cycles \
ubench512 umults512

1024: costpoly1024.out checkct1024untuned bench1024untuned test1024.out \
checkct1024mults bench1024mults \
checkct1024cycles bench1024cycles \
ubench1024 umults1024

timecop: \
checkct505untuned checkct505mults checkct505cycles \
checkct512untuned checkct512mults checkct512cycles \
checkct1024untuned checkct1024mults checkct1024cycles \

	valgrind ./checkct505untuned
	valgrind ./checkct505mults
	valgrind ./checkct505cycles
	valgrind ./checkct512untuned
	valgrind ./checkct512mults
	valgrind ./checkct512cycles
	valgrind ./checkct1024untuned
	valgrind ./checkct1024mults
	valgrind ./checkct1024cycles

# ----- benchmarks:

checkct505cycles: checkct.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_tunecycles505.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct505cycles checkct.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_tunecycles505.a libtest.a

checkct512cycles: checkct.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_tunecycles512.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct512cycles checkct.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_tunecycles512.a libtest.a

checkct1024cycles: checkct.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunecycles1024.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct1024cycles checkct.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunecycles1024.a libtest.a

checkct505mults: checkct.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_tunemults505.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct505mults checkct.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_tunemults505.a libtest.a

checkct512mults: checkct.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_tunemults512.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct512mults checkct.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_tunemults512.a libtest.a

checkct1024mults: checkct.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunemults1024.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct1024mults checkct.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunemults1024.a libtest.a

checkct505untuned: checkct.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct505untuned checkct.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

checkct512untuned: checkct.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct512untuned checkct.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

checkct1024untuned: checkct.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o checkct1024untuned checkct.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

bench505cycles: bench.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_tunecycles505.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench505cycles bench.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_tunecycles505.a libtest.a

bench512cycles: bench.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_tunecycles512.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench512cycles bench.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_tunecycles512.a libtest.a

bench1024cycles: bench.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunecycles1024.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench1024cycles bench.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunecycles1024.a libtest.a

bench505mults: bench.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_tunemults505.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench505mults bench.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_tunemults505.a libtest.a

bench512mults: bench.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_tunemults512.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench512mults bench.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_tunemults512.a libtest.a

bench1024mults: bench.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunemults1024.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench1024mults bench.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_tunemults1024.a libtest.a

bench505untuned: bench.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench505untuned bench.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

bench512untuned: bench.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench512untuned bench.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

bench1024untuned: bench.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o bench1024untuned bench.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

# ----- tests run automatically:

testrandom.out: testrandom
	[ -f testrandom.time ] && cat testrandom.time || :
	time ./testrandom > testrandom.out
	cmp testrandom.out testrandom.exp

testrandom: testrandom.o libhighctidh_base.a libtest.a
	$(CC) -o testrandom testrandom.o libhighctidh_base.a libtest.a

testrandom.o: testrandom.c random.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c testrandom.c

test505: test.c \
libhighctidh_505.a libkat505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o test505 test.c \
		libhighctidh_505.a libkat505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

test512: test.c \
libhighctidh_512.a libkat512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o test512 test.c \
		libhighctidh_512.a libkat512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

test1024: test.c \
libhighctidh_1024.a libkat1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o test1024 test.c \
		libhighctidh_1024.a libkat1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

test505.out: test505 test505.exp
	[ -f test505.time ] && cat test505.time || :
	time ./test505 > test505.out
	cmp test505.out test505.exp

test512.out: test512 test512.exp
	[ -f test512.time ] && cat test512.time || :
	time ./test512 > test512.out
	cmp test512.out test512.exp

test1024.out: test1024 test1024.exp
	[ -f test1024.time ] && cat test1024.time || :
	time ./test1024 > test1024.out
	cmp test1024.out test1024.exp

# ----- microbenchmarks (some run automatically):

costpoly505.out: costpoly505
	./costpoly505 > costpoly505.out
	cmp costpoly505.out costpoly.py

costpoly512.out: costpoly512
	./costpoly512 > costpoly512.out
	cmp costpoly512.out costpoly.py

costpoly1024.out: costpoly1024
	./costpoly1024 > costpoly1024.out
	cmp costpoly1024.out costpoly.py

costpoly505: costpoly.c \
libhighctidh_505.a libhighctidh_base.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o costpoly505 costpoly.c \
		libhighctidh_505.a libhighctidh_base.a libtest.a

costpoly512: costpoly.c \
libhighctidh_512.a libhighctidh_base.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o costpoly512 costpoly.c \
		libhighctidh_512.a libhighctidh_base.a libtest.a

costpoly1024: costpoly.c \
libhighctidh_1024.a libhighctidh_base.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o costpoly1024 costpoly.c \
		libhighctidh_1024.a libhighctidh_base.a libtest.a

umults505: umults.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o umults505 umults.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

umults512: umults.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o umults512 umults.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

umults1024: umults.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o umults1024 umults.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

ubench505: ubench.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o ubench505 ubench.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

ubench512: ubench.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o ubench512 ubench.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

ubench1024: ubench.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o ubench1024 ubench.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

# ----- libhighctidh_tune{mults,cycles}:

libhighctidh_tunemults505.a: steps_tunemults505.o Makefile
	ar cr libhighctidh_tunemults505.a steps_tunemults505.o
	ranlib libhighctidh_tunemults505.a

libhighctidh_tunemults512.a: steps_tunemults512.o Makefile
	ar cr libhighctidh_tunemults512.a steps_tunemults512.o
	ranlib libhighctidh_tunemults512.a

libhighctidh_tunemults1024.a: steps_tunemults1024.o Makefile
	ar cr libhighctidh_tunemults1024.a steps_tunemults1024.o
	ranlib libhighctidh_tunemults1024.a

steps_tunemults505.o: steps_tunemults505.c steps.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps_tunemults505.c

steps_tunemults512.o: steps_tunemults512.c steps.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps_tunemults512.c

steps_tunemults1024.o: steps_tunemults1024.c steps.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps_tunemults1024.c

steps_tunemults505.c: tunemults505.out tune2c Makefile
	./tune2c < tunemults505.out > steps_tunemults505.c

steps_tunemults512.c: tunemults512.out tune2c Makefile
	./tune2c < tunemults512.out > steps_tunemults512.c

steps_tunemults1024.c: tunemults1024.out tune2c Makefile
	./tune2c < tunemults1024.out > steps_tunemults1024.c

tunemults505.out: tunemults505 Makefile
	[ -f tunemults505.time ] && cat tunemults505.time || :
	time ./tunemults505 > tunemults505.out

tunemults512.out: tunemults512 Makefile
	[ -f tunemults512.time ] && cat tunemults512.time || :
	time ./tunemults512 > tunemults512.out

tunemults1024.out: tunemults1024 Makefile
	[ -f tunemults1024.time ] && cat tunemults1024.time || :
	time ./tunemults1024 > tunemults1024.out

tunemults505: tunemults.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o tunemults505 tunemults.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

tunemults512: tunemults.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o tunemults512 tunemults.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

tunemults1024: tunemults.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o tunemults1024 tunemults.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

libhighctidh_tunecycles505.a: steps_tunecycles505.o Makefile
	ar cr libhighctidh_tunecycles505.a steps_tunecycles505.o
	ranlib libhighctidh_tunecycles505.a

libhighctidh_tunecycles512.a: steps_tunecycles512.o Makefile
	ar cr libhighctidh_tunecycles512.a steps_tunecycles512.o
	ranlib libhighctidh_tunecycles512.a

libhighctidh_tunecycles1024.a: steps_tunecycles1024.o Makefile
	ar cr libhighctidh_tunecycles1024.a steps_tunecycles1024.o
	ranlib libhighctidh_tunecycles1024.a

steps_tunecycles505.o: steps_tunecycles505.c steps.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps_tunecycles505.c

steps_tunecycles512.o: steps_tunecycles512.c steps.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps_tunecycles512.c

steps_tunecycles1024.o: steps_tunecycles1024.c steps.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps_tunecycles1024.c

steps_tunecycles505.c: tunecycles505.out tune2c Makefile
	./tune2c < tunecycles505.out > steps_tunecycles505.c

steps_tunecycles512.c: tunecycles512.out tune2c Makefile
	./tune2c < tunecycles512.out > steps_tunecycles512.c

steps_tunecycles1024.c: tunecycles1024.out tune2c Makefile
	./tune2c < tunecycles1024.out > steps_tunecycles1024.c

tunecycles505.out: tunecycles505 Makefile
	[ -f tunecycles505.time ] && cat tunecycles505.time || :
	time ./tunecycles505 > tunecycles505.out

tunecycles512.out: tunecycles512 Makefile
	[ -f tunecycles512.time ] && cat tunecycles512.time || :
	time ./tunecycles512 > tunecycles512.out

tunecycles1024.out: tunecycles1024 Makefile
	[ -f tunecycles1024.time ] && cat tunecycles1024.time || :
	time ./tunecycles1024 > tunecycles1024.out

tunecycles505: tunecycles.c \
libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o tunecycles505 tunecycles.c \
		libhighctidh_505.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

tunecycles512: tunecycles.c \
libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o tunecycles512 tunecycles.c \
		libhighctidh_512.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

tunecycles1024: tunecycles.c \
libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o tunecycles1024 tunecycles.c \
		libhighctidh_1024.a libhighctidh_base.a libhighctidh_untuned.a libtest.a

# ----- libhighctidh_{505,512,1024}, size-dependent functions:

libhighctidh_505.a: uintbig505.o fp505.o fp_inv505.o fp_sqrt505.o primes505.o poly505.o mont505.o elligator505.o skgen505.o validate505.o csidh505.o Makefile
	ar cr libhighctidh_505.a uintbig505.o fp505.o fp_inv505.o fp_sqrt505.o primes505.o poly505.o mont505.o elligator505.o skgen505.o validate505.o csidh505.o
	ranlib libhighctidh_505.a

libhighctidh_512.a: uintbig512.o fp512.o fp_inv512.o fp_sqrt512.o primes512.o poly512.o mont512.o elligator512.o skgen512.o validate512.o csidh512.o Makefile
	ar cr libhighctidh_512.a uintbig512.o fp512.o fp_inv512.o fp_sqrt512.o primes512.o poly512.o mont512.o elligator512.o skgen512.o validate512.o csidh512.o
	ranlib libhighctidh_512.a

libhighctidh_1024.a: uintbig1024.o fp1024.o fp_inv1024.o fp_sqrt1024.o primes1024.o poly1024.o mont1024.o elligator1024.o skgen1024.o validate1024.o csidh1024.o Makefile
	ar cr libhighctidh_1024.a uintbig1024.o fp1024.o fp_inv1024.o fp_sqrt1024.o primes1024.o poly1024.o mont1024.o elligator1024.o skgen1024.o validate1024.o csidh1024.o
	ranlib libhighctidh_1024.a

csidh505.o: csidh.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h int64mask.h elligator.h elligator_namespace.h random.h random_namespace.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o csidh505.o -c csidh.c

csidh512.o: csidh.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h int64mask.h elligator.h elligator_namespace.h random.h random_namespace.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o csidh512.o -c csidh.c

csidh1024.o: csidh.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h int64mask.h elligator.h elligator_namespace.h random.h random_namespace.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o csidh1024.o -c csidh.c

validate505.o: validate.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o validate505.o -c validate.c

validate512.o: validate.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o validate512.o -c validate.c

validate1024.o: validate.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o validate1024.o -c validate.c

skgen505.o: skgen.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h random.h random_namespace.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o skgen505.o -c skgen.c

skgen512.o: skgen.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h random.h random_namespace.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o skgen512.o -c skgen.c

skgen1024.o: skgen.c csidh.h uintbig.h uintbig_namespace.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont.h proj.h mont_namespace.h primes.h primes_namespace.h csidh_namespace.h random.h random_namespace.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o skgen1024.o -c skgen.c

elligator505.o: elligator.c crypto_declassify.h elligator.h proj.h fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h elligator_namespace.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o elligator505.o -c elligator.c

elligator512.o: elligator.c crypto_declassify.h elligator.h proj.h fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h elligator_namespace.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o elligator512.o -c elligator.c

elligator1024.o: elligator.c crypto_declassify.h elligator.h proj.h fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h elligator_namespace.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o elligator1024.o -c elligator.c

mont505.o: mont.c steps.h steps_namespace.h mont.h uintbig.h uintbig_namespace.h proj.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont_namespace.h poly.h poly_namespace.h int64mask.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o mont505.o -c mont.c

mont512.o: mont.c steps.h steps_namespace.h mont.h uintbig.h uintbig_namespace.h proj.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont_namespace.h poly.h poly_namespace.h int64mask.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o mont512.o -c mont.c

mont1024.o: mont.c steps.h steps_namespace.h mont.h uintbig.h uintbig_namespace.h proj.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont_namespace.h poly.h poly_namespace.h int64mask.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o mont1024.o -c mont.c

poly505.o: poly.c mont.h uintbig.h uintbig_namespace.h proj.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont_namespace.h poly.h poly_namespace.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o poly505.o -c poly.c

poly512.o: poly.c mont.h uintbig.h uintbig_namespace.h proj.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont_namespace.h poly.h poly_namespace.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o poly512.o -c poly.c

poly1024.o: poly.c mont.h uintbig.h uintbig_namespace.h proj.h fp.h fp_namespace.h randombytes.h crypto_declassify.h mont_namespace.h poly.h poly_namespace.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-o poly1024.o -c poly.c

fp_inv505.o: fp_inv505.c fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h crypto_declassify.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp_inv505.c

fp_inv512.o: fp_inv512.c fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h crypto_declassify.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp_inv512.c

fp_inv1024.o: fp_inv1024.c fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h crypto_declassify.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp_inv1024.c

fp_sqrt505.o: fp_sqrt505.c fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h crypto_declassify.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp_sqrt505.c

fp_sqrt512.o: fp_sqrt512.c fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h crypto_declassify.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp_sqrt512.c

fp_sqrt1024.o: fp_sqrt1024.c fp.h uintbig.h uintbig_namespace.h fp_namespace.h randombytes.h crypto_declassify.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp_sqrt1024.c

primes505.o: primes505.c primes.h primes_namespace.h Makefile
	$(CC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c primes505.c

primes512.o: primes512.c primes.h primes_namespace.h Makefile
	$(CC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c primes512.c

primes1024.o: primes1024.c primes.h primes_namespace.h Makefile
	$(CC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c primes1024.c

fp505.o: fp505.S fp.h fp_namespace.h uintbig.h uintbig_namespace.h Makefile
	$(SCC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp505.S

fp512.o: fp512.S fp.h fp_namespace.h uintbig.h uintbig_namespace.h Makefile
	$(SCC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp512.S

fp1024.o: fp1024.S fp.h fp_namespace.h uintbig.h uintbig_namespace.h Makefile
	$(SCC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c fp1024.S

uintbig505.o: uintbig505.S uintbig.h uintbig_namespace.h Makefile
	$(SCC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c uintbig505.S

uintbig512.o: uintbig512.S uintbig.h uintbig_namespace.h Makefile
	$(SCC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c uintbig512.S

uintbig1024.o: uintbig1024.S uintbig.h uintbig_namespace.h Makefile
	$(SCC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c uintbig1024.S

kat505.o: kat.h kat_namespace.h primes.h fp.h Makefile
	$(SCC) -DBITS=505 -D'NAMESPACEBITS(x)=highctidh_505_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c kat505.c

kat512.o: kat.h kat_namespace.h primes.h fp.h Makefile
	$(SCC) -DBITS=512 -D'NAMESPACEBITS(x)=highctidh_512_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c kat512.c

kat1024.o: kat.h kat_namespace.h primes.h fp.h Makefile
	$(SCC) -DBITS=1024 -D'NAMESPACEBITS(x)=highctidh_1024_##x' -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c kat1024.c

# ----- libhighctidh_untuned, size-independent but normally replaced by tuned functions:

libhighctidh_untuned.a: steps_untuned.o Makefile
	ar cr libhighctidh_untuned.a steps_untuned.o
	ranlib libhighctidh_untuned.a

steps_untuned.o: steps_untuned.c steps.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps_untuned.c

# ----- libhighctidh_base, size-independent functions:

libhighctidh_base.a: steps.o random.o Makefile
	ar cr libhighctidh_base.a steps.o random.o
	ranlib libhighctidh_base.a

steps.o: steps.c steps.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c steps.c

random.o: random.c random.h int32_sort.h randombytes.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c random.c

# ----- functions that libhighctidh wants from a core crypto library:

libtest.a: crypto_classify.o crypto_declassify.o randombytes.o int32_sort.o Makefile
	ar cr libtest.a crypto_classify.o crypto_declassify.o randombytes.o int32_sort.o
	ranlib libtest.a

randombytes.o: randombytes.c randombytes.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c randombytes.c

int32_sort.o: int32_sort.c int32_sort.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c int32_sort.c

crypto_declassify.o: crypto_declassify.c crypto_declassify.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c crypto_declassify.c

crypto_classify.o: crypto_classify.c crypto_classify.h Makefile
	$(CC) -D'NAMESPACEGENERIC(x)=highctidh_##x' \
		-c crypto_classify.c

# ----- known answer tests libraries:
libkat505.a: kat505.o Makefile
	ar cr libkat505.a kat505.o
	ranlib libkat505.a

libkat512.a: kat512.o Makefile
	ar cr libkat512.a kat512.o
	ranlib libkat512.a

libkat1024.a: kat1024.o Makefile
	ar cr libkat1024.a kat1024.o
	ranlib libkat1024.a


.PHONY:

clean:
	rm bench505cycles bench505mults bench505untuned checkct505cycles checkct505mults checkct505untuned costpoly505 *505.o bench512cycles bench512mults bench512untuned checkct512cycles checkct512mults checkct512untuned costpoly512 *512.o bench1024cycles bench1024mults bench1024untuned checkct1024cycles checkct1024mults checkct1024untuned costpoly1024 *1024.o crypto_classify.o crypto_declassify.o int32_sort.o randombytes.o random.o steps.o steps_untuned.o testrandom.o
