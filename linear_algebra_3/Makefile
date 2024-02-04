TARGETS=  basic_matmul_3-1 bit_matmul_3-2 faster_matmul_3-3 fast_transpose_matmul_3-4+5 strassen_boolean boolean_matinv strassen_gfp
TESTS=  basic_matmul_3-1.test bit_matmul_3-2.test faster_matmul_3-3.test fast_transpose_matmul_3-4+5.test

all:	$(TARGETS)

test: $(TESTS)

boolean_matinv.test: boolean_matinv
	cat Rand2.mat | ./$< > $@ ;	diff $@ Rand2_Inv.mat

%.test: %
	cat Rand.mat Rand.mat | ./$< > $@ ;	diff $@ Rand_Squared.mat

extern_matmul.o: fast_transpose_matmul_3-4+5.c
	$(CC) $(CFLAGS) -DNO_MAIN -c -o extern_matmul.o fast_transpose_matmul_3-4+5.c

strassen_boolean: strassen_boolean.c extern_matmul.o
	$(CC) $(CFLAGS) -o $@ $^

testclean:
	rm -f *.test

clean: testclean
	rm -f *.o $(TARGETS)
