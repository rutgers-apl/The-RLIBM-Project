all: test


test: sinf.c test.c
	gcc -g -Wall -std=c11  test.c sinf.c -mfma -o sin_test -lm

clean:
	rm -rf sin_test *~ *.bak
