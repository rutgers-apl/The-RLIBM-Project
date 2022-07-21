all: test


test: sinf.c test.c
	gcc -g -Wall -std=c11  test.c sinf.c -mfma -o sin_test -lm

exhaustive: exhaustive_test.c
	gcc -g -Wall -std=c11 sinf.c exhaustive_test.c -mfma -o exhaustive -lm

clean:
	rm -rf sin_test *~ *.bak exhaustive
