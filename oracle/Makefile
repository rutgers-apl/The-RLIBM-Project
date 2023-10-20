CC = gcc

funcNames = Log2 Log Sin

src = $(patsubst %, %.c, $(funcNames))
obj = $(patsubst %, %, $(funcNames))

all: $(obj)

%: %.c
	$(CC) -static  -frounding-math -fsignaling-nans -O3 -std=c11 $^ -lm -lmpfr -lgmp -o $@


clean:
	rm -rf $(obj) *.dSYM *~
