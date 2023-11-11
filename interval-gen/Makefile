CC=g++

funNames = Log2 Log  Log10 Log10-new

src = $(patsubst %, %.cpp, $(funcNames))
obj = $(patsubst %, %, $(funcNames))

all: Log2 Log Log10 Log10-new Exp10 Exp Exp2


%: %.cpp
	$(CC) -O3 -static  -frounding-math -fsignaling-nans $^ -lm -lmpfr -lgmp -o $@

clean:
	rm -rf Log Log2 Log10 LogGenTables LnGenTables Log-256 $(obj) *~ \#* Exp Exp2 Exp10
