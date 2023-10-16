CC=g++

funNames = Log2

src = $(patsubst %, %.cpp, $(funcNames))
obj = $(patsubst %, %, $(funcNames))

all: Log2 LogGenTables


%: %.cpp
	$(CC) -g -static  -frounding-math -fsignaling-nans -O3  $^ -lm -lmpfr -lgmp -o $@

clean:
	rm -rf Log2 LogGenTables $(obj) *~ \#*
