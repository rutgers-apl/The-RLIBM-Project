all: gen_hull polygen 

gen_hull: gen_hull.cpp
	g++ -O3 gen_hull.cpp  -o gen_hull  -lgmp -lmpfr -lz -lm


polygen: polygen.cpp polyeval.cpp
	@test $${SOPLEX_INCLUDE?Please set environment variable SOPLEX_INCLUDE}
	@test $${SOPLEX_LIB?Please set environment variable SOPLEX_LIB}
	g++ -O3 -I$(SOPLEX_INCLUDE) -std=c++11 polygen.cpp polyeval.cpp $(SOPLEX_LIB) -o polygen -lgmp -lz -lm

clean:
	rm -rf *~ gen_hull polygen \#*
