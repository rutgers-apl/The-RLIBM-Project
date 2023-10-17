all: polygen


polygen: polygen.h polygen.cpp
	@test $${SOPLEX_INCLUDE?Please set environment variable SOPLEX_INCLUDE}
	@test $${SOPLEX_LIB?Please set environment variable SOPLEX_LIB}
	g++ -O2 -I$(SOPLEX_INCLUDE) -std=c++11 polygen.cpp polyeval.cpp  $(SOPLEX_LIB) -o polygen -lgmp -lz -lm

clean:
	rm -f *~ polygen *.bak

