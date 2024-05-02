all: maxfs clarkson_maxfs

maxfs: maxfs.cpp Utilities.cpp
	@test $${SOPLEX_INCLUDE?Please set environment variable SOPLEX_INCLUDE}
	@test $${SOPLEX_LIB?Please set environment variable SOPLEX_LIB}	
	g++ -O2 -I$(SOPLEX_INCLUDE) -I../include/ -I../ -std=c++11 maxfs.cpp ../polyeval.cpp Utilities.cpp $(SOPLEX_LIB) -o maxfs -lgmp -lz -lm

clarkson_maxfs: clarkson_maxfs.cpp
	@test $${SOPLEX_INCLUDE?Please set environment variable SOPLEX_INCLUDE}
	@test $${SOPLEX_LIB?Please set environment variable SOPLEX_LIB}	
	g++ -O2 -I$(SOPLEX_INCLUDE) -I../include/ -I../ -std=c++11 clarkson_maxfs.cpp ../polyeval.cpp $(SOPLEX_LIB) -o clarkson_maxfs -lgmp -lz -lm

clean:
	rm -rf *~ maxfs clarkson_maxfs
