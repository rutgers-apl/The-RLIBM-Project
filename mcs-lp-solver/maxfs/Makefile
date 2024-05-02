all: maxfs clarkson_maxfs

maxfs: maxfs.cpp Utilities.cpp
	g++ -O2 -I/filer/tmp1/sn349/soplex-4.0.1/src/ -I../include/ -I../ -std=c++11 maxfs.cpp ../polyeval.cpp Utilities.cpp /filer/tmp1/sn349/soplex-4.0.1/build/lib/libsoplex.a -o maxfs -lgmp -lz -lm

clarkson_maxfs: clarkson_maxfs.cpp
	g++ -O2 -I/filer/tmp1/sn349/soplex-4.0.1/src/ -I../include/ -I../ -std=c++11 clarkson_maxfs.cpp ../polyeval.cpp /filer/tmp1/sn349/soplex-4.0.1/build/lib/libsoplex.a -o clarkson_maxfs -lgmp -lz -lm

clean:
	rm -rf *~ maxfs clarkson_maxfs
