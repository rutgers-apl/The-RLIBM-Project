all: gen_hull polygen 

gen_hull: gen_hull.cpp
	g++ -O3 gen_hull.cpp  -o gen_hull  -lgmp -lmpfr -lz -lm


polygen: polygen.cpp polyeval.cpp
	g++ -O3 -I/filer/tmp1/sn349/soplex-4.0.1/src/ -std=c++11 polygen.cpp polyeval.cpp /filer/tmp1/sn349/soplex-4.0.1/build/lib/libsoplex.a -o polygen -lgmp -lz -lm

clean:
	rm -rf *~ gen_hull polygen \#*
