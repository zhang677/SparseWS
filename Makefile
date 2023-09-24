all: test

CXX = icpc
OPTS = -O3 -DNDEBUG -march=native -mtune=native -ffast-math -fopenmp
EIGEN_INC = -I /home/zgh23/code/eigen-3.4.0
BOOST_INC = -I /home/zgh23/tools/boost_1_83_0
WCAP = 0

test: test.o Makefile
	g++ -o test-coord -g ${EIGEN_INC} ${BOOST_INC} test.o -liomp5 -lpthread -lm -ldl -DCAP=${WCAP}

test.o: main.cpp ./utils/lib.h Makefile
	g++ -o test.o -g -c -fpermissive -std=c++17 -pthread ${EIGEN_INC} ${BOOST_INC} -DCAP=${WCAP} main.cpp

.PHONY: clean run

clean:
	rm -f test*
	rm -f *.o

run:
	./test
