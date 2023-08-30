all: test

CXX = icpc
OPTS = -O3 -DNDEBUG -march=native -mtune=native -ffast-math -fopenmp
EIGEN_INC = -I /home/nfs_data/zhanggh/mytaco/learn-taco/zghshared/eigen-3.4.0
MKL_INSTALL_PATH = /home/eva_share/opt/intel/oneapi/
MKLLIBS = -L$(MKL_INSTALL_PATH)
LDMKL = -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread

test: test.o Makefile
	g++ -o test-coord -g ${EIGEN_INC} test.o -liomp5 -lpthread -lm -ldl ${LDMKL} -DCAP=${WCAP}

test.o: main.cpp ./utils/lib.h Makefile
	g++ -o test.o -g -c -fpermissive -std=c++17 -pthread ${EIGEN_INC} ${MKLLIBS} -DCAP=${WCAP} main.cpp

.PHONY: clean run

clean:
	rm -f test*
	rm -f *.o

run:
	./test