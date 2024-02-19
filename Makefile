all: test

CXX = icpc
OPTS = -O3 -DNDEBUG -march=native -mtune=native -ffast-math -fopenmp
EIGEN_INC = -I /home/zgh23/code/eigen-3.4.0
BOOST_INC = -I /home/zgh23/tools/boost_1_83_0
SPLATT_INC = -I /home/zgh23/code/splatt/include
SPLATT_LIB = -L /home/zgh23/code/splatt/build/Linux-x86_64/lib
WCAP = 0
ALGPREFIX = SPWS
ALGNAME = COORDCF
ALG = ${ALGPREFIX}${ALGNAME}
CPPFILE = bench_taco_outer_cc_noT
EXEPREFIX = bench-cc-noT

test: test.o Makefile
	g++ -o ${EXEPREFIX}-${ALGNAME} -g -fopenmp -DCAP=${WCAP} -D${ALG} ${EIGEN_INC} ${BOOST_INC} test.o ${SPLATT_INC}  -Wl,-Bstatic ${SPLATT_LIB} -lsplatt -Wl,-Bdynamic -liomp5 -lpthread -lm -ldl -L/home/zgh23/code/OpenBLAS/lib/lib -lopenblas 

test.o: ./utils/lib.h Makefile
	g++ -o test.o -g -c -fpermissive -std=c++17 -pthread -fopenmp ${EIGEN_INC} ${BOOST_INC} ${SPLATT_INC} -DCAP=${WCAP} -D${ALG} ${CPPFILE}.cpp

.PHONY: clean run

clean:
	rm -f *.o

run:
	./test
