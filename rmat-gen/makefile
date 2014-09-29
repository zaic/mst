CC=icc
CXX=icpc
CFLAGS=-O3 -w -std=c99 -fopenmp
CXXFLAGS=-O3 -w -std=c++11 -fopenmp
LDFLAGS=-lm -lrt -pthread -fopenmp

all: gen_RMAT_serial convert

gen_RMAT_serial: gen_RMAT_serial.o sssp_serial_reference.o
	$(CXX) $^ -o $@ ${LDFLAGS}

convert: convert.o
	${CXX} $^ -o $@ ${LDFLAGS}

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o gen_RMAT_serial convert
