#if 0
CXX=icpc
CXXFLAGS=-O3 -g -lrt -march=native -std=c++11 -Wall -Wextra -Wshadow -Wno-unused-result
#CXXFLAGS=-O0 -g -std=c++11 -Wall -Wextra -Wshadow -Wno-unused-result
EXE=gen_simple.out gen_cube.out reference.out boruvka_simple.out boruvka_el.out boruvka_el_seq.out

all: ${EXE}

gen_simple.out: gen.o gen_simple.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

gen_cube.out: gen.o gen_cube.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

reference.out: gen.o reference.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_simple.out: gen.o boruvka_simple.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_el.out: gen.o boruvka_el.cpp makefile.h
	${CXX} ${CXXFLAGS} -fopenmp $^ -o $@

boruvka_el_seq.out: gen.o boruvka_el_seq.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

gen.o: gen.h gen.cpp makefile
	${CXX} ${CXXFLAGS} gen.cpp -c -o $@

clean:
	rm -rf *.o ${EXE}
#endif
