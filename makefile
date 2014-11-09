#if 0
CXX=icpc
CXXFLAGS=-O3 -g -lrt -lpthread -march=native -std=c++11 -Wall -Wextra -Wshadow -Wno-unused-result -fopenmp -DON_NUMA -pg
#CXXFLAGS=-O0 -g -lrt -lpthread -march=native -std=c++11 -Wall -Wextra -Wshadow -Wno-unused-result -fopenmp -DON_NUMA
EXE=gen_simple.out gen_cube.out turboboost_stub.out reference.out boruvka_simple.out boruvka_el.out boruvka_el_uma.out boruvka_el_seq.out boruvka_al_merge.out boruvka_al_copy.out

all: ${EXE}

gen_simple.out: gen.o gen_simple.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

gen_cube.out: gen.o gen_cube.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

turboboost_stub.out: gen.o turboboost_stub.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

reference.out: gen.o reference.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_simple.out: gen.o boruvka_simple.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_el.out: gen.o boruvka_el.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_el_uma.out: gen.o boruvka_el.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_el_seq.out: gen.o boruvka_el_seq.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_al_merge.out: gen.o boruvka_al_merge.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_al_copy.out: gen.o boruvka_al_copy.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

gen.o: gen.h gen.cpp makefile
	${CXX} ${CXXFLAGS} gen.cpp -c -o $@

clean:
	rm -rf *.o ${EXE}
#endif
