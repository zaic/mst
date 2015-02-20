#if 0
CXX=icpc
CXXFLAGS=-O3 -g -lrt -lpthread -march=native -std=c++11 -Wall -Wextra -Wshadow -Wno-unused-result -fopenmp -DON_NUMA -pg -DUSE_EDGE_STRUCT -DON_HOME -DUSE_REORDER_BFS -DUSE_SKIP_LOOPSq=5
EXE=vector_test.out graphviz.out gen_simple.out gen_cube.out turboboost_stub.out reference.out boruvka_simple.out boruvka_el.out boruvka_el_uma.out boruvka_el_seq.out boruvka_al_merge.out boruvka_al_copy_dfs.out boruvka_al_copy_bfs.out boruvka_fl_bfs_list.out boruvka_fl_bfs_vector.out boruvka_fl_pj.out boruvka_el_offset.out

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
	${CXX} ${CXXFLAGS} -DUSE_FAST_REDUCTION -DUSE_REDUCTION_SIMPLE -DUSE_SKIP_LAST_ITER $^ -o $@

boruvka_el_offset.out: gen.o boruvka_el_offset.cpp makefile.h
	${CXX} ${CXXFLAGS} -DUSE_SKIP_LOOPS=6 -DUSE_REDUCTION_TREE $^ -o $@

boruvka_el_uma.out: gen.o boruvka_el_uma.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_el_seq.out: gen.o boruvka_el_seq.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_al_merge.out: gen.o boruvka_al_merge.cpp makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@

boruvka_al_copy_dfs.out: gen.o boruvka_al_copy_dfs.cpp makefile.h
	${CXX} ${CXXFLAGS} -DUSE_SMALL_VECTOR $^ -o $@

boruvka_al_copy_bfs.out: gen.o boruvka_al_copy_bfs.cpp vector.h makefile.h
	${CXX} ${CXXFLAGS} -DAL_BFS__SMPCOPY -DUSE_SMALL_VECTOR -DUSE_BARANCING $^ -o $@

boruvka_fl_bfs_list.out: gen.o boruvka_fl_bfs_list.cpp vector.h makefile.h
	${CXX} ${CXXFLAGS} $^ -DUSE_BARANCING -DUSE_SMALL_VECTOR -o $@

boruvka_fl_bfs_vector.out: gen.o boruvka_fl_bfs_vector.cpp vector.h makefile.h
	${CXX} ${CXXFLAGS} -DUSE_BOUND -DUSE_SMALL_VECTOR $^ -o $@

boruvka_fl_pj.out: gen.o boruvka_fl_pj.cpp vector.h makefile.h
	${CXX} ${CXXFLAGS} -DUSE_BOUND $^ -o $@



gen.o: gen.h gen.cpp stat.h makefile
	${CXX} ${CXXFLAGS} gen.cpp -c -o $@

vector_test.out: vector.h vector.cpp makefile
	${CXX} ${CXXFLAGS} -DTEST vector.cpp -o $@

graphviz.out: graphviz.cpp gen.o makefile.h
	${CXX} ${CXXFLAGS} $^ -o $@



clean:
	rm -rf *.o ${EXE}
#endif
