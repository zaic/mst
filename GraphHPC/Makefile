# defines
CC=icc
MPICC=mpicc
CXX=icpc
MPICXX=mpicxx
CFLAGS= -g -O3 -Wall -std=gnu99 -openmp -I~/nsu/graph/mst
CXXFLAGS= -g -O3 -Wall -I~/nsu/graph/mst
LDFLAGS= -g -O3 -lrt 


TARGET = gen_valid_info validation gen_RMAT gen_SSCA2 mst_reference mst_reference_boruvka mst_reference_mpi mst mst_zaic

all: $(TARGET)

# your own implementation, executable must called mst
mst: main.o mst.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

mst_reference_mpi: main_mpi.mpi.o graph_tools.mpi.o gen_RMAT_mpi.mpi.o mst_reference_mpi.mpi.o gen_SSCA2_mpi.mpi.o 
	$(MPICXX) $^ -o $@ $(LDFLAGS)

mst_reference_boruvka: main.o mst_reference_boruvka.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

#mst_zaic: main.o graph_tools.o ../gen.o ../boruvka_el.o
#	${CXX} ${LDFLAGS} $^ -o $@ -lpthread -lrt -fopenmp

mst_zaic: main_zaic.cpp graph_tools.o ../gen.o
	${CXX} ${CXXFLAGS} -lpthread -lrt -fopenmp $^ -o $@


# reference implementation	
mst_reference: main.o mst_reference.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_RMAT: gen_RMAT.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_SSCA2: gen_SSCA2.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

validation: validation.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_valid_info: graph_tools.o mst_reference.o gen_valid_info.o
	$(CXX) $^ -o $@ $(LDFLAGS)

%.mpi.o: %.cpp
	$(MPICXX) -DUSE_MPI $(CXXFLAGS) -o $@ -c $<

%.mpi.o: %.c
	$(MPICC) -DUSE_MPI $(CFLAGS) -o $@ -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)

