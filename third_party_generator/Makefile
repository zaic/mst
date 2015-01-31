CC=icc
CXX=icpc
CFLAGS=-O3 -Wall -Wextra -Wshadow -std=gnu99 -openmp 
CXXFLAGS=-O3 -Wall -Wextra -Wshadow -std=c++11
LDFLAGS=-O3 -lrt 

TARGET=gen_RMAT gen_SSCA2

all: $(TARGET)

gen_RMAT: gen_RMAT.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_SSCA2: gen_SSCA2.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)

