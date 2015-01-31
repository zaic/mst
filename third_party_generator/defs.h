#ifndef __GRAPH_HPC_DEFS_H
#define __GRAPH_HPC_DEFS_H
#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <vector>

#define DEFAULT_ARITY 16
#define SMALL_COMPONENT_EDGES_THRESHOLD   2
#define FNAME_LEN   256
#define WEIGHT_ERROR 0.0001

using namespace std;

typedef uint32_t vertex_id_t;
typedef uint64_t edge_id_t;
typedef double weight_t;

/* The graph data structure*/
typedef struct
{
    /***
     The minimal graph repesentation consists of:
     n        -- the number of vertices
     m        -- the number of edges
     endV     -- an array of size 2*m (m for directed graphs) that stores the 
                 destination ID of an edge <src->dest>.
     rowsIndices -- an array of size n+1 that stores the degree 
                 (out-degree in case of directed graphs) and pointers to
                 the endV array. The degree of vertex i is given by 
                 rowsIndices[i+1]-rowsIndices[i], and the edges out of i are
                 stored in the contiguous block endV[rowsIndices[i] .. rowsIndices[i+1]-1].
     Vertices are ordered from 0 in our internal representation
     ***/
    int scale; /* log2 of vertices number */
    vertex_id_t n;
    edge_id_t m;
    int avg_vertex_degree; /* relation m / n */
    bool directed; 

    double a, b, c; /* RMAT graph parameters */

    bool permute_vertices;

    edge_id_t* rowsIndices;
    vertex_id_t* endV;
    
    /* Edge weights */
    weight_t* weights;

    weight_t min_weight, max_weight;

    /* Search parameters */
    uint32_t nRoots;
    vertex_id_t *roots;
    uint64_t *numTraversedEdges;

    /* Distributed version parameters */
    int nproc, rank;
    vertex_id_t local_n; /* local vertices number */
    edge_id_t local_m; /* local edges number */
    char filename[FNAME_LEN]; /*filename for output graph*/
} graph_t;

typedef struct
{
    vertex_id_t numTrees;
    edge_id_t numEdges;
    edge_id_t* p_edge_list;
    edge_id_t* edge_id;

} forest_t;


/* write graph to file */
void writeGraph(graph_t *G, char *filename);
void writeBinaryGraph(graph_t *G, char *filename);
void writeTextGraph_MPI(graph_t *G);
void writeBinaryGraph_MPI(graph_t *G);

/* read graph from file */
void readGraph(graph_t *G, char *filename);
void readGraph_rankFiles_MPI(graph_t *G, char *filename);
void readGraph_singleFile_MPI(graph_t *G, char *filename);

/* free graph memory */
void freeGraph(graph_t *G);


#ifdef __cplusplus
/* Single Source Shortest Path */
extern "C" void MST (graph_t *G, vector < vector < edge_id_t > > &trees);
/* initialize algorithm memory */
extern "C" void init_mst(graph_t *G);
extern "C" void finalize_mst();
extern "C" void gen_SSCA2_graph_MPI(graph_t *G);
extern "C" void gen_RMAT_graph_MPI(graph_t *G);
#else
void MST (graph_t *G, vector < vector < edge_id_t > > &trees);
void init_mst(graph_t *G);
void finalize_mst();

void gen_SSCA2_graph_MPI(graph_t *G);
void gen_RMAT_graph_MPI(graph_t *G);

#endif

//void sssp_mpi(vertex_id_t root, graph_t *G, weight_t *dist);

/* macroses for obtaining vertices distribution between nodes */
extern int rank, size;
extern uint32_t TotVertices; /*number of vertices in graph*/
extern int lgsize; /* log2 (number of processes) */
//#define SIZE_MUST_BE_A_POWER_OF_TWO
#ifdef SIZE_MUST_BE_A_POWER_OF_TWO 
#define MOD_SIZE(v) ((v) & ((1 << lgsize) - 1))
#define DIV_SIZE(v) ((v) >> lgsize)
#define MUL_SIZE(x) ((x) << lgsize)
#define VERTEX_OWNER(v) ((int)(MOD_SIZE(v)))
#define VERTEX_LOCAL(v) ((vertex_id_t)(DIV_SIZE(v)))
#define VERTEX_TO_GLOBAL(i) ((vertex_id_t)(MUL_SIZE((uint64_t)i) + (int)(rank)))
#else
#define MOD_SIZE(v) ((v) % size)
#define DIV_SIZE(v) ((v) / size)
#define MUL_SIZE(x) ((x) * size)
#define VERTEX_OWNER(v) ((int)(v/DIV_SIZE(TotVertices)))
#define VERTEX_LOCAL(v) ((vertex_id_t)(v%DIV_SIZE(TotVertices) ))
#define VERTEX_TO_GLOBAL(v_local) ((vertex_id_t)(DIV_SIZE(TotVertices)*rank + (int)(v_local)))
#endif 



#endif
