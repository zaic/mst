#ifndef __GRAPH_HPC_DEFS_H
#define __GRAPH_HPC_DEFS_H

#include <stdint.h>
#include <stdbool.h>

#define DEFAULT_ARITY 16
#define SMALL_COMPONENT_EDGES_THRESHOLD   2

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
    
    edge_id_t **mRowsIndices;
    vertex_id_t *mEndV;
    weight_t *mWeights;
    
    /* Edge weights */
    weight_t* weights;

    /* Reversed edages */
    edge_id_t *revRowsIndices;
    vertex_id_t *revStartV;
    weight_t *revWeights;

    weight_t min_weight, max_weight;

    /* Search parameters */
    uint32_t nRoots;
    vertex_id_t *roots;
    uint64_t *numTraversedEdges;

    /* Distributed version parameters */
    int nproc, rank;
    vertex_id_t local_n; /* local vertices number */
    edge_id_t local_m; /* local edges number */

} graph_t;

/* write graph to file */
void writeGraph(graph_t *G, char *filename);

/* read graph from file */
void readGraph(graph_t *G, char *filename);

/* free graph memory */
void freeGraph(graph_t *G);

/* Single Source Shortest Path */
#ifdef __cplusplus
extern "C" void sssp(vertex_id_t root, graph_t *G, weight_t *dist, uint64_t *traversed_edges);
#else
void sssp(vertex_id_t root, graph_t *G, weight_t *dist, uint64_t *traversed_edges);
#endif

void bfs(vertex_id_t root, graph_t *G, uint64_t *traversed_edges);

//void sssp_mpi(vertex_id_t root, graph_t *G, weight_t *dist);

extern int lgsize; /* log2 (number of processes) */
#define MOD_SIZE(v) ((v) & ((1 << lgsize) - 1))
#define DIV_SIZE(v) ((v) >> lgsize)
#define MUL_SIZE(x) ((x) << lgsize)

/* macroses for obtaining vertices distribution between nodes */
#define VERTEX_OWNER(v) ((int)(MOD_SIZE(v)))
#define VERTEX_LOCAL(v) ((vertex_id_t)(DIV_SIZE(v)))
#define VERTEX_TO_GLOBAL(r, i) ((vertex_id_t)(MUL_SIZE((vertex_id_t)i) + (int)(r)))



#endif
