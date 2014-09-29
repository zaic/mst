#include <stdio.h>
#include <stdlib.h>
#include <error.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <thread>
#include "defs.h"
#include "common.h"

char* inFilename;
char* outFilename;
uint32_t rootNumberToValidate;

void usage(int , char **argv)
{
    printf("Usage:\n");
    printf("    %s -in <input> [options]\n", argv[0]);
    printf("Options:\n");
    printf("    -in <input> -- input graph filename\n");
    printf("    -out <output> -- output filename (distances from root vertex). By default output is '<input>.v'\n");
    printf("    -root <root> -- root number for validation. Default is the first vertex\n");
    exit(1);
}

void init (int argc, char** argv, graph_t* )
{
    inFilename = outFilename = NULL;
    rootNumberToValidate = 0;
    for (int i = 1; i < argc; ++i) {
   		if (!strcmp(argv[i], "-in")) {
            inFilename = argv[++i];
        }
   		if (!strcmp(argv[i], "-out")) {
            outFilename = argv[++i];
        }
		if (!strcmp(argv[i], "-root")) {
			rootNumberToValidate = (int) atoi(argv[++i]);
        }
    }
    if (!inFilename) usage(argc, argv);
    if (!outFilename) {
        outFilename = (char *)malloc((strlen(inFilename) + 3) * sizeof(char));
        sprintf(outFilename, "%s.v", inFilename);
    }
}

void readGraph(graph_t *G, char *filename)
{
    edge_id_t arity;
    uint8_t align;
    FILE *F = fopen(filename, "rb");
    if (!F) error(EXIT_FAILURE, 0, "Error in opening file %s", filename);

    assert(fread(&G->n, sizeof(vertex_id_t), 1, F) == 1);
    G->scale = log(G->n) / log (2);   
    assert(fread(&arity, sizeof(edge_id_t), 1, F) == 1); 
    assert(fread(&G->directed, sizeof(bool), 1, F) == 1);
    assert(fread(&align, sizeof(uint8_t), 1, F) == 1);

    G->m = G->n * arity;

	G->rowsIndices = (edge_id_t *)malloc((G->n+1) * sizeof(edge_id_t));
    assert(G->rowsIndices);

	assert(fread(G->rowsIndices, sizeof(edge_id_t), G->n+1, F) == (G->n+1));
    
    G->endV = (vertex_id_t *)malloc((G->rowsIndices[G->n] + 123) * sizeof(vertex_id_t));
    assert(G->endV);
    memset(G->endV, 0, (G->rowsIndices[G->n] + 123) * sizeof(vertex_id_t));

    assert(fread(G->endV, sizeof(vertex_id_t), G->rowsIndices[G->n], F) == G->rowsIndices[G->n]);

    assert(fread(&G->nRoots, sizeof(uint32_t), 1, F) == 1);
    
    G->roots = (vertex_id_t *)malloc(G->nRoots * sizeof(vertex_id_t));
    assert(G->roots);
    G->numTraversedEdges = (edge_id_t *)malloc(G->nRoots * sizeof(edge_id_t));
    assert(G->numTraversedEdges);

    assert(fread(G->roots, sizeof(vertex_id_t), G->nRoots, F) == G->nRoots);
    assert(fread(G->numTraversedEdges, sizeof(edge_id_t), G->nRoots, F) == G->nRoots);

    G->weights = (weight_t *)malloc(G->m * sizeof(weight_t));
    assert(G->weights);

    assert(fread(G->weights, sizeof(weight_t), G->m, F) == G->m);

    fclose(F);
}

struct Edge {
    int32_t dest;
    double weight;
};

void saveGraph(graph_t *G, char *filename) {
    FILE *f = fopen(filename, "wb");
    assert(f);

    int32_t tmp;

    tmp = G->n;
    fwrite(&tmp, sizeof(int32_t), 1, f);
    tmp = G->m;
    fwrite(&tmp, sizeof(int32_t), 1, f);

    for (int i = 0; i <= G->n; ++i) {
        tmp = G->rowsIndices[i];
        fwrite(&tmp, sizeof(int32_t), 1, f);
    }
    for (int i = 0; i < G->m; ++i) {
        Edge q;
        q.dest = G->endV[i];
        q.weight = G->weights[i];
        fwrite(&q, sizeof(Edge), 1, f);
    }

    fclose(f);
}

void addReversedEdges(graph_t *G) {
    G->revRowsIndices = (edge_id_t *)malloc((G->n+1) * sizeof(edge_id_t));
    G->revStartV = (vertex_id_t *)malloc(G->rowsIndices[G->n] * sizeof(vertex_id_t));
    G->revWeights = (weight_t *)malloc(G->m * sizeof(weight_t));
    uint32_t *in_degree = (uint32_t*)calloc(G->n, sizeof(uint32_t));

    assert(G->revRowsIndices);
    assert(G->revStartV);
    assert(G->revWeights);
    assert(in_degree);

    for (vertex_id_t v = 0; v < G->n; ++v) {
        for(edge_id_t i = G->rowsIndices[v]; i < G->rowsIndices[v + 1]; ++i) {
            vertex_id_t u = G->endV[i];
            ++in_degree[u];
        }
    }

    G->revRowsIndices[0] = 0;
    for (vertex_id_t v = 0; v < G->n; ++v) {
        G->revRowsIndices[v + 1] = G->revRowsIndices[v] + in_degree[v];
    }

    for (vertex_id_t v = 0; v < G->n; ++v) {
        for(edge_id_t i = G->rowsIndices[v]; i < G->rowsIndices[v + 1]; ++i) {
            vertex_id_t u = G->endV[i];
            edge_id_t row_number = G->revRowsIndices[u] + (--in_degree[u]);
            G->revStartV[row_number] = v;
            G->revWeights[row_number] = G->weights[i];
        }
    }

    G->mRowsIndices = new edge_id_t*[getConcurrencyLevel()];
    G->mEndV = new vertex_id_t[G->m + 100];
    G->mWeights = new weight_t[G->m + 100];
    edge_id_t edges_sum = 0;
    vertex_id_t cur_vertex = 0;
    edge_id_t last_edge = 0;
    for (int thread_number = 0; thread_number < getConcurrencyLevel(); ++thread_number) {
        edge_id_t expected_sum = (G->m) * (thread_number + 1) / getConcurrencyLevel();
        for(int i=0; i<G->n; ++i) in_degree[i] = 0;
        vertex_id_t start_vertex = cur_vertex;
        for (; edges_sum < expected_sum; ++cur_vertex) {
            edges_sum += G->revRowsIndices[cur_vertex + 1] - G->revRowsIndices[cur_vertex];
        }

        for (vertex_id_t u = start_vertex; u < cur_vertex; ++u) {
            for (edge_id_t i = G->revRowsIndices[u]; i < G->revRowsIndices[u + 1]; ++i) {
                vertex_id_t v = G->revStartV[i];
                ++in_degree[v];
            }
        }

        //E(thread_number); E(edges_sum); Eo(last_edge); 
        G->mRowsIndices[thread_number] = new edge_id_t[G->n + 1];
        G->mRowsIndices[thread_number][0] = last_edge;
        last_edge = edges_sum;
        for (vertex_id_t v = 0; v < G->n; ++v) {
            G->mRowsIndices[thread_number][v + 1] = G->mRowsIndices[thread_number][v] + in_degree[v];
            // E(thread_number); E(v); Eo(G->mRowsIndices[thread_number][v + 1]);
        }

        for (vertex_id_t u = start_vertex; u < cur_vertex; ++u) {
            for (edge_id_t i = G->revRowsIndices[u]; i < G->revRowsIndices[u + 1]; ++i) {
                vertex_id_t v = G->revStartV[i];
                assert(in_degree[v] > 0);
                edge_id_t row_number = G->mRowsIndices[thread_number][v] + (--in_degree[v]);
                G->mEndV[row_number] = u;
                G->mWeights[row_number] = G->revWeights[i];
            }
        }
    }

    free(in_degree);
    free(G->revRowsIndices);
    free(G->revStartV);
    free(G->revWeights);
}

/* write distances from root vertex to each others to output file. -1 = infinity */
void writeDistance(char* filename, weight_t *dist, vertex_id_t n)
{
    FILE *F = fopen(filename, "wb");
    assert(fwrite(dist, sizeof(weight_t), n, F) == n);
    fclose(F);
}

void freeGraph(graph_t *G)
{
    free(G->rowsIndices);
    free(G->endV);
    free(G->weights);
    free(G->roots);
    free(G->numTraversedEdges);
}

int main(int argc, char **argv) 
{
    weight_t* dist;
    graph_t g;
    struct timespec start_ts, finish_ts;
    double *sssp_perf;

    /* initializing and reading the graph */
    init(argc, argv, &g); 
    readGraph(&g, inFilename);
    saveGraph(&g, outFilename);

    return 0;
}
