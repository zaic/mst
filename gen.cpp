#include "gen.h"
#include <cassert>
#include <cstdio>

const weight_t MAX_WEIGHT = 1.0;

vid_t vertexCount;
eid_t edgesCount;
eid_t *edgesIds;
Edge *edges;

bool EdgeDestCmp::operator()(const Edge& a, const Edge& b) const {
    if (a.dest != b.dest) return a.dest < b.dest;
    return a.weight < b.weight;
}

bool EdgeWeightCmp::operator()(const Edge& a, const Edge& b) const {
    return std::make_pair(a.weight, a.dest) < std::make_pair(b.weight, b.dest);
}

void readAll(char *filename) {
    FILE *f = fopen(filename, "rb");
    assert(f);

    fread(&vertexCount, sizeof(vid_t), 1, f);
    edgesIds = new eid_t[vertexCount + 1];

    fread(&edgesCount, sizeof(eid_t), 1, f);
    edges = new Edge[edgesCount];

    fread(edgesIds, sizeof(eid_t), vertexCount + 1, f);
    fread(edges, sizeof(Edge), edgesCount, f);

    fclose(f);
}
