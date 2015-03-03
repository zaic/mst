#pragma once

#include "GraphHPC/defs.h"
#include <cinttypes>
#include <utility>
#include <iostream>

#define E(x) { std::cerr << #x << " = " << (x) << "   "; }
#define Eo(x) { std::cerr << #x << " = " << (x) << std::endl; }
#define EO(x) Eo(x)


typedef int32_t vid_t;
typedef int64_t eid_t;
typedef double weight_t;
typedef uint64_t thread_vector_t;

typedef std::pair<eid_t, int>      pei;
typedef std::pair<eid_t, vid_t>    pev;
typedef std::pair<vid_t, int>      pvi;
typedef std::pair<vid_t, vid_t>    pvv;
typedef std::pair<weight_t, vid_t> pwv;

struct Edge {
    vid_t dest;
    int origOffset;
    weight_t weight;
};

struct EdgeDestCmp { bool operator()(const Edge& a, const Edge& b) const; };
struct EdgeWeightCmp { bool operator()(const Edge& a, const Edge& b) const; }; // TODO inline

const int kMaxThreads = 20;
const int kMaxIterations = 30;

extern vid_t vertexCount;
extern eid_t edgesCount;
extern eid_t *edgesIds;
extern Edge *edges;
extern int threadsCount;
extern int iterationNumber;

extern vid_t *graphEdgesTo;
extern weight_t *graphEdgesWeight;

//extern bool *componentEnd;
extern bool *isCoolEdge;

extern const weight_t MAX_WEIGHT;

void readAll(char *filename);
void convertAll(graph_t *G);

extern vid_t *rev; // original ordering, required to restore answer
extern vid_t *que; // original ordering, required to restore answer
//eid_t *origEdgeIds; // original data

void doReorderBfs();
void doReorderSimple();
static void doReorder() {
#if defined(USE_REORDER_BFS)
    doReorderBfs();
#elif defined(USE_REORDER_SIMPLE)
    doReorderSimple();
#else
    Eo("no reorder");
    rev = new vid_t[vertexCount];
    que = new vid_t[vertexCount];
#pragma omp parallel for
    for (vid_t i = 0; i < vertexCount; ++i) rev[i] = que[i] = i;
#endif
}

inline eid_t vertexDegree(vid_t v) {
    return edgesIds[v + 1] - edgesIds[v];
}

template<typename T>
inline T bit(int shift) {
    return T(1) << shift;
}

int64_t currentNanoTime();
int stickThisThreadToCore(int coreId);

struct RDTSC {
    double timers[1024][64]; // ToDo fix array location on NUMA
    double oneSecond;

    RDTSC();

    double get();
    void start(int timerId);
    double end(int timerId);
};
extern RDTSC rdtsc;
