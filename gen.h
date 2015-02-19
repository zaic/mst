#pragma once

#include <cinttypes>
#include <utility>
#include <iostream>

#define E(x) { std::cerr << #x << " = " << (x) << "   "; }
#define Eo(x) { std::cerr << #x << " = " << (x) << std::endl; }
#define EO(x) Eo(x)


typedef int32_t vid_t;
typedef int64_t eid_t;
typedef double weight_t;

typedef std::pair<vid_t, vid_t> pvv;
typedef std::pair<eid_t, int> pei;
typedef std::pair<weight_t, vid_t> pwv;
typedef std::pair<vid_t, int> pvi;

struct Edge {
    vid_t dest;
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

extern bool *componentEnd;

extern const weight_t MAX_WEIGHT;

void readAll(char *filename);

void doReorderBfs();
void doReorderSimple();
static void doReorder() {
#if defined(USE_REORDER_BFS)
    doReorderBfs();
#elif defined(USE_REDUCTION_SIMPLE)
    doReorderSimple();
#else
    // do nothing
#endif
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
