#include "gen.h"
#include <numeric>
#include <set>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <algorithm>
#include <vector>
#ifdef __clang__
#include "omp.h"
#else
#include <omp.h>
#endif

const weight_t MAX_WEIGHT = 1.0;

//
// Global data
//
vid_t vertexCount;
eid_t edgesCount;
eid_t *edgesIds;
Edge *edges;
bool *isCoolEdge;
int threadsCount;
int iterationNumber;
vid_t *rev;
vid_t *que;

//
// BFS-reorder specific variable
//
bool *componentEnd;

//
// Edge comparators
//
bool EdgeDestCmp::operator()(const Edge& a, const Edge& b) const {
    if (a.dest != b.dest) return a.dest < b.dest;
    return a.weight < b.weight;
}

bool EdgeWeightCmp::operator()(const Edge& a, const Edge& b) const {
    return a.weight < b.weight;
}

//
// Reorder functions
//
void doReorderBfs() {
    stickThisThreadToCore(0);

    //componentEnd = new bool[vertexCount]();
    bool *visit = new bool[vertexCount]();
    que = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount));
    rev = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount));
    memset(que, 0xc0, sizeof(vid_t) * vertexCount);

    Eo(vertexDegree(0)); // TODO read below
    std::vector<pev> largeVertexes;
    vid_t cntr = 0;
#if 1
    for (vid_t v = 0; v < vertexCount; ++v) if (vertexDegree(v) > 99) {
        ++cntr;
        visit[v] = true;
        //que[v] = v;
        largeVertexes.push_back(pev(vertexDegree(v), v));
    }
    sort(largeVertexes.begin(), largeVertexes.end());
    const vid_t threadOffset = vertexCount / threadsCount;
    for (vid_t i = 0; i < largeVertexes.size(); i += threadsCount) {
        const vid_t from = i;
        const vid_t to = std::min<vid_t>(largeVertexes.size(), i + threadsCount);
        std::random_shuffle(largeVertexes.begin() + from, largeVertexes.begin() + to);
    }
    for (vid_t i = 0; i < largeVertexes.size(); ++i) {
        int toThread = i % threadsCount;
        vid_t toPos = i / threadsCount;
        vid_t pos = toThread * threadOffset + toPos;
        assert(que[pos] < 0);
        assert(pos < vertexCount);
        que[pos] = largeVertexes[i].second;
    }
#endif
    Eo(cntr);

    std::vector<vid_t> vertexByDegree(vertexCount, 0);
    std::iota(vertexByDegree.begin(), vertexByDegree.end(), 0);
    std::sort(vertexByDegree.begin(), vertexByDegree.end(), [&](vid_t a, vid_t b) {
            return vertexDegree(a) < vertexDegree(b);
            });

    vid_t fr = 0, bc = 0;
    for (vid_t ii = 0; ii < vertexCount; ++ii) { // start from vertex with lower degree
        const vid_t i = vertexByDegree[ii];
        if (visit[i]) continue;
        while (bc < vertexCount && que[bc] >= 0) ++bc;
        que[bc++] = i;
        visit[i] = true;
        while (fr < bc) {
            const vid_t v = que[fr++];
            for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e) {
                const vid_t u = edges[e].dest;
                if (visit[u]) continue;
                visit[u] = true;
                while (bc < vertexCount && que[bc] >= 0) ++bc;
                que[bc++] = u;
            }
        }
        //componentEnd[bc - 1] = true;
    }
    while (bc < vertexCount && que[bc] >= 0) ++bc;
    E(largeVertexes.size()); Eo(double(largeVertexes.size()) / vertexCount);
    E(fr); E(bc); Eo(vertexCount);
    assert(bc == vertexCount);
    delete[] visit;

#pragma omp parallel for
    for (vid_t i = 0; i < vertexCount; ++i)
        rev[que[i]] = i;

    eid_t *nextEdgesIds = static_cast<eid_t*>(malloc(sizeof(eid_t) * (vertexCount + 1)));
    Edge *nextEdges = static_cast<Edge*>(malloc(sizeof(Edge) * edgesCount));
    eid_t *sumEdges = new eid_t[threadsCount];
    memset(sumEdges, 0, sizeof(eid_t) * threadsCount);
    nextEdgesIds[0] = 0;
    for (int i = 0; i < threadsCount; ++i) {
        stickThisThreadToCore(i);
        const vid_t vertexBegin = int64_t(vertexCount) * (i + 0) / threadsCount;
        const vid_t vertexEnd   = int64_t(vertexCount) * (i + 1) / threadsCount;
        for (vid_t v = vertexBegin; v < vertexEnd; ++v) {
            const vid_t nextv = que[v];
            nextEdgesIds[v + 1] = nextEdgesIds[v] + edgesIds[nextv + 1] - edgesIds[nextv];
            sumEdges[i] += edgesIds[nextv + 1] - edgesIds[nextv];
        }
    }
    std::sort(sumEdges, sumEdges + threadsCount);
    double disbalanceFactor = double(sumEdges[threadsCount - 1]) / sumEdges[0];
    E(sumEdges[0]); E(sumEdges[threadsCount - 1]);Eo(disbalanceFactor);
    if (disbalanceFactor > 1.5 && false) {
        free(nextEdges);
        free(nextEdgesIds);
        goto clean;
    }
#pragma omp parallel
    {
        const int i = omp_get_thread_num();
        stickThisThreadToCore(i);
        const vid_t vertexBegin = int64_t(vertexCount) * i / threadsCount;
        const vid_t vertexEnd   = int64_t(vertexCount) * (i + 1) / threadsCount;
        for (vid_t v = vertexBegin; v < vertexEnd; ++v) {
            const vid_t nextv = que[v];
            for (eid_t e = edgesIds[nextv]; e < edgesIds[nextv + 1]; ++e) {
                nextEdges[nextEdgesIds[v] + e - edgesIds[nextv]] = edges[e];
                nextEdges[nextEdgesIds[v] + e - edgesIds[nextv]].dest = rev[edges[e].dest];
            }
        }
    }

    free(edges);
    free(edgesIds);
    edges = nextEdges;
    edgesIds = nextEdgesIds;
clean:
    //free(que);
    //free(rev);
    delete[] sumEdges;
}

void doReorderSimple() {
    stickThisThreadToCore(0);

    bool *visit = new bool[vertexCount]();
    que = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount));
    rev = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount));
    vid_t pos = 0;
    for (vid_t v = 0; v < vertexCount; ++v) if (!visit[v]) {
        visit[v] = true;
        que[pos++] = v;
        for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e) {
            const vid_t u = edges[e].dest;
            if (visit[u]) continue;
            visit[u] = true;
            que[pos++] = u;
        }
    }
    assert(pos == vertexCount);
    delete[] visit;

#pragma omp parallel for
    for (vid_t i = 0; i < vertexCount; ++i)
        rev[que[i]] = i;

    eid_t *nextEdgesIds = static_cast<eid_t*>(malloc(sizeof(eid_t) * (vertexCount + 1)));
    Edge *nextEdges = static_cast<Edge*>(malloc(sizeof(Edge) * edgesCount));
    nextEdgesIds[0] = 0;
    for (int i = 0; i < threadsCount; ++i) {
        stickThisThreadToCore(i);
        const vid_t vertexBegin = int64_t(vertexCount) * (i + 0) / threadsCount;
        const vid_t vertexEnd   = int64_t(vertexCount) * (i + 1) / threadsCount;
        for (vid_t v = vertexBegin; v < vertexEnd; ++v) {
            const vid_t nextv = que[v];
            nextEdgesIds[v + 1] = nextEdgesIds[v] + edgesIds[nextv + 1] - edgesIds[nextv];
        }
    }
#pragma omp parallel
    {
        const int i = omp_get_thread_num();
        stickThisThreadToCore(i);
        const vid_t vertexBegin = int64_t(vertexCount) * i / threadsCount;
        const vid_t vertexEnd   = int64_t(vertexCount) * (i + 1) / threadsCount;
        for (vid_t v = vertexBegin; v < vertexEnd; ++v) {
            const vid_t nextv = que[v];
            for (eid_t e = edgesIds[nextv]; e < edgesIds[nextv + 1]; ++e) {
                nextEdges[nextEdgesIds[v] + e - edgesIds[nextv]] = edges[e];
                nextEdges[nextEdgesIds[v] + e - edgesIds[nextv]].dest = rev[edges[e].dest];
            }
        }
    }
    free(que);
    //free(rev);
    free(edges);
    free(edgesIds);
    edges = nextEdges;
    edgesIds = nextEdgesIds;
}

//
// Read input data
//
void readAll(char *filename) {
    iterationNumber = 0;
#pragma omp parallel
    {
#pragma omp master
        {
            threadsCount = omp_get_num_threads();
        }
    }
    Eo(threadsCount);

    FILE *f = fopen(filename, "rb");
    assert(f);

    fread(&vertexCount, sizeof(vid_t), 1, f);
    fread(&edgesCount, sizeof(eid_t), 1, f);

    // NUMA-specific part
    // Each core allocate memory and read data, which will processed on this core
    edgesIds = (eid_t*)malloc(sizeof(eid_t) * (vertexCount + 1));
    edges = (Edge*)malloc(sizeof(Edge) * (edgesCount));

    for (int i = 0; i < threadsCount; ++i) {
        stickThisThreadToCore(i);
        const int vertexBegin = int64_t(vertexCount + 1) * i / threadsCount;
        const int vertexEnd   = int64_t(vertexCount + 1) * (i + 1) / threadsCount;
        fread(edgesIds + vertexBegin, sizeof(eid_t), vertexEnd - vertexBegin, f);
    }

    eid_t loopsCount = 0;
    for (int i = 0; i < threadsCount; ++i) {
        stickThisThreadToCore(i);
        const vid_t vertexBegin = int64_t(vertexCount) * i / threadsCount;
        const vid_t vertexEnd   = int64_t(vertexCount) * (i + 1) / threadsCount;
        const eid_t edgeBegin   = edgesIds[vertexBegin];
        const eid_t edgeEnd     = edgesIds[vertexEnd];
        fread(edges + edgeBegin, sizeof(Edge), edgeEnd - edgeBegin, f);
        for (vid_t v = vertexBegin; v < vertexEnd; ++v) {
            for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e)
                if (edges[e].dest == v) {
                    assert(vertexDegree(v) > 1);
                    edges[e].weight = MAX_WEIGHT;
                    ++loopsCount;
                }
        }
    }
    Eo(loopsCount);
    fclose(f);
}

void convertAll(graph_t *G) {
    iterationNumber = 0;
#pragma omp parallel
    {
#pragma omp master
        {
            threadsCount = omp_get_num_threads();
        }
    }

    //std::set<weight_t> allWeight;

    // TODO NUMA
    vertexCount = G->n;
    edgesCount = G->m;
    edgesIds = (eid_t*)malloc(sizeof(eid_t) * (vertexCount + 1));
    edges = (Edge*)malloc(sizeof(Edge) * (edgesCount));
#pragma omp parallel
    {
        stickThisThreadToCore(omp_get_thread_num());
#pragma omp for
        for (vid_t i = 0; i <= vertexCount; ++i)
            edgesIds[i] = static_cast<eid_t>(G->rowsIndices[i]);
#pragma omp for
        for (vid_t v = 0; v < vertexCount; ++v)
            for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e) {
                edges[e].dest = G->endV[e];
                edges[e].weight = G->weights[e];
                edges[e].origOffset = e - edgesIds[v];
                // remove loops
                if (edges[e].dest == v) {
                    assert(vertexDegree(v) > 1);
                    edges[e].weight = MAX_WEIGHT;
                }
            }
    }
    //Eo(allWeight.size());
}

//
// Timer functions
//
int64_t currentNanoTime() {
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return int64_t(ts.tv_sec) * int64_t(1e9) + ts.tv_nsec;
}

RDTSC::RDTSC() {
    double t1, t2;
    t1 = get();
    sleep(1);
    t2 = get();
    oneSecond = (t2 - t1);
}

double RDTSC::get() {
    unsigned int time_edx, time_eax;

#ifdef ON_JSCC
    unsigned a, d;
    asm volatile("rdtsc" : "=a" (a), "=d" (d));  
    return ((unsigned long long)a) | (((unsigned long long)d) << 32);
#endif

    asm volatile (  "rdtscp\n\t"
            "mov %%edx, %0\n\t"
            "mov %%eax, %1\n\t"
            //                    "cpuid\n\t"
            : "=r"(time_edx), "=r"(time_eax) ::
            "%rax", "%rbx", "%rcx", "%rdx");

    return (double)(((unsigned long long)time_edx) << 32 | (unsigned long long)time_eax);
}

void RDTSC::start(int timerId) {
    timers[timerId][0] = get();
}

double RDTSC::end(int timerId) {
    double res = (get() - timers[timerId][0]) / oneSecond;
    return res;
}

RDTSC rdtsc;

int stickThisThreadToCore(int coreId) {
    const int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
#ifdef ON_HOME
    coreId *= 2;
    coreId += 1;
#endif /* ON_HOME */
#ifdef USE_HYPERTHREADING
    int hardwareCoreId = coreId / 2;
    int htOffset = num_cores / 2;
    coreId = hardwareCoreId + htOffset * (coreId % 2);
#endif /* USE_HYPERTHREADING */
    if (coreId >= num_cores) return 0;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(coreId, &cpuset);

    pthread_t current_thread = pthread_self();
    return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}
