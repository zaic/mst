#include "gen.h"
#include <cassert>
#include <cstdio>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#ifdef __clang__
#include "omp.h"
#else
#include <omp.h>
#endif

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
    return a.weight < b.weight;
}

void readAll(char *filename) {
    FILE *f = fopen(filename, "rb");
    assert(f);

    fread(&vertexCount, sizeof(vid_t), 1, f);
    fread(&edgesCount, sizeof(eid_t), 1, f);

    // NUMA-specific part
    // Each core allocate memory and read data, which will processed on this core
    int curThreadsCount;
#pragma omp parallel
    {
#pragma omp master
        {
            curThreadsCount = omp_get_num_threads();
        }
    }
    Eo(curThreadsCount);

    edgesIds = (eid_t*)malloc(sizeof(eid_t) * (vertexCount + 1));
    edges = (Edge*)malloc(sizeof(Edge) * (edgesCount));

    for (int i = 0; i < curThreadsCount; ++i) {
        stickThisThreadToCore(i);
        const int vertexBegin = int64_t(vertexCount + 1) * i / curThreadsCount;
        const int vertexEnd   = int64_t(vertexCount + 1) * (i + 1) / curThreadsCount;
        fread(edgesIds + vertexBegin, sizeof(eid_t), vertexEnd - vertexBegin, f);
    }

    for (int i = 0; i < curThreadsCount; ++i) {
        stickThisThreadToCore(i);
        const int vertexBegin = int64_t(vertexCount) * i / curThreadsCount;
        const int vertexEnd   = int64_t(vertexCount) * (i + 1) / curThreadsCount;
        const int edgeBegin   = edgesIds[vertexBegin];
        const int edgeEnd     = edgesIds[vertexEnd];
        fread(edges + edgeBegin, sizeof(Edge), edgeEnd - edgeBegin, f);
    }

#ifdef USE_REORDER
    stickThisThreadToCore(0);

    bool *visit = new bool[vertexCount];
    vid_t *que = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount));
    vid_t *rev = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount));
    vid_t fr = 0, bc = 0;
    for (vid_t i = 0; i < vertexCount; ++i) if (!visit[i]) {
        que[bc++] = i;
        while (fr < bc) {
            const vid_t v = que[fr++];
            visit[v] = true;
            for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e) {
                const vid_t u = edges[e].dest;
                if (visit[u]) continue;
                visit[u] = true;
                que[bc++] = u;
            }
        }
    }
    assert(bc == vertexCount);
    delete[] visit;

#pragma omp for
    for (vid_t i = 0; i < vertexCount; ++i)
        rev[que[i]] = i;

    eid_t *nextEdgesIds = static_cast<eid_t*>(malloc(sizeof(eid_t) * (vertexCount + 1)));
    Edge *nextEdges = static_cast<Edge*>(malloc(sizeof(Edge) * edgesCount));
    nextEdgesIds[0] = 0;
    for (int i = 0; i < curThreadsCount; ++i) {
        stickThisThreadToCore(i);
        const vid_t vertexBegin = int64_t(vertexCount + 1) * i / curThreadsCount;
        const vid_t vertexEnd   = int64_t(vertexCount + 1) * (i + 1) / curThreadsCount;
        for (vid_t v = vertexBegin; v < vertexEnd; ++v) {
            const vid_t nextv = que[v];
            nextEdgesIds[v + 1] = nextEdgesIds[v] + edgesIds[nextv + 1] - edgesIds[nextv];
        }
    }
#pragma omp parallel
    {
        const int i = omp_get_thread_num();
        stickThisThreadToCore(i);
        const vid_t vertexBegin = int64_t(vertexCount + 1) * i / curThreadsCount;
        const vid_t vertexEnd   = int64_t(vertexCount + 1) * (i + 1) / curThreadsCount;
        for (vid_t v = vertexBegin; v < vertexEnd; ++v) {
            const vid_t nextv = que[v];
            for (eid_t e = edgesIds[nextv]; e < edgesIds[nextv + 1]; ++e) {
                nextEdges[nextEdgesIds[v] + e - edgesIds[nextv]] = edges[e];
                nextEdges[nextEdgesIds[v] + e - edgesIds[nextv]].dest = rev[edges[e].dest];
            }
        }
    }
    free(que);
    free(edges);
    free(edgesIds);
    edges = nextEdges;
    edgesIds = nextEdgesIds;
#endif /* USE_REORDER */

    fclose(f);
}

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
