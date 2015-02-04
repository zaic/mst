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
#ifdef ON_HOME
    coreId *= 1;
#endif
    const int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (coreId >= num_cores) return 0;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(coreId, &cpuset);

    pthread_t current_thread = pthread_self();
    return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}
