#include "gen.h"
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <set>
#include <tuple>
#include <vector>
#include <numeric>
#ifdef __clang__
#include "omp.h"
#else
#include <omp.h>
#endif
using namespace std;

typedef pair<weight_t, eid_t> pwe;

vid_t *comp;
weight_t taskResult;

int threadsCount;
#ifndef ON_NUMA
#error NUMA only
#else
vid_t vertexesPerThread;
#endif
eid_t *startEdgesIds;
int iterationNumber;

struct Result {
    weight_t weight;
    vid_t destComp, from, to;

    bool operator<(const Result& o) const {
        if (weight != o.weight) return weight < o.weight;
        return destComp < o.destComp;
    }
};

Result *bestResult;

const int kMaxIterations = 30;
const int kMaxThreads = 20;
double times[kMaxIterations][kMaxThreads][40];



bool doAll() {
    Eo(iterationNumber);
    int updated = 0; // reduce stage 
    weight_t tmpTaskResult = 0.0; // merge stage

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);

        //
        // find min
        // 
        rdtsc.start(threadId);
        //auto result = localResult[threadId];

        const vid_t startVertex = int64_t(vertexCount) * threadId / threadsCount;
        for (vid_t vo = 0; vo < vertexCount; ++vo) { // iterate over all vertexes with specified offset
            const vid_t v = vo; // TODO (startVertex + vo) % vertexCount;
            const vid_t cv = comp[v];
            const int cvThread = cv / vertexesPerThread;
            if (cvThread != threadId) continue; // skip vertexes, which belong to another thread

            const eid_t edgesStart = startEdgesIds[v];
            const eid_t edgesEnd = edgesIds[v + 1];
            if (edgesStart >= edgesEnd) continue;
            
            weight_t startWeight = edges[edgesStart].weight;
            eid_t newEdgesStart = edgesStart;

            for (eid_t e = edgesStart; e < edgesEnd; ++e) {
                weight_t weight = edges[e].weight;
                if (weight > bestResult[cv].weight) break;

                const vid_t u = edges[e].dest;
                const vid_t cu = comp[u];
                if (cu == cv) continue;

                if (weight < bestResult[cv].weight || (weight == bestResult[cv].weight && cu < bestResult[cv].destComp)) {
                    bestResult[cv] = Result{edges[e].weight, cu, v, u};
                }
                if (weight > startWeight) {
                    startWeight = weight;
                    newEdgesStart = e;
                }
            }

            if (newEdgesStart != edgesStart) startEdgesIds[v] = newEdgesStart;
        }
        times[iterationNumber][threadId][0] = rdtsc.end(threadId);

#pragma omp barrier

        //
        // reduce min
        //
        rdtsc.start(threadId);
#if 0
        for (int treeIteration = 0; (1 << treeIteration) < threadsCount; treeIteration++) {
            const int reduceTo = ((threadId >> (treeIteration + 1)) << (treeIteration + 1));
            const int reduceFrom = reduceTo + (1 << treeIteration);
            const int threadsPerComp = (1 << (treeIteration + 1));
            const vid_t reduceStartComp = int64_t(vertexCount) * (threadId & ((1 << (treeIteration + 1)) - 1)) / threadsPerComp;
            vid_t reduceEndCompPre = int64_t(vertexCount) * ((threadId + 1) & ((1 << (treeIteration + 1)) - 1)) / threadsPerComp;
            const vid_t reduceEndComp = (reduceEndCompPre == 0 ? vertexCount : reduceEndCompPre);

            for (vid_t i = reduceStartComp; i < reduceEndComp; ++i) {
                if (comp[i] == i) {
                    if (localResult[reduceFrom][i] < localResult[reduceTo][i])
                        localResult[reduceTo][i] = localResult[reduceFrom][i];
                    localResult[reduceFrom][i].weight = MAX_WEIGHT + 0.1;
                }
            }
#pragma omp barrier
        }
#endif

#pragma omp for reduction(+:updated) nowait
        for (vid_t i = 0; i < vertexCount; ++i) {
            if (comp[i] == i) {
                //bestResult[i] = localResult[0][i];
                if (bestResult[i].weight <= MAX_WEIGHT) {
                    comp[i] = bestResult[i].destComp;
                    updated = 1;
                }
                //bestResult[i].weight = MAX_WEIGHT + 0.1;
            }
        }

        times[iterationNumber][threadId][1] = rdtsc.end(threadId);
        // if (!updated) return false; TODO

#pragma omp barrier

        //
        // merge components
        //
        rdtsc.start(threadId);
#pragma omp for reduction(+:tmpTaskResult) nowait
        for (vid_t i = 0; i < vertexCount; ++i) {
            Result best = bestResult[i];
            if (best.weight > MAX_WEIGHT) continue;
            vid_t oc = best.destComp;

            if (comp[oc] == i) {
                if (i < oc) {
                    comp[i] = i;
                } else {
                    comp[i] = oc;
                    bestResult[i].weight = MAX_WEIGHT + 0.1;
                    tmpTaskResult += best.weight;
                }
                bestResult[i].weight = MAX_WEIGHT + 0.1;
            } else {
                tmpTaskResult += best.weight;
                comp[i] = oc;
                bestResult[i].weight = MAX_WEIGHT + 0.1;
            }
        }
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);
    }

    // pointer jumping
    int changed = 100500;
    while (changed) {
        changed = 0;
#pragma omp parallel 
        {
            const int threadId = omp_get_thread_num();
            rdtsc.start(threadId);
#pragma omp for reduction(+:changed) nowait
            for (vid_t i = 0; i < vertexCount; ++i) {
                vid_t myComp = comp[i];
                if (myComp == i) continue;
                vid_t parentComp = comp[myComp];
                if (parentComp == i) {
                    comp[i] = min(i, myComp);
                    changed = 1;
                } else if (myComp != parentComp) {
                    comp[i] = parentComp;
                    changed = 1;
                }
            }
            times[iterationNumber][threadId][3] = rdtsc.end(threadId);
        }
    }

    taskResult += tmpTaskResult;
    ++iterationNumber;
    return updated;
}

void doPrepare() {
#ifndef ON_NUMA
    vertexIds = new vid_t[vertexCount + 1];
    vertexIds[0] = 0;
#endif

    comp = new vid_t[vertexCount];
    // TODO omp for
    for (vid_t i = 0; i < vertexCount; ++i)
        comp[i] = i;

    bestResult = new Result[vertexCount];
    startEdgesIds = new eid_t[vertexCount];

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        const int curThreadsCount = omp_get_num_threads();
        // TODO stick thread
        if (!threadId) {
            threadsCount = curThreadsCount;
            vertexesPerThread = vertexCount / threadsCount;
            assert(vertexesPerThread * threadsCount == vertexCount);
        }
#pragma omp barrier

#ifndef ON_NUMA
        eid_t degreeEnd = int64_t(edgesCount) * (threadId + 1) / threadsCount;
        eid_t degreeSum = 0;
        for (vid_t i = 0; i < vertexCount; ++i) {
            eid_t diff = edgesIds[i + 1] - edgesIds[i];
            degreeSum += diff;
            if (degreeSum >= degreeEnd) {
                vertexIds[threadId + 1] = i + 1;
                break;
            }
        }
#else
        //vertexIds[threadId + 1] = int64_t(vertexCount) * (threadId + 1) / threadsCount;
#endif

#pragma omp for nowait
        for (vid_t i = 0; i < vertexCount; ++i)
            bestResult[i] = Result{MAX_WEIGHT + 0.1, 0, 0, 0};

#pragma omp for nowait
        for (vid_t i = 0; i < vertexCount; ++i) {
            sort(edges + edgesIds[i], edges + edgesIds[i + 1], EdgeWeightCmp());
            startEdgesIds[i] = edgesIds[i];
        }

        //E(threadId); Eo(vertexIds[threadId + 1]);
    }
}

void warmup() {
    const int64_t iterationCount = 1e9;
#pragma omp parallel
    {
        volatile int64_t fpre = 0, fcur = 1;
        stickThisThreadToCore(omp_get_thread_num());
        for (int64_t i = iterationCount; i > 0; --i) {
            int64_t fnext = fpre + fcur;
            fpre = fcur;
            fcur = fnext;
        }
    }
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }
    readAll(argv[1]);
    warmup();
    fprintf(stderr, "Done\n");

    int64_t prepareTime = -currentNanoTime();
    doPrepare();
    prepareTime += currentNanoTime();

    int64_t calcTime = -currentNanoTime();
    while (doAll());
    calcTime += currentNanoTime();

    printf("%.10lf\n", double(taskResult));

    fprintf(stderr, "%.3lf\n%.3lf\n", double(prepareTime) / 1e9, double(calcTime) / 1e9);

#if 1
    for (int i = 0; i < iterationNumber; ++i) {
        fprintf(stderr, "iteration %2d\n", i);
        for (int j = 0; j < threadsCount; ++j) {
            fprintf(stderr, "%02d:   ", j);
            for (int k = 0; k < 4; ++k)
                fprintf(stderr, "%.6lf ", times[i][j][k]);
            fputs("\n", stderr);
        }
    }
#endif

    return 0;
}
