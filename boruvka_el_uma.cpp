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

vid_t *vertexIds;
eid_t *startEdgesIds;

struct Result {
    weight_t weight;
    vid_t destComp, from, to;

    bool operator<(const Result& o) const {
        if (weight != o.weight) return weight < o.weight;
        return destComp < o.destComp;
    }
};

Result **localResult, *bestResult;

double times[kMaxThreads][kMaxIterations][40];



bool doAll() {
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
        auto result = localResult[threadId];

#if 0
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            vid_t icomp = comp[i];
            result[icomp] = make_pair(MAX_WEIGHT + 0.1, 0);
        }
#elif 0
        for (vid_t i = 0; i < vertexCount; ++i) {
            result[i].weight = MAX_WEIGHT + 0.1;
        }
#endif
        for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) {
            const eid_t edgesStart = startEdgesIds[v];
            const eid_t edgesEnd = edgesIds[v + 1];
            if (edgesStart >= edgesEnd) continue;
            
            const vid_t cv = comp[v];
            weight_t startWeight = edges[edgesStart].weight;
            eid_t newEdgesStart = edgesStart;

            for (eid_t e = edgesStart; e < edgesEnd; ++e) {
                weight_t weight = edges[e].weight;
                if (weight > result[cv].weight) break;

                const vid_t u = edges[e].dest;
                const vid_t cu = comp[u];
                if (cu == cv) continue;

                if (weight < result[cv].weight || (weight == result[cv].weight && cu < result[cv].destComp)) {
                    result[cv] = Result{edges[e].weight, cu, v, u};
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
#pragma omp for reduction(+:updated) nowait
        for (vid_t i = 0; i < vertexCount; ++i) {
            if (comp[i] == i) {
                Result best{MAX_WEIGHT + 0.1, 0, 0, 0};
                for (int j = 0; j < threadsCount; ++j) {
                    best = min(best, localResult[j][i]); // !
                    localResult[j][i].weight = MAX_WEIGHT + 0.1;
                }
                bestResult[i] = best;
                if (best.weight <= MAX_WEIGHT) {
                    comp[i] = best.destComp;
                    updated = 1;
                }
            } else {
                //bestResult[i].weight = MAX_WEIGHT + 0.1;
            }
        }
#else
        for (int treeIteration = 0; (1 << treeIteration) < threadsCount; treeIteration++) {
            const int reduceTo = ((threadId >> (treeIteration + 1)) << (treeIteration + 1));
            const int reduceFrom = reduceTo + (1 << treeIteration);
            const int threadsPerComp = (1 << (treeIteration + 1));
            const vid_t reduceStartComp = int64_t(vertexCount) * (threadId & ((1 << (treeIteration + 1)) - 1)) / threadsPerComp;
            vid_t reduceEndCompPre = int64_t(vertexCount) * ((threadId + 1) & ((1 << (treeIteration + 1)) - 1)) / threadsPerComp;
            const vid_t reduceEndComp = (reduceEndCompPre == 0 ? vertexCount : reduceEndCompPre);

#if 0
#pragma omp master
            {
                E("==================");
                Eo(treeIteration);
            }
#pragma omp barrier
#pragma omp critical
            {
                E(threadId); 
                E(reduceStartComp); E(reduceEndComp);
                E(reduceTo); Eo(reduceFrom);
            }
#endif

            for (vid_t i = reduceStartComp; i < reduceEndComp; ++i) {
                if (comp[i] == i) {
                    if (localResult[reduceFrom][i] < localResult[reduceTo][i])
                        localResult[reduceTo][i] = localResult[reduceFrom][i];
                    localResult[reduceFrom][i].weight = MAX_WEIGHT + 0.1;
                }
            }
#pragma omp barrier
        }

#pragma omp for reduction(+:updated) nowait
        for (vid_t i = 0; i < vertexCount; ++i) {
            if (comp[i] == i) {
                bestResult[i] = localResult[0][i];
                if (localResult[0][i].weight <= MAX_WEIGHT) {
                    comp[i] = localResult[0][i].destComp;
                    updated = 1;
                }
                localResult[0][i].weight = MAX_WEIGHT + 0.1;
            }
        }
#endif
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

    // TODO reduce edges?

    taskResult += tmpTaskResult;
    ++iterationNumber;
    return updated;
}

void doPrepare() {
    vertexIds = new vid_t[vertexCount + 1];
    vertexIds[0] = 0;

    comp = new vid_t[vertexCount];
    for (vid_t i = 0; i < vertexCount; ++i)
        comp[i] = i;

    bestResult = new Result[vertexCount];
    startEdgesIds = new eid_t[vertexCount];

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        const int curThreadsCount = omp_get_num_threads();
        if (!threadId) {
            threadsCount = curThreadsCount;
            localResult = new Result*[threadsCount];
        }
#pragma omp barrier

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

        localResult[threadId] = new Result[vertexCount];
        for (vid_t i = 0; i < vertexCount; ++i)
            localResult[threadId][i] = Result{MAX_WEIGHT + 0.1, 0, 0, 0};

        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
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

#if 0
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
