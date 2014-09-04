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

#include <iostream>
#define E(x) { cerr << #x << " = " << (x) << "   "; }
#define Eo(x) { cerr << #x << " = " << (x) << endl; }
#define EO(x) Eo(x)

typedef pair<weight_t, eid_t> pwe;

vid_t *comp;
weight_t taskResult;

int threadsCount;
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



bool doAll() {
    // find min
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
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
            const vid_t cv = comp[v];
            weight_t startWeight = edges[startEdgesIds[v]].weight;
            for (eid_t e = startEdgesIds[v]; e < edgesIds[v + 1]; ++e) {
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
                    startEdgesIds[v] = e;
                }
            }
        }
    }

    // reduce min
    int updated = 0;
#pragma omp parallel for reduction(+:updated)
    for (vid_t i = 0; i < vertexCount; ++i) {
        Result best{MAX_WEIGHT + 0.1, 0, 0, 0};
        if (comp[i] == i) {
            for (int j = 0; j < threadsCount; ++j) {
                best = min(best, localResult[j][i]); // !
                localResult[j][i].weight = MAX_WEIGHT + 0.1;
            }
        }
        bestResult[i] = best;
        if (best.weight <= MAX_WEIGHT) {
            comp[i] = best.destComp;
            updated = 1;
        }
    }
    if (!updated) return false;

    // merge components
    weight_t tmpTaskResult = 0.0;
#pragma omp parallel for reduction(+:tmpTaskResult)
    for (vid_t i = 0; i < vertexCount; ++i) {
        Result best = bestResult[i];
        if (best.weight > MAX_WEIGHT) continue;
        vid_t oc = best.destComp;

        if (comp[oc] == i) {
            if (i < oc) {
                comp[i] = i;
            } else {
                comp[i] = oc;
                tmpTaskResult += best.weight;
            }
        } else {
            tmpTaskResult += best.weight;
            comp[i] = oc;
        }
    }
    taskResult += tmpTaskResult;

    // pointer jumping
    int changed = 100500;
    while (changed) {
        changed = 0;
#pragma omp parallel for reduction(+:changed)
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
    }

    // TODO reduce edges?

    return true;
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


int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }
    readAll(argv[1]);

    int64_t prepareTime = -currentNanoTime();
    doPrepare();
    prepareTime += currentNanoTime();

    int64_t calcTime = -currentNanoTime();
    while (doAll());
    calcTime += currentNanoTime();

    printf("%.10lf\n", double(taskResult));

    fprintf(stderr, "%.3lf\n%.3lf\n", double(prepareTime) / 1e9, double(calcTime) / 1e9);

    return 0;
}
