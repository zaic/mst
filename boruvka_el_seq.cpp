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
    {
        const int threadId = 0;//omp_get_thread_num();
        auto result = localResult[threadId];

#if 0
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            vid_t icomp = comp[i];
            result[icomp] = make_pair(MAX_WEIGHT + 0.1, 0);
        }
#else
        for (vid_t i = 0; i < vertexCount; ++i) {
            result[i].weight = MAX_WEIGHT + 0.1;
        }
#endif
        for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) {
            const vid_t cv = comp[v];
            for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e) {
                if (edges[e].weight > result[cv].weight) continue;
                const vid_t u = edges[e].dest;
                const vid_t cu = comp[u];
                if (cu == cv) continue;
                if (edges[e].weight < result[cv].weight || (edges[e].weight == result[cv].weight && cu < result[cv].destComp)) {
                    result[cv] = Result{edges[e].weight, cu, v, u};
                }
            }
        }
    }

    // reduce min
    int updated = 0;
    for (vid_t i = 0; i < vertexCount; ++i) {
        Result best{MAX_WEIGHT + 0.1, 0, 0, 0};
        for (int j = 0; j < threadsCount; ++j) {
            best = min(best, localResult[j][i]); // !
        }
        bestResult[i] = best;
        if (best.weight <= MAX_WEIGHT) {
            comp[i] = best.destComp;
            updated = 1;
        }
    }
    if (!updated) return false;

    // merge components
    for (vid_t i = 0; i < vertexCount; ++i) {
        Result best = bestResult[i];
        if (best.weight > MAX_WEIGHT) continue;
        vid_t oc = best.destComp;

        if (comp[oc] == i) {
            if (i < oc) {
                comp[i] = i;
            } else {
                comp[i] = oc;
            }
        } else {
            taskResult += best.weight;
            comp[i] = oc;
        }
    }

    // pointer jumping
    int changed = 100500;
    while (changed) {
        changed = 0;
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

//#pragma omp parallel
    {
        const int threadId = 0;//omp_get_thread_num();
        const int curThreadsCount = 1;//omp_get_num_threads();
        if (!threadId) {
            threadsCount = curThreadsCount;
            localResult = new Result*[threadsCount];
        }
//#pragma omp barrier

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

        /*
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i)
            sort(edges + edgesIds[i], edges + edgesIds[i + 1], EdgeDestCmp());
            */
    }
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }
    readAll(argv[1]);

    doPrepare();
    while (doAll());

    printf("%.10lf\n", double(taskResult));

    return 0;
}
