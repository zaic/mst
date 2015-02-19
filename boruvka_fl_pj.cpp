#include "gen.h"
#include <cstring>
#include <cassert>
#include <cstdio>
#include <queue>
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
#include "vector.h"
using namespace std;

vid_t *comp;
weight_t taskResult;

#ifdef ON_NUMA
vid_t vertexesPerThread;
#else
#define vertexesPerThread #error numa_only
#endif
vid_t *vertexIds;

double times[kMaxThreads][kMaxIterations][40];


// temporary min-edge data
eid_t *bestEid; // TODO very temporary array
vid_t *bestComp; // component, which is our pointer is pointing to

/*
 *  FL data
 */
namespace FlData {
    struct ListElement {
        vid_t vertex;       // 4 bytes   0... 7
        eid_t startEdgeId;  // 8 bytes   8...15
        eid_t endEdgeId;    // 8 bytes  16...23
    };
    typedef Vector<ListElement, 10> SmallVector;
#ifdef USE_SMALL_VECTOR
    SmallVector *list;
#else
    std::vector<ListElement> *list;
#endif

#ifdef USE_BOUND
    const int kVertexBound = 100000;
    const int kBoundQueue = 100;

    std::vector<vid_t> parallelProcess[kMaxThreads];
    eid_t bestEidPerThread[kMaxThreads][kMaxThreads][kVertexBound]; // [owner thread][processing thread][vertex id]
#endif /* USE_BOUND */



    void allocFalData() {
#ifdef USE_SMALL_VECTOR
        list = new SmallVector[vertexCount];
#else
        list = new std::vector<ListElement>[vertexCount];
#endif
    }
}




bool doAll() {
    int updated = 0; // reduce stage 
    weight_t tmpTaskResult = 0.0; // merge stage
    Eo(iterationNumber);

#pragma omp parallel reduction(+:updated) reduction(+:tmpTaskResult)
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);

        //
        // find minimum edge
        //
        rdtsc.start(threadId);
#ifdef USE_SKIP_LOOPS
        if (iterationNumber > USE_SKIP_LOOPS) {
            for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) {
                const eid_t edgesStart = startEdgesIds[v];
                const eid_t edgesEnd = edgesIds[v + 1];
                if (edgesStart >= edgesEnd) continue;

                const vid_t cv = comp[v];
                eid_t newEdgesStart = edgesStart;

                for (; newEdgesStart < edgesEnd; ++newEdgesStart) {
                    const vid_t u = edges[newEdgesStart].dest;
                    const vid_t cu = comp[u];
                    if (cu == cv) continue;
                    break;
                }
                if (newEdgesStart != edgesStart) startEdgesIds[v] = newEdgesStart;
            }
        }
        times[iterationNumber][threadId][0] = rdtsc.end(threadId);
#pragma omp barrier
#endif

        //
        // find minimum edge
        //
        rdtsc.start(threadId);

#ifdef USE_BOUND
        FlData::parallelProcess[threadId].clear();
#endif /* USE_BOUND */
        const vid_t vto = vertexIds[threadId + 1];
        for (vid_t v = vertexIds[threadId]; v < vto; ++v) if (comp[v] == v) {
            if (FlData::list[v].size() >= FlData::kVertexBound) {
                FlData::parallelProcess[threadId].push_back(v);
                continue;
            }
            weight_t curBestWeight = MAX_WEIGHT + 1e-3;
            eid_t curBestEid = -1;
            
#ifdef USE_SMALL_VECTOR
            auto& data = FlData::list[v];
#else
            vector<FlData::ListElement>& data = FlData::list[v];
#endif /* USE_SMALL_VECTOR */
            for (size_t i = 0; i < data.size(); ) {
                const vid_t curVertex = data[i].vertex;
                eid_t curEdgeId = data[i].startEdgeId;
                const eid_t endEdgeId = data[i].endEdgeId;
                while (curEdgeId < endEdgeId && comp[edges[curEdgeId].dest] == v) ++curEdgeId;
                if (curEdgeId == endEdgeId) {
                    // remove
#ifdef USE_SMALL_VECTOR
                    data.removeAt(i);
#else
                    data[i] = data.back();
                    data.pop_back();
#endif /* USE_SMALL_VECTOR */
                    continue;
                } else {
                    data[i].startEdgeId = curEdgeId; // ToDo omptimize: write only when changed
                }
                if (edges[curEdgeId].weight < curBestWeight) {
                    curBestWeight = edges[curEdgeId].weight;
                    curBestEid = curEdgeId;
                }
                ++i;
            }
            bestEid[v] = curBestEid;
            bestComp[v] = (curBestEid == -1 ? -1 : comp[edges[curBestEid].dest]);
        } else {
            bestEid[v] = -1;
        }
        assert(FlData::parallelProcess[threadId].size() < FlData::kBoundQueue);

        times[iterationNumber][threadId][1] = rdtsc.end(threadId);
#pragma omp barrier
        rdtsc.start(threadId);

#ifdef USE_BOUND
        for (int t = 0; t < threadsCount; ++t) for (size_t vid = 0; vid < FlData::parallelProcess[t].size(); ++vid) {
            const vid_t v = FlData::parallelProcess[t][vid];
            weight_t curBestWeight = MAX_WEIGHT + 1e-3;
            eid_t curBestEid = -1;
#ifdef USE_SMALL_VECTOR
            auto& data = FlData::list[v];
#else
            vector<FlData::ListElement>& data = FlData::list[v];
#endif /* USE_SMALL_VECTOR */

            size_t ifrom = data.size() * (threadId + 0) / threadsCount;
            size_t ito   = data.size() * (threadId + 1) / threadsCount;
            for (size_t i = ifrom; i < ito; ++i) {
                const vid_t curVertex = data[i].vertex;
                eid_t curEdgeId = data[i].startEdgeId;
                const eid_t endEdgeId = data[i].endEdgeId;
                while (curEdgeId < endEdgeId && comp[edges[curEdgeId].dest] == v) ++curEdgeId;
                if (curEdgeId == endEdgeId) {
                    // skip (instead of remove)
                    // TODO: set startEdgeId
                    continue;
                } else {
                    data[i].startEdgeId = curEdgeId; // ToDo omptimize: write only when changed. or remove
                }
                if (edges[curEdgeId].weight < curBestWeight) {
                    curBestWeight = edges[curEdgeId].weight;
                    curBestEid = curEdgeId;
                }
            }
            FlData::bestEidPerThread[t][threadId][vid] = curBestEid;
        }
#endif /* USE_BOUND */
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);
#pragma omp barrier
        rdtsc.start(threadId);

        for (size_t vid = 0; vid < FlData::parallelProcess[threadId].size(); ++vid) {
            const vid_t v= FlData::parallelProcess[threadId][vid];
            weight_t curBestWeight = MAX_WEIGHT + 1e-3;
            eid_t curBestEid = -1;
            for (int t = 0; t < threadsCount; ++t) {
                eid_t teid = FlData::bestEidPerThread[threadId][t][vid];
                if (teid != -1 && edges[teid].weight < curBestWeight) {
                    curBestWeight = edges[teid].weight;
                    curBestEid = teid;
                }
            }
            bestEid[v] = curBestEid;
            bestComp[v] = (curBestEid == -1 ? -1 : comp[edges[curBestEid].dest]);
        }

        times[iterationNumber][threadId][2] += rdtsc.end(threadId);

#pragma omp barrier

        //
        // recalc answer, add eges to MST ("merge")
        //
        rdtsc.start(threadId);
        
#pragma omp for nowait
        for (vid_t v = 0; v < vertexCount; ++v) if (bestEid[v] != -1) {
            comp[v] = bestComp[v];
        }

        times[iterationNumber][threadId][3] += rdtsc.end(threadId);
#pragma omp barrier
        rdtsc.start(threadId);

#pragma omp for nowait
        for (vid_t i = 0; i < vertexCount; ++i) {
            const eid_t edge = bestEid[i];
            if (edge == -1) continue;
            updated = 1;
            //Result best = bestResult[i];
            //if (best.weight > MAX_WEIGHT) continue;
            vid_t oc = bestComp[i]; // best.destComp;

            if (comp[oc] == i) {
                //if (i < oc) {
                bool crazeCond = ((i^oc)&1);
                if ((crazeCond && i < oc) || (!crazeCond && i > oc)) {
                    comp[i] = i;
                } else {
                    comp[i] = oc;
                    tmpTaskResult += edges[edge].weight;
                }
            } else {
                tmpTaskResult += edges[edge].weight;
                comp[i] = oc;
            }
            bestEid[i] = -1; // TODO remove ?
        }

        times[iterationNumber][threadId][3] += rdtsc.end(threadId);
    }

    //
    // pointer jumping
    // 
    int changed = 100500;
    while (changed) {
        changed = 0;
#pragma omp parallel 
        {
            const int threadId = omp_get_thread_num();
            stickThisThreadToCore(threadId);
            rdtsc.start(threadId);
#pragma omp for reduction(+:changed) nowait
            for (vid_t i = 0; i < vertexCount; ++i) {
                vid_t myComp = comp[i];
                if (myComp == i) continue;
                vid_t parentComp = comp[myComp];
                if (parentComp == i) {
#if 0
                    comp[i] = min(i, myComp);
#else
                    if (min(i, myComp)&1)
                        comp[i] = min(i, myComp);
                    else
                        comp[i] = max(i, myComp);
#endif
                    changed = 1;
                } else if (myComp != parentComp) {
                    comp[i] = parentComp;
                    changed = 1;
                }
            }
            times[iterationNumber][threadId][4] += rdtsc.end(threadId);
        }
    }
    
    //
    // merge vertexes lists
    //
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);
        rdtsc.start(threadId);

        const vid_t myIdsFrom = vertexIds[threadId];
        const vid_t myIdsToto = vertexIds[threadId + 1];
        // TODO create bitvector readyForCopy ?
        if (iterationNumber == 0) {
            for (vid_t u = 0; u < vertexCount; ++u) {
                const vid_t v = comp[u];
                if (v == u || v < myIdsFrom || v >= myIdsToto) continue;
#ifdef USE_SMALL_VECTOR
                FlData::list[v].push_back(FlData::list[u].onStack[0]);
#else
                FlData::list[v].push_back(FlData::list[u][0]);
                FlData::list[u].clear();
#endif /* USE_SMALL_VECTOR */
            }
        } else {
            for (vid_t u = 0; u < vertexCount; ++u) { // iterate over all vertexes
                const vid_t v = comp[u];
                if (v == u) continue;
                //if (v / vertexesPerThread != threadId) continue; // check if this vertex points to our thread
                if (v < myIdsFrom || v >= myIdsToto) continue;
                if (FlData::list[u].empty()) continue;

                // merge vectors
#ifdef USE_SMALL_VECTOR
                // TODO implement swap for vector
                FlData::list[v].append(FlData::list[u]);
#else
                if (FlData::list[v].size() > FlData::list[u].size()) {
                    FlData::list[v].insert(FlData::list[v].end(), FlData::list[u].begin(), FlData::list[u].end());
                } else {
                    FlData::list[u].insert(FlData::list[u].end(), FlData::list[v].begin(), FlData::list[v].end());
                    FlData::list[u].swap(FlData::list[v]);
                }

                FlData::list[u].clear();
#endif /* USE_SMALL_VECTOR */
            }
        }

        times[iterationNumber][threadId][5] = rdtsc.end(threadId);
    }

    taskResult += tmpTaskResult;
    ++iterationNumber;
    return updated;
}

void doPrepare() {
    comp = new vid_t[vertexCount];
    for (vid_t i = 0; i < vertexCount; ++i)
        comp[i] = i;

    bestEid = new eid_t[vertexCount];
    bestComp = new vid_t[vertexCount];

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        const int curThreadsCount = omp_get_num_threads();
        stickThisThreadToCore(threadId);
#pragma omp master
        {
            threadsCount = curThreadsCount;
            assert(threadsCount < kMaxThreads);

            vertexIds = new vid_t[threadsCount + 1];
            vertexIds[0] = 0;

            FlData::allocFalData();
#ifdef ON_NUMA
            vertexesPerThread = vertexCount / threadsCount;
            assert(vertexesPerThread * threadsCount == vertexCount);
#endif /* ON_NUMA */
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
        vertexIds[threadId + 1] = int64_t(vertexCount) * (threadId + 1) / threadsCount;
#endif /* ON_NUMA */
#pragma omp barrier

        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
#if 1
            FlData::list[i].reserve(10);
#endif
            FlData::list[i].push_back(FlData::ListElement{i, edgesIds[i], edgesIds[i + 1]});

            std::sort(edges + edgesIds[i], edges + edgesIds[i + 1], EdgeWeightCmp());
        }
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
    //warmup();
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
            for (int k = 0; k < 6; ++k) // skip-loops, find-min-local, find-min-large-parallel, add edges to mst, pj, merge
                fprintf(stderr, "%.6lf ", times[i][j][k]);
            fputs("\n", stderr);
        }
    }
#endif

    return 0;
}
