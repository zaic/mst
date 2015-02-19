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

typedef pair<weight_t, eid_t> pwe;

vid_t *comp;
weight_t taskResult;

#ifdef ON_NUMA
vid_t vertexesPerThread;
#else
#define vertexesPerThread #error numa_only
#endif
vid_t *vertexIds;

double times[kMaxThreads][kMaxIterations][40];


/*
 *  temporary bfs data
 */
#ifdef USE_SMALL_VECTOR
Vector<vid_t, 100> *gComp;
#else
vector<vid_t> *gComp; // adjacency list
#endif
vector<pvv> graph_messages[kMaxThreads][kMaxThreads];

vector<vid_t> *fullComp; // contains entire new component
eid_t *bestEid; // TODO very temporary array
const int kMaxBfsQueue = 10000;
vid_t bfs_queue[kMaxThreads][kMaxBfsQueue];
signed char *bfs_visited[kMaxThreads];

/*
 *  FL data
 */
namespace FlData {
    pvv *list;
    vid_t *listBegin, *listEnd;
    eid_t *startEdge;

    void allocFalData() {
        list      = (pvv*)  malloc(vertexCount * sizeof(pvv));  
        listBegin = (vid_t*)malloc(vertexCount * sizeof(vid_t)); // vertex -> list
        listEnd   = (vid_t*)malloc(vertexCount * sizeof(vid_t));
        startEdge = (eid_t*)malloc(vertexCount * sizeof(eid_t));
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

#if 0
#ifdef USE_EDGE_STRUCT
        auto curEdges = edgesByThread[threadId];
#endif
        auto curEdgesIds = edgesIdsByThread[threadId];
        for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) {
            // E(v); E(curEdgesIds[v]); Eo(curEdgesIds[v + 1]);
            if (curEdgesIds[v] >= curEdgesIds[v + 1]) { // there is no edges from vertex v on this thread
                bestEid[v] = -1;
                continue;
            }
            assert(comp[v] == v);
            eid_t bestId = curEdgesIds[v];
            for (eid_t e = curEdgesIds[v] + 1; e < curEdgesIds[v + 1]; ++e) {
#ifdef USE_EDGE_STRUCT
                if (curEdges[e] < curEdges[bestId])
#else
                if (ExtEdge::compareEdge(threadId, e, threadId, bestId))
#endif
                    bestId = e;
            }
#ifdef USE_EDGE_STRUCT
            bestEid[v] = curEdges[bestId].origId;
#else
            bestEid[v] = ExtEdge::origId[threadId][bestId];
#endif
        }
#endif

        const vid_t vto = vertexIds[threadId + 1];
        for (vid_t v = vertexIds[threadId]; v < vto; ++v) if (comp[v] == v) {
            weight_t curBestWeight = MAX_WEIGHT;
            eid_t curBestEid = -1;
            
            vid_t prevComp = -1;
            for (vid_t curComp = FlData::listBegin[v]; curComp != -1; ) {
                vid_t curVertex = FlData::list[curComp].first;
                eid_t curEdgeId = FlData::startEdge[curVertex];
                while (curEdgeId < edgesIds[curVertex + 1] && comp[edges[curEdgeId].dest] == v) ++curEdgeId;
                if (curEdgeId == edgesIds[curVertex + 1]) {
                    // remove
                    //E(v); E(curComp); Eo("remove");
                    const vid_t nextComp = FlData::list[curComp].second;
                    //prevComp = curComp;
                    //curComp = nextComp;
                    if (prevComp == -1) {
                        FlData::listBegin[v] = nextComp;
                    } else {
                        FlData::list[prevComp].second = nextComp;
                    }
                    curComp = nextComp;
                    continue;
                } else {
                    FlData::startEdge[curVertex] = curEdgeId; // ToDo omptimize: write only when changed
                }
                if (edges[curEdgeId].weight < curBestWeight) {
                    curBestWeight = edges[curEdgeId].weight;
                    curBestEid = curEdgeId;
                }

                prevComp = curComp;
                curComp = FlData::list[curComp].second;
                assert(prevComp != curComp);
            }
            FlData::listEnd[v] = prevComp;
            //E(v); Eo(curBestEid);
            bestEid[v] = curBestEid;
        } else {
            bestEid[v] = -1;
        }

        times[iterationNumber][threadId][0] = rdtsc.end(threadId);
#pragma omp barrier

        //
        // renum component
        //
        rdtsc.start(threadId);

#pragma omp for
        for (vid_t i = 0; i < vertexCount; ++i) {
            gComp[i].clear();
            fullComp[i].clear();
        }

        for (int i = 0; i < threadsCount; ++i) {
            graph_messages[threadId][i].clear();
#ifdef USE_SMALL_VECTOR
            if (iterationNumber == 0)
                graph_messages[threadId][i].reserve(vertexCount * 4 / 2 / threadsCount / threadsCount);
#endif
        }
#pragma omp barrier
        const vid_t ito = vertexIds[threadId + 1];
        for (vid_t i = vertexIds[threadId]; i < ito; ++i) {
            const eid_t curBestEid = bestEid[i];
            if (curBestEid == -1) continue;
            const vid_t toi = (iterationNumber == 0 ? edges[curBestEid].dest : comp[edges[curBestEid].dest]); // TODO improve pref using destComp from ExtEdge
            const vid_t toi_thread = toi / vertexesPerThread;
            graph_messages[threadId][toi_thread].push_back(pvv(toi, i));
            gComp[i].push_back(toi);
        }
#pragma omp barrier
        for (int i = 0; i < threadsCount; ++i)
            for (const pvv e : graph_messages[i][threadId])
                gComp[e.first].push_back(e.second);
        times[iterationNumber][threadId][1] = rdtsc.end(threadId);
#pragma omp barrier

        rdtsc.start(threadId);
        {
            eid_t edgesByThreads[threadsCount];
            memset(edgesByThreads, 0, sizeof(edgesByThreads));
            eid_t edgesCountOnNextIter = 0;

            signed char *visited = bfs_visited[threadId];
            for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) if (!gComp[i].empty() && visited[i] < iterationNumber) {
                // bfs
                vector<vid_t> bfs_component(1, i);
                {
                    vid_t *que = bfs_queue[threadId];
                    int from = 0, to = 1;
                    que[from] = i;
                    visited[i] = iterationNumber;
                    while (from < to) {
                        const vid_t wave_from = que[from++];
                        for (int64_t wave_to_id = 0; wave_to_id < gComp[wave_from].size(); ++wave_to_id) {
                            vid_t wave_to = gComp[wave_from].at(wave_to_id);
                        //for (vid_t wave_to : gComp[wave_from]) 
                            if (visited[wave_to] < iterationNumber) {
                                visited[wave_to] = iterationNumber;
                                //assert(to < max_bfs_queue_len);
                                que[to++] = wave_to;
                                bfs_component.push_back(wave_to);
                            }
                        }
                    }
                }
                sort(bfs_component.begin(), bfs_component.end());
                vid_t sum = std::accumulate(bfs_component.begin(), bfs_component.end(), vid_t(0));
                vid_t transfer = bfs_component[sum % bfs_component.size()];
                if (transfer / vertexesPerThread != threadId) continue;

                if (bfs_component.size() > 1)
                    updated = 1;
                fullComp[i].swap(bfs_component);

                set<pvv> compEdges;
                for (vid_t jc : fullComp[i])
                    for (int64_t to_id = 0; to_id < gComp[jc].size(); ++to_id) {
                        vid_t to = gComp[jc].at(to_id);
                    //for (vid_t to : gComp[jc])
                        if (jc < to && !compEdges.count(pvv(jc, to))) {
                            compEdges.insert(pvv(jc, to));
                            if (comp[edges[bestEid[jc]].dest] == to)
                                tmpTaskResult += edges[bestEid[jc]].weight;
                            else
                                tmpTaskResult += edges[bestEid[to]].weight;
                        }
                    }
                //assert(transfer != i);
                for (vid_t j : fullComp[i])
                    comp[j] = transfer;

                if (i != transfer) fullComp[i].swap(fullComp[transfer]);
            }
        }
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);
    }

#ifdef USE_BARANCING_disabled
    // barancing
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);
        rdtsc.start(threadId);

        eid_t edgesWorkPerThreads[threadsCount];
        memset(edgesWorkPerThreads, 0, sizeof(edgesWorkPerThreads));
        vector<pvv> swaps;

        const vid_t vto = vertexIds[threadId + 1];
        for (vid_t i = vertexIds[threadId]; i < vto; ++i) if (!fullComp[i].empty()) {
            // fullComp contains many vertexes from one thread
            // time may be reduced if generate something line fullCompThreads, which will
            // store only unique threads ids
            random_shuffle(fullComp[i].begin(), fullComp[i].end());
            vid_t moveWorkTo = *min_element(fullComp[i].begin(), fullComp[i].end(), [&](vid_t a, vid_t b) -> bool {
                    int ta = a / vertexesPerThread;
                    int tb = b / vertexesPerThread;
                    return edgesWorkPerThreads[ta] < edgesWorkPerThreads[tb];
                    });
            int moveThreadId = moveWorkTo / vertexesPerThread;
            eid_t compEdgesCount = accumulate(fullComp[i].begin(), fullComp[i].end(), eid_t(0), [&](eid_t sum, vid_t v) -> eid_t {
                        const int vThreadId = v / vertexesPerThread;
                        return sum + edgesIdsByThread[vThreadId][v + 1] - edgesIdsByThread[vThreadId][v];
                        });
            if (compEdgesCount > 1000000) {
                Eo(compEdgesCount);
            }
            if (moveThreadId != threadId) {
                swaps.push_back(make_pair(i, moveWorkTo));
                comp[i] = moveWorkTo;
                comp[moveWorkTo] = moveWorkTo;
                edgesWorkPerThreads[moveThreadId] += compEdgesCount;
            } else {
                edgesWorkPerThreads[threadId] += compEdgesCount;
            }
        }

#pragma omp barrier
        for (const pvv& i : swaps)
            fullComp[i.first].swap(fullComp[i.second]);

        times[iterationNumber][threadId][4] += rdtsc.end(threadId);
    }
#endif

#pragma omp barrier // ???

    // pointer jumping
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
                    comp[i] = min(i, myComp);
                    changed = 1;
                } else if (myComp != parentComp) {
                    comp[i] = parentComp;
                    changed = 1;
                }
            }
            times[iterationNumber][threadId][3] += rdtsc.end(threadId);
        }
    }
    

    // merge edges lists
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);
        rdtsc.start(threadId);

        for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) { // iterate over all vertexes in this thread
            if (comp[v] == v && fullComp[v].size() > 1) { // if this vertex represents component and have edges to another component
                // union components
                for (vid_t u : fullComp[v]) if (u != v) { // iterate over all vertexes in new component
                    assert(FlData::listEnd[v] != -1);
                    assert(FlData::list[FlData::listEnd[v]].second == -1);
                    assert(FlData::listEnd[u] != -1);
                    FlData::list[FlData::listEnd[v]].second = FlData::listBegin[u];

                    assert(FlData::list[FlData::listEnd[u]].second == -1);
                    FlData::listEnd[v] = FlData::listEnd[u];
                    assert(FlData::listEnd[v] != -1);
                    assert(FlData::list[FlData::listEnd[v]].second == -1);
                }
            }
        }

        times[iterationNumber][threadId][5] = rdtsc.end(threadId);
#pragma omp barrier

    }

    taskResult += tmpTaskResult;
    ++iterationNumber;
    return updated;
}

void doPrepare() {
    comp = new vid_t[vertexCount];
    for (vid_t i = 0; i < vertexCount; ++i)
        comp[i] = i;

#ifdef USE_SMALL_VECTOR
    gComp = new Vector<vid_t, 100>[vertexCount];
#else
    gComp = new vector<vid_t>[vertexCount];
#endif
    fullComp = new vector<vid_t>[vertexCount];
    bestEid = new eid_t[vertexCount];

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        const int curThreadsCount = omp_get_num_threads();
#pragma omp master
        {
            threadsCount = curThreadsCount;

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
            FlData::list[i] = pvv(i, -1);
            FlData::listBegin[i] = i;
            FlData::listEnd[i] = i;
            FlData::startEdge[i] = edgesIds[i];

            std::sort(edges + edgesIds[i], edges + edgesIds[i + 1], EdgeWeightCmp());
        }

        // bfs data
        bfs_visited[threadId] = new signed char[vertexCount];
        std::fill(bfs_visited[threadId], bfs_visited[threadId] + vertexCount, -1);
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
            for (int k = 0; k < 6; ++k)
                fprintf(stderr, "%.6lf ", times[i][j][k]);
            fputs("\n", stderr);
        }
    }
#endif

    return 0;
}
