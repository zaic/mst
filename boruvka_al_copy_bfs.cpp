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
const int kMaxIterations = 40;
const int kMaxThreads = 20;

vid_t *comp;
weight_t taskResult;

int threadsCount;
#ifdef ON_NUMA
vid_t vertexesPerThread;
#else
#define vertexesPerThread #error numa_only
#endif
vid_t *vertexIds;
int iterationNumber;

#ifdef USE_EDGE_STRUCT
struct ExtEdge {
    weight_t weight; // 8
    eid_t origId;    // 8
    vid_t destComp;  // 4

    bool operator<(const ExtEdge& o) const {
        if (weight != o.weight) return weight < o.weight;
        return destComp < o.destComp;
    }
};
ExtEdge **edgesByThread;
#else
namespace ExtEdge {
    weight_t *weight[kMaxThreads];
    eid_t    *origId[kMaxThreads];
    vid_t    *destComp[kMaxThreads];

    inline bool compareEdge(int ta, eid_t a, int tb, eid_t b) {
        const weight_t wa = weight[ta][a];
        const weight_t wb = weight[tb][b];
        if (wa != wb) return wa < wb;
        return destComp[ta][a] < destComp[tb][b];
    }

    void allocEdges() {
        for (int i = 0; i < threadsCount; ++i) {
            Eo(edgesCount);
            weight[i] = (weight_t*)malloc(edgesCount * sizeof(weight_t));
            origId[i] = (eid_t*)malloc(edgesCount * sizeof(eid_t));
            destComp[i] = (vid_t*)malloc(edgesCount * sizeof(vid_t));
        }
    }
} /* namespace ExtEdge */
#endif /* USE_EDGE_STRUCT */

double times[kMaxThreads][kMaxIterations][40];


/*
 *  temporary dfs data
 */
#ifdef USE_SMALL_VECTOR
Vector<vid_t, 100> *gComp;
#else
vector<vid_t> *gComp; // adjacency list
#endif
vector<pvv> graph_messages[kMaxThreads][kMaxThreads];

vector<vid_t> *fullComp; // contains entire new component
eid_t *bestEid; // TODO very temporary array
/*
 *  temporary bfs data
 */
const int kMaxBfsQueue = 10000;
vid_t bfs_queue[kMaxThreads][kMaxBfsQueue];
signed char *bfs_visited[kMaxThreads];

/*
 *  edges lists
 */
eid_t *sumOfAllEdges; // TODO make 2-dim: [arrayFrom][arrayTo + gap]
eid_t **edgesIdsByThread; // indexes in previous array




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
            gComp[i].pushBack(toi);
        }
#pragma omp barrier
        for (int i = 0; i < threadsCount; ++i)
            for (const pvv e : graph_messages[i][threadId])
                gComp[e.first].pushBack(e.second);
        times[iterationNumber][threadId][1] = rdtsc.end(threadId);
#pragma omp barrier

        rdtsc.start(threadId);
//#pragma omp master
        {
            eid_t edgesByThreads[threadsCount];
            memset(edgesByThreads, 0, sizeof(edgesByThreads));
            eid_t edgesCountOnNextIter = 0;

            signed char *visited = bfs_visited[threadId];
            //for (vid_t i = 0; i < vertexCount; ++i) if (!gComp[i].empty() && visited[i] < iterationNumber) {
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

            sumOfAllEdges[threadId] = edgesCountOnNextIter;
        }
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);
    }

#ifdef USE_BARANCING
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
#if 0
#pragma omp critical
        {
        E(threadId); Eo(edgesWorkPerThreads[0]);
        E(threadId); Eo(edgesWorkPerThreads[1]);
        E(threadId); Eo(edgesWorkPerThreads[2]);
        E(threadId); Eo(edgesWorkPerThreads[3]);
        }
#endif

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

        eid_t sumFutureEdges = 0;
        // for (int i = 0; i < threadsCount; ++i) sumFutureEdges += sumOfAllEdges[i][threadId];
#ifdef USE_EDGE_STRUCT
        //ExtEdge *nextIterEdges = new ExtEdge[edgesCount]; // TODO fix: using sumFutureEdges
        ExtEdge *nextIterEdges = (ExtEdge*)malloc(sizeof(ExtEdge) * edgesCount);
#else
        weight_t *nextWeight   = (weight_t*)malloc(edgesCount * sizeof(weight_t));
        eid_t    *nextOrigId   = (eid_t*)   malloc(edgesCount * sizeof(eid_t));
        vid_t    *nextDestComp = (vid_t*)   malloc(edgesCount * sizeof(vid_t));
#endif
        //eid_t *nextEdgesIds = new eid_t[vertexCount + 1]; // TODO array size can be reduced
        eid_t *nextEdgesIds = (eid_t*)malloc(sizeof(eid_t) * (vertexCount + 1));
        nextEdgesIds[0] = 0;
        eid_t edgeId = 0;

        for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) { // iterate over all vertexes in this thread
            if (comp[v] == v && fullComp[v].size() > 1) { // if this vertex represents component and have edges to another component
                // copy all edges
                for (vid_t u : fullComp[v]) { // iterate over all vertexes in new component
                    int ownedThread = u / vertexesPerThread;
                    for (eid_t e = edgesIdsByThread[ownedThread][u]; e < edgesIdsByThread[ownedThread][u + 1]; ++e) {
#ifdef USE_EDGE_STRUCT
                        if (comp[edgesByThread[ownedThread][e].destComp] == v) continue; // skip loops
                        nextIterEdges[edgeId] = edgesByThread[ownedThread][e];
#else
                        const vid_t edgeDest = ExtEdge::destComp[ownedThread][e];
                        if (comp[edgeDest] == v) continue; // skip loops
                        const eid_t edgeOrigId = ExtEdge::origId[ownedThread][e];
                        const weight_t edgeWeidght = ExtEdge::weight[ownedThread][e];

                        nextWeight[edgeId] = edgeWeidght;
                        nextOrigId[edgeId] = edgeOrigId;
                        nextDestComp[edgeId] = edgeDest;
#endif
                        ++edgeId;
                    }
                }
            }
            nextEdgesIds[v + 1] = edgeId;
        }

        times[iterationNumber][threadId][5] = rdtsc.end(threadId);
#pragma omp barrier

#ifdef USE_EDGE_STRUCT
        free(edgesByThread[threadId]);
        edgesByThread[threadId] = nextIterEdges;
#else
        free(ExtEdge::weight[threadId]);
        ExtEdge::weight[threadId] = nextWeight;
        free(ExtEdge::origId[threadId]);
        ExtEdge::origId[threadId] = nextOrigId;
        free(ExtEdge::destComp[threadId]);
        ExtEdge::destComp[threadId] = nextDestComp;
#endif
        free(edgesIdsByThread[threadId]);
        edgesIdsByThread[threadId] = nextEdgesIds;
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

            sumOfAllEdges = new eid_t[threadsCount];
            edgesIdsByThread = new eid_t*[threadsCount];
#ifdef USE_EDGE_STRUCT
            edgesByThread = new ExtEdge*[threadsCount];
#else
            ExtEdge::allocEdges();
#endif /* USE_EDGE_STRUCT */
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

        sumOfAllEdges[threadId] = 0;
        edgesIdsByThread[threadId] = new eid_t[vertexCount + 1]; // TODO: it's possible to reduce size of this array
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            const eid_t curEdgesCount = edgesIds[i + 1] - edgesIds[i];
            sumOfAllEdges[threadId] += curEdgesCount;
        }
#ifdef USE_EDGE_STRUCT
        edgesByThread[threadId] = new ExtEdge[sumOfAllEdges[threadId]];
        auto myEdges = edgesByThread[threadId];
#endif
        eid_t myEdgeId = 0;
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            edgesIdsByThread[threadId][i] = myEdgeId;
            for (eid_t j = edgesIds[i]; j < edgesIds[i + 1]; ++j) if (edges[j].dest != i) {
#ifdef USE_EDGE_STRUCT
                myEdges[myEdgeId] = ExtEdge{edges[j].weight, j, edges[j].dest};
#else
                ExtEdge::weight[threadId][myEdgeId] = edges[j].weight;
                ExtEdge::origId[threadId][myEdgeId] = j;
                ExtEdge::destComp[threadId][myEdgeId] = edges[j].dest;
#endif
                ++myEdgeId;
            }
        }
        edgesIdsByThread[threadId][vertexIds[threadId + 1]] = myEdgeId;

        // bfs data
        bfs_visited[threadId] = new signed char[vertexCount];
        for (vid_t i = 0; i < vertexCount; ++i)
            bfs_visited[threadId][i] = -1;
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
