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
using namespace std;

typedef pair<weight_t, eid_t> pwe;

vid_t *comp;
weight_t taskResult;

int threadsCount;
vid_t *vertexIds;
int iterationNumber;

struct ExtEdge {
    weight_t weight;
    eid_t origId;
    vid_t destComp;

    bool operator<(const ExtEdge& o) const {
        if (weight != o.weight) return weight < o.weight;
        return destComp < o.destComp;
    }
};

const int kMaxIterations = 40;
const int kMaxThreads = 20;
double times[kMaxThreads][kMaxIterations][40];


/*
 *  temporary dfs data
 */
vector<vid_t> *gComp; // adjacency list
vector<vid_t> *fullComp; // contains entire new component
eid_t *bestEid; // TODO very temporary array
/*
 *  temporary bfs data
 */
const int max_bfs_queue_len = 1000;
vid_t bfs_queue[max_bfs_queue_len];
signed char *bfs_visited;

/*
 *  edges lists
 */
eid_t *sumOfAllEdges; // TODO make 2-dim: [arrayFrom][arrayTo + gap]
ExtEdge **edgesByThread;
eid_t **edgesIdsByThread; // indexes in previous array




bool doAll() {
    int updated = 0; // reduce stage 
    weight_t tmpTaskResult = 0.0; // merge stage
    Eo(iterationNumber);

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);

        //
        // find minimum edge
        //
        rdtsc.start(threadId);

        auto curEdges = edgesByThread[threadId];
        auto curEdgesIds = edgesIdsByThread[threadId];
        for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) {
            // E(v); E(curEdgesIds[v]); Eo(curEdgesIds[v + 1]);
            if (curEdgesIds[v] >= curEdgesIds[v + 1]) { // there is no edges from vertex v on this thread
                bestEid[v] = -1;
                continue;
            }
            assert(comp[v] == v);
            eid_t bestId = curEdgesIds[v];
            for (eid_t e = curEdgesIds[v] + 1; e < curEdgesIds[v + 1]; ++e)
                if (curEdges[e] < curEdges[bestId])
                    bestId = e;
            bestEid[v] = curEdges[bestId].origId;
            // if (iterationNumber > 0) { E(v); E(curEdgesIds[v]); Eo(curEdgesIds[v + 1]); Eo(bestEid[v]); }
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

#pragma omp for
        for (vid_t i = 0; i < vertexCount; ++i) if (comp[i] == i && bestEid[i] != -1) {
            vid_t toi = comp[edges[bestEid[i]].dest]; // TODO improve pref using destComp from ExtEdge
            //assert(comp[toi] != comp[i]);
            gComp[i].push_back(toi);
            //gComp[toi].push_back(i);
        }
#pragma omp master
        {
            // build graph where vertexes are components, which were constructed on prev iteration
            for (vid_t i = 0; i < vertexCount; ++i) if (comp[i] == i && bestEid[i] != -1) {
                vid_t toi = comp[edges[bestEid[i]].dest]; // TODO improve pref using destComp from ExtEdge
                //assert(comp[toi] != comp[i]);
                //gComp[i].push_back(toi);
                gComp[toi].push_back(i);
            }
        }
        times[iterationNumber][threadId][1] = rdtsc.end(threadId);
#pragma omp barrier

        rdtsc.start(threadId);
#pragma omp master
        {
            eid_t edgesByThreads[threadsCount];
            memset(edgesByThreads, 0, sizeof(edgesByThreads));

            for (vid_t i = 0; i < vertexCount; ++i) if (!gComp[i].empty() && bfs_visited[i] < iterationNumber) {
                // bfs
                {
                    int from = 0, to = 1;
                    bfs_queue[from] = i;
                    bfs_visited[i]  = iterationNumber;
                    vector<vid_t> bfs_component(1, i);
                    while (from < to) {
                        const vid_t wave_from = bfs_queue[from++];
                        for (vid_t wave_to : gComp[wave_from]) if (bfs_visited[wave_to] < iterationNumber) {
                            bfs_visited[wave_to] = iterationNumber;
                            //assert(to < max_bfs_queue_len);
                            bfs_queue[to++] = wave_to;
                            bfs_component.push_back(wave_to);
                        }
                    }
                    //fullComp[i] = std::move(vector<vid_t>(bfs_queue, bfs_queue + to));
                    fullComp[i].swap(bfs_component);
                }

                // old part begins
                if (fullComp[i].size() > 1)
                    updated = 1;
                vector<vid_t> theirComps;
                for (int cur : fullComp[i]) {
                    int ownedThread = 0;
                    while (vertexIds[ownedThread + 1] <= cur) ++ownedThread;
                    theirComps.push_back(ownedThread);
                }
                //E(i); Eo(fullComp[i].size());

                const eid_t sumEdgesCount = accumulate(fullComp[i].begin(), fullComp[i].end(), 0, [](eid_t res, vid_t cur) -> eid_t { 
                    int ownedThread = 0;
                    while (vertexIds[ownedThread + 1] <= cur) ++ownedThread;
                    return res + edgesIdsByThread[ownedThread][cur + 1] - edgesIdsByThread[ownedThread][cur];
#if 0
                        int curThreadId = 0;
                        for (; vertexIds[curThreadId] > cur; ++curThreadId);
                        //return res + sumOfAllEdges;  TODO calc edges count for vertex "cur"
#endif
                });

                //vid_t transfer = fullComp[i][(i * 37 + iterationNumber * 47) % fullComp[i].size()];
                vid_t transfer = 0;
                for (int nt = 1; nt < fullComp[i].size(); ++nt) if (edgesByThreads[theirComps[nt]] < edgesByThreads[theirComps[transfer]])
                    transfer = nt;
                edgesByThreads[theirComps[transfer]] += sumEdgesCount;
                transfer = fullComp[i][transfer];

                /*
                int bestThreadId = upper_bound(vertexIds, vertexIds + threadsCount + 1, transfer) - vertexIds - 1;
                for (vid_t j : fullComp[i]) {
                    const int qthreadId = lower_bound(vertexIds, vertexIds + threadsCount + 1, j) - vertexIds - 1;
                    if (edgesByThreads[qthreadId] < edgesByThreads[bestThreadId]) {
                        transfer = j;
                        bestThreadId = qthreadId;
                    }
                }
                */
                //edgesByThreads[bestThreadId] += sumEdgesCount;

                set<pvv> compEdges;
                for (vid_t jc : fullComp[i])
                    for (vid_t to : gComp[jc])
                        if (jc < to && !compEdges.count(pvv(jc, to))) {
                            compEdges.insert(pvv(jc, to));
                            if (comp[edges[bestEid[jc]].dest] == to)
                                tmpTaskResult += edges[bestEid[jc]].weight;
                            else
                                tmpTaskResult += edges[bestEid[to]].weight;
                        }
                //assert(transfer != i);
                for (vid_t j : fullComp[i])
                    comp[j] = transfer;

                if (i != transfer) fullComp[i].swap(fullComp[transfer]);
            }

            //for (int i = 0; i < threadsCount; ++i) { E(i); Eo(edgesByThreads[i]); }
        }
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);

#if 0
#pragma omp for reduction(+:tmpTaskResult) nowait
        for (vid_t i = 0; i < vertexCount; ++i) {
            if (adjacencyLists[i].listSize == 0) continue;

            vid_t oc = adjacencyLists[i].edges[0].destComp;

            if (comp[oc] == i) {
                if (i < oc) {
                    comp[i] = i;
                } else {
                    comp[i] = oc;
                    tmpTaskResult += adjacencyLists[i].edges[0].weight;
                }
            } else {
                tmpTaskResult += adjacencyLists[i].edges[0].weight;
                comp[i] = oc;
            }
        }
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);
#endif
    }

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
            times[iterationNumber][threadId][3] = rdtsc.end(threadId);
        }
    }

    // merge edges lists
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId * 2);
        rdtsc.start(threadId);

        eid_t sumFutureEdges = 0;
        // for (int i = 0; i < threadsCount; ++i) sumFutureEdges += sumOfAllEdges[i][threadId];
        ExtEdge *nextIterEdges = new ExtEdge[edgesCount]; // TODO fix: using sumFutureEdges
        eid_t *nextEdgesIds = new eid_t[vertexCount + 1]; // TODO array size can be reduced
        nextEdgesIds[0] = 0;
        eid_t edgeId = 0;

        for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) { // iterate over all vertexes in this thread
            if (comp[v] == v && fullComp[v].size() > 1) { // if this vertex represents component and have edges to another component
                // copy all edges
                for (vid_t u : fullComp[v]) { // iterate over all vertexes in new component
                    int ownedThread = 0;
                    while (vertexIds[ownedThread + 1] <= u) ++ownedThread;
                    for (eid_t e = edgesIdsByThread[ownedThread][u]; e < edgesIdsByThread[ownedThread][u + 1]; ++e) {
                        if (comp[edgesByThread[ownedThread][e].destComp] == v) continue; // skip loops
                        nextIterEdges[edgeId] = edgesByThread[ownedThread][e];
                        ++edgeId;
                    }
                }
            }
            nextEdgesIds[v + 1] = edgeId;
        }

        times[iterationNumber][threadId][4] = rdtsc.end(threadId);
#pragma omp barrier

        delete[] edgesByThread[threadId];
        edgesByThread[threadId] = nextIterEdges;
        delete[] edgesIdsByThread[threadId];
        edgesIdsByThread[threadId] = nextEdgesIds;
    }

    taskResult += tmpTaskResult;
    ++iterationNumber;
    return updated;
}

void doPrepare() {
    vertexIds = new vid_t[vertexCount + 1]; // TODO fix: array size is threads_count + 1
    vertexIds[0] = 0;

    comp = new vid_t[vertexCount];
    for (vid_t i = 0; i < vertexCount; ++i)
        comp[i] = i;

    gComp = new vector<vid_t>[vertexCount];
    fullComp = new vector<vid_t>[vertexCount];
    bestEid = new eid_t[vertexCount];

    // bfs data
    bfs_visited = new signed char[vertexCount];
    for (vid_t i = 0; i < vertexCount; ++i)
        bfs_visited[i] = -1;

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        const int curThreadsCount = omp_get_num_threads();
#pragma omp master
        {
            threadsCount = curThreadsCount;
            sumOfAllEdges = new eid_t[threadsCount];
            edgesByThread = new ExtEdge*[threadsCount];
            edgesIdsByThread = new eid_t*[threadsCount];
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
#endif
#pragma omp barrier

        sumOfAllEdges[threadId] = 0;
        edgesIdsByThread[threadId] = new eid_t[vertexCount + 1]; // TODO: it's possible to reduce size of this array
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            const eid_t curEdgesCount = edgesIds[i + 1] - edgesIds[i];
            sumOfAllEdges[threadId] += curEdgesCount;
        }
        edgesByThread[threadId] = new ExtEdge[sumOfAllEdges[threadId]];
        auto myEdges = edgesByThread[threadId];
        eid_t myEdgeId = 0;
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            edgesIdsByThread[threadId][i] = myEdgeId;
            for (eid_t j = edgesIds[i]; j < edgesIds[i + 1]; ++j) if (edges[j].dest != i) {
                myEdges[myEdgeId] = ExtEdge{edges[j].weight, j, edges[j].dest};
                if (myEdgeId % 100000 == 0) {
                    //E(myEdgeId); E(i); E(j); E(edgesByThread[0][myEdgeId].weight); E(edgesByThread[0][myEdgeId].origId); Eo(edgesByThread[0][myEdgeId].weight);
                }
                ++myEdgeId;
            }
        }
        Eo(myEdgeId);
        edgesIdsByThread[threadId][vertexIds[threadId + 1]] = myEdgeId;
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
            for (int k = 0; k < 5; ++k)
                fprintf(stderr, "%.6lf ", times[i][j][k]);
            fputs("\n", stderr);
        }
    }
#endif

    return 0;
}
