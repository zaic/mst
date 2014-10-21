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

struct EdgeList {
    ExtEdge *edges;
    eid_t listSize;

    EdgeList() {

    }

    EdgeList(eid_t len, int cpuLocation = 0) {
        edges = new ExtEdge[len];
        listSize = len;
    }

    void clean() {
        delete[] edges;
        edges = NULL;
        listSize = 0;
    }

    ~EdgeList() {
        clean();
    }

    EdgeList& operator=(EdgeList&& o) {
        if (this != &o) {
            clean();
            swap(edges, o.edges);
            swap(listSize, o.listSize);
        }
        return *this;
    }

    void merge(const EdgeList& listOne, const EdgeList& listTwo) {
        eid_t posOne = 0, posTwo = 0;
        // ToDo filter edges
        while (posOne < listOne.listSize && posTwo < listTwo.listSize) {
            edges[listSize++] = (listOne.edges[posOne] < listTwo.edges[posTwo]) ?
                listOne.edges[posOne++] :
                listTwo.edges[posTwo++];
        }
        while (posOne < listOne.listSize)
            edges[listSize++] = listOne.edges[posOne++];
        while (posTwo < listTwo.listSize)
            edges[listSize++] = listTwo.edges[posTwo++];
    }
};

EdgeList *adjacencyLists;

vector<vid_t> *gComp;
bool *visitComp;
vector<vid_t> *fullComp;
void dfs(vid_t v, vector<vid_t>& tmpComp) {
    visitComp[v] = true;
    tmpComp.push_back(v);
    for (vid_t u : gComp[v]) if (!visitComp[u]) dfs(u, tmpComp);
}

bool doAll() {
    int updated = 0; // reduce stage 
    weight_t tmpTaskResult = 0.0; // merge stage
    Eo(iterationNumber);

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);

        //
        // renum component
        //
        rdtsc.start(threadId);

#pragma omp for
        for (vid_t i = 0; i < vertexCount; ++i) {
            gComp[i].clear();
            fullComp[i].clear();
        }
#pragma omp master
        {
            for (vid_t i = 0; i < vertexCount; ++i) if (comp[i] == i && adjacencyLists[i].listSize > 0) {
                vid_t toi = adjacencyLists[i].edges[0].destComp;
                assert(comp[toi] != comp[i]);
                gComp[i].push_back(toi);
                gComp[toi].push_back(i);
            }
            memset(visitComp, 0, vertexCount);
            for (vid_t i = 0; i < vertexCount; ++i) if (!gComp[i].empty() && !visitComp[i]) {
                dfs(i, fullComp[i]);
                for (vid_t j : fullComp[i])
                    comp[j] = fullComp[i].back();
                if (fullComp[i].size() > 1)
                    updated = 1;
                //E(i); Eo(fullComp[i].size());
                set<pvv> compEdges;
                for (vid_t jc : fullComp[i])
                    for (vid_t to : gComp[jc])
                        if (jc < to && !compEdges.count(pvv(jc, to))) {
                            compEdges.insert(pvv(jc, to));
                            if (adjacencyLists[jc].edges[0].destComp == to)
                                tmpTaskResult += adjacencyLists[jc].edges[0].weight;
                            else
                                tmpTaskResult += adjacencyLists[to].edges[0].weight;
                        }
                vid_t transfer = fullComp[i].back();
                assert(transfer != i);
                fullComp[i].swap(fullComp[transfer]);
            }
        }
        times[iterationNumber][threadId][0] = rdtsc.end(threadId);

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
            times[iterationNumber][threadId][1] = rdtsc.end(threadId);
        }
    }

    // merged edges lists
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);
        rdtsc.start(threadId);
#pragma omp for
        for (vid_t i = 0; i < vertexCount; ++i) if (comp[i] == i) {
            //Eo(i);
            vector<eid_t> curPos(fullComp[i].size(), 0); // position in edges lists for each merged component
            eid_t have = accumulate(fullComp[i].begin(), fullComp[i].end(), 0, [](eid_t a, eid_t b) ->  eid_t { return a + adjacencyLists[b].listSize; }); // sum of all merged edges
            EdgeList mergedList(have); // new list
            eid_t copyTo = 0; // position in new list
            //Eo(have); Eo(fullComp[i].size());

            while (have > 0) {
                vid_t best = vertexCount; // index in fullComp and curPos vectors
                eid_t bestPos = 0; // index in edges list for best component
                vid_t bestComp = 0; // real vertex index for best component
                for (vid_t j = 0; j < fullComp[i].size(); ++j) {
                    const vid_t jc = fullComp[i][j];
                    const eid_t ps = curPos[j];
                    if (ps < adjacencyLists[jc].listSize && (best == vertexCount || adjacencyLists[jc].edges[ps] < adjacencyLists[bestComp].edges[bestPos])) {
                        best = j;
                        bestPos = ps;
                        bestComp = jc;
                    }
                }
                //E(best); E(vertexCount); Eo(have);
                assert(best < vertexCount);
                --have;

                bool ok = true;
                for (vid_t jc : fullComp[i]) // check for loops
                    if (adjacencyLists[bestComp].edges[bestPos].destComp == jc) {
                        ok = false;
                        break;
                    }
                if (ok) {
                    mergedList.edges[copyTo] = adjacencyLists[bestComp].edges[bestPos];
                    mergedList.edges[copyTo].destComp = comp[mergedList.edges[copyTo].destComp];
                    copyTo++;
                    // TODO fix dest comp
                }
                curPos[best]++;
            }
            mergedList.listSize = copyTo;
            // remove old lists
            for (vid_t jc : fullComp[i])
                adjacencyLists[jc].clean();
            adjacencyLists[comp[i]] = move(mergedList);
        }
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);
    }

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

    gComp = new vector<vid_t>[vertexCount];
    visitComp = new bool[vertexCount];
    fullComp = new vector<vid_t>[vertexCount];

    adjacencyLists = new EdgeList[vertexCount];

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        const int curThreadsCount = omp_get_num_threads();
#pragma omp master
        {
            threadsCount = curThreadsCount;
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

        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            sort(edges + edgesIds[i], edges + edgesIds[i + 1], EdgeWeightCmp());
            adjacencyLists[i] = EdgeList(edgesIds[i + 1] - edgesIds[i]);
            eid_t copyTo = 0;
            for (eid_t j = 0; j < adjacencyLists[i].listSize; ++j) {
                const eid_t eid = edgesIds[i] + j;
                if (edges[eid].dest == i) continue;
                adjacencyLists[i].edges[copyTo++] = ExtEdge{edges[eid].weight, eid, edges[eid].dest};
            }
            adjacencyLists[i].listSize = copyTo;
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
            for (int k = 0; k < 3; ++k)
                fprintf(stderr, "%.6lf ", times[i][j][k]);
            fputs("\n", stderr);
        }
    }
#endif

    return 0;
}
