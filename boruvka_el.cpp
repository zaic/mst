#include "gen.h"
#include <cmath>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <set>
#include <tuple>
#include <map>
#include <vector>
#include <numeric>
#include "stat.h"
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
    vid_t destComp;
#ifdef USE_RESULT_VERTEX
    vid_t from;
#else
    eid_t edgeId;
#endif

    bool operator<(const Result& o) const {
        return weight < o.weight;
    }
};

Result **localResult, *bestResult;

double times[kMaxIterations][kMaxThreads][40];

bool haveOuterComps[kMaxThreads];
#ifdef USE_SKIP_LOOPS
bool useSkipLoops = true;
#endif /* USE_SKIP_LOOPS */
eid_t loopsViewEdges, loopsSkipEdges; // USE_SKIP_LOOPS
#ifdef USE_SKIP_LAST_ITER
vid_t aliveComponents;
#endif /* USE_SKIP_LAST_ITER */
#ifdef USE_COMPRESS
int64_t *generatedComps;
vid_t lastUsedVid;
vid_t prevUsedVid;
#endif /* USE_COMPRESS */
#ifdef USE_ANSWER_IN_VECTOR
namespace Answer {
    eid_t **data;
    vid_t **pos;

    void init() {
        data = new eid_t*[threadsCount];
        pos = new vid_t*[threadsCount];
#pragma omp parallel
        {
            const int threadId = omp_get_thread_num();
            data[threadId] = static_cast<eid_t*>(malloc(sizeof(eid_t) * vertexCount));
            pos[threadId] = static_cast<vid_t*>(malloc(sizeof(vid_t) * 16));
            pos[threadId][0] = 0;
            const int pageSize = 4 * 1024;
            const int elementsPerPage = pageSize / sizeof(eid_t);
            for (vid_t i = vertexCount / threadsCount - 1; i >= 0; i -= elementsPerPage)
                data[threadId][i] = 0;
        }
    }

    void reset() {
        for (int i = 0; i < threadsCount; ++i) pos[i][0] = 0;
    }
}
#endif /* USE_ANSWER_IN_VECTOR */
vid_t compBound[kMaxThreads][16];
pvv *vertexBound;

//Stat<vid_t> actulVers;
//Stat<eid_t> less010, more010;
Stat<vid_t> diedComps;
Stat<vid_t> catchComps;
Stat<vid_t> potentComps;
vid_t maxWindowSize;
bool useMagic = true;



template<int definedIter>
bool doAll() {
    Eo(iterationNumber);
    int updated = 0; // reduce stage 
    int diedComponents = 0; // merge stage
    weight_t tmpTaskResult = 0.0; // merge stage

    const bool usePrefetchStartEdge = (iterationNumber < 3 || double(vertexCount) / generatedComps[threadsCount - 1] > pow(10.0, iterationNumber));
    if (iterationNumber > 4) useMagic = false; // TODO template parameter
    if (iterationNumber < 2) useMagic = false;
    if (iterationNumber == 2 && maxWindowSize < vertexCount / 2) useMagic = true;
    Eo(useMagic);

#pragma omp parallel reduction(+:diedComponents) reduction(+:loopsViewEdges) reduction(+:loopsSkipEdges)
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);

        //
        // find min
        // 
        rdtsc.start(threadId);
        auto result = (definedIter == 0 ? bestResult : localResult[threadId]);
#ifdef USE_COMPRESS
        const weight_t MAX_DROPPED_WEIGHT = MAX_WEIGHT + 0.3 * iterationNumber;
#else
        const weight_t MAX_DROPPED_WEIGHT = MAX_WEIGHT + 0.1;
#endif /* USE_COMPRESS */


#ifdef USE_SKIP_LOOPS
        if (iterationNumber >= USE_SKIP_LOOPS + 1 && useSkipLoops) {
            for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1] - 1; ++v) {
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
                if (iterationNumber == USE_SKIP_LOOPS + 1) {
                    loopsViewEdges += newEdgesStart - edgesStart + 1;
                    loopsSkipEdges += newEdgesStart - edgesStart;
                }
                if (newEdgesStart != edgesStart) startEdgesIds[v] = newEdgesStart;
            }
        }
        times[iterationNumber][threadId][4] = rdtsc.end(threadId);
#endif /* USE_SKIP_LOOPS */

        rdtsc.start(threadId);
        bool outerComps = false;
        vid_t minv = vertexCount * 2, maxv = 0;
        if (definedIter == 0) {
            const int64_t vkf = 98; // TODO fix on SSCA2
            const vid_t vfrom = vertexIds[threadId] * vkf / 100;
            const vid_t vto = (threadId == threadsCount - 1 ? vertexIds[threadsCount] : vertexIds[threadId + 1] * vkf / 100);
            //for (vid_t v = vertexIds[threadId]; v < vertexIds[threadId + 1]; ++v) {
            for (vid_t v = vfrom; v < vto; ++v) {
                if (usePrefetchStartEdge) {
                    __builtin_prefetch(edges + startEdgesIds[v + PREFETCH_START_EDGE]);
                }
                const eid_t edgesStart = startEdgesIds[v];
                const eid_t edgesEnd = edgesIds[v + 1];

                if (edgesStart >= edgesEnd) {
                    result[v].weight = MAX_WEIGHT + 0.1;
                    comp[v] = v;
                } else {
                    const vid_t u = edges[edgesStart].dest;
                    /*
#ifdef USE_RESULT_VERTEX
                    result[v] = Result{edges[edgesStart].weight, u, v};
#else
                    result[v] = Result{edges[edgesStart].weight, u, edgesStart};
#endif
*/
                    result[v].weight = MAX_DROPPED_WEIGHT;
                    comp[v] = u;
                    assert(comp[v] != v);
                    startEdgesIds[v] = edgesStart + 1; // TODO fix start edge on initialization
                }
            }

        } else {
            const vid_t vBegin = vertexIds[threadId], vEnd = vertexIds[threadId + 1];
            vid_t currCompBegin = vBegin, currCompEnd = vBegin + 1;
            vid_t currCompValue = comp[vBegin];
            if (useMagic) {
                while (currCompBegin > 0 && currCompBegin > currCompEnd - MAGIC_BOUND && comp[currCompBegin - 1] == currCompValue) --currCompBegin;
            }
            for (vid_t v = vBegin; v < vEnd; ++v) {
                if (usePrefetchStartEdge) {
                    __builtin_prefetch(edges + startEdgesIds[v + PREFETCH_START_EDGE]);
                    if (!useMagic) __builtin_prefetch(comp + edges[startEdgesIds[v + PREFETCH_START_EDGE / 2]].dest); // TODO
                }
                if (useMagic) {
                    if (comp[v] != currCompValue) {
                        currCompValue = comp[v];
                        currCompBegin = v;
                        currCompEnd = v + 1;
                    }
                    while (currCompEnd < currCompBegin + MAGIC_BOUND && currCompEnd < vertexCount && comp[currCompEnd] == currCompValue)
                        ++currCompEnd;
                }

                const eid_t edgesStart = startEdgesIds[v];
                const eid_t edgesEnd = edgesIds[v + 1];
                if (edgesStart >= edgesEnd) continue;

                const vid_t cv = comp[v];
                minv = std::min(minv, cv);
                maxv = std::max(maxv, cv);
                //assert(prevUsedVid <= v && v < lastUsedVid);
                eid_t newEdgesStart = edgesStart;
#ifdef USE_FAST_REDUCTION
                if (cv < vertexIds[threadId] || cv >= vertexIds[threadId + 1]) outerComps = true;
#endif /* USE_FAST_REDUCTION */

                if (useMagic) {
                    if (currCompBegin <= vertexBound[v].first && vertexBound[v].second < currCompEnd) {
                        //catchComps.add(threadId, 1);
                        startEdgesIds[v] = edgesEnd;
                        continue;
                    }
                }

                if (useMagic) {
                    for (; newEdgesStart < edgesEnd; ++newEdgesStart) {
                        const vid_t u = edges[newEdgesStart].dest;
                        if (currCompBegin <= u && u < currCompEnd) continue;
                        const vid_t cu = comp[u];
                        if (cu != cv) break;
                    }
                } else {
#pragma unroll(2)
                    for (; newEdgesStart < edgesEnd; ++newEdgesStart) {
                        const vid_t u = edges[newEdgesStart].dest;
                        const vid_t cu = comp[u];
                        if (cu != cv) break;
                    }
                }

                if (newEdgesStart < edgesEnd) {
                    const weight_t weight = edges[newEdgesStart].weight;
                    if (weight < result[cv].weight) {
                        const vid_t u = edges[newEdgesStart].dest;
                        const vid_t cu = comp[u];
#ifdef USE_RESULT_VERTEX
                        result[cv] = Result{edges[newEdgesStart].weight, cu, v};
#else
                        result[cv] = Result{edges[newEdgesStart].weight, cu, newEdgesStart};
#endif
                    }
                }

                if (newEdgesStart != edgesStart) startEdgesIds[v] = newEdgesStart;
            }
            haveOuterComps[threadId] = outerComps;
            compBound[threadId][0] = minv;
            compBound[threadId][1] = maxv;
        } /* definedIter */
        times[iterationNumber][threadId][0] = rdtsc.end(threadId);
#pragma omp barrier

        //
        // reduce min
        //
        rdtsc.start(threadId);

#ifdef USE_FAST_REDUCTION
        bool doFastReduction = true;
        for (int i = 0; i < threadsCount; ++i) if (haveOuterComps[i]) {
            doFastReduction = false;
            break;
        }
#else
        const bool doFastReduction = false;
#endif

        if (definedIter == 0) {
            if (!threadId) updated = 1;

        } else if (doFastReduction) {
            int localUpdated = 0;
            const vid_t vto = vertexIds[threadId + 1];
            for (vid_t v = vertexIds[threadId]; v < vto; ++v) {
                if (comp[v] == v && localResult[threadId][v].weight <= MAX_WEIGHT) {
                    bestResult[v] = localResult[threadId][v];
                    //if (localResult[threadId][v].weight <= MAX_WEIGHT) {
                        comp[v] = localResult[threadId][v].destComp;
                        //isCoolEdge[localResult[threadId][v].edgeId] = true;
                        localUpdated = 1;
                        localResult[threadId][v].weight = MAX_WEIGHT + 0.1;
                    //}
                } else {
                    //assert(localResult[threadId][v].weight > MAX_WEIGHT);
                    bestResult[v].weight = MAX_WEIGHT + 0.1;
                }
            }

            if (localUpdated)
                updated = localUpdated;

        } else {

#if defined(USE_REDUCTION_SIMPLE)
#ifdef USE_COMPRESS
            int localUpdated = 0;
            const vid_t reduceInterval = lastUsedVid - prevUsedVid;
            const vid_t reduceFrom = prevUsedVid + int64_t(reduceInterval) * (threadId + 0) / threadsCount;
            const vid_t reduceTo   = prevUsedVid + int64_t(reduceInterval) * (threadId + 1) / threadsCount;
//#pragma omp for nowait
#if 1
            for (vid_t i = reduceFrom; i < reduceTo; ++i) {
		    thread_vector_t haveResult = 0;
		    for (int j = 0; j < threadsCount; ++j)
			    if (localResult[j][i].weight <= MAX_WEIGHT)
				    haveResult |= bit<thread_vector_t>(j);

		    int bestThread = -1;
		    for (int j = 0; j < threadsCount; ++j) if (haveResult & bit<thread_vector_t>(j)) {
			    if (bestThread == -1 || localResult[j][i].weight < localResult[bestThread][i].weight)
				    bestThread = j;
		    }

		    if (bestThread == -1) {
			    bestResult[i].weight = MAX_WEIGHT + 0.1;
		    } else {
			    bestResult[i] = localResult[bestThread][i];
			    comp[i] = localResult[bestThread][i].destComp;
			    //isCoolEdge[localResult[bestThread][i].edgeId] = true;
			    localUpdated = 1;
		    }

		    //#pragma ivdep
		    for (int j = 0; j < threadsCount; ++j) if (haveResult & bit<thread_vector_t>(j)) {
			    localResult[j][i].weight = MAX_WEIGHT + 0.1;
		    }
            }
#else
	    //memcpy(bestResult + reduceFrom, localResult[threadId] + reduceFrom, sizeof(Result) * (reduceTo - reduceFrom));
            for (vid_t i = reduceFrom; i < reduceTo; ++i) {
                bestResult[i] = localResult[threadId][i];
                if (localResult[threadId][i].weight <= MAX_WEIGHT)
                    localResult[threadId][i].weight = MAX_WEIGHT + 0.1;
            }

            for (int t = 0; t < threadsCount; ++t) if (t != threadId) {
                for (vid_t i = reduceFrom; i < reduceTo; ++i) {
                    if (localResult[t][i].weight < bestResult[i].weight) {
                        bestResult[i] = localResult[t][i];
                    }
                    localResult[t][i].weight = MAX_WEIGHT + 0.1;
                }
            }

            for (vid_t i = reduceFrom; i < reduceTo; ++i) if (bestResult[i].weight <= MAX_WEIGHT) {
                comp[i] = bestResult[i].destComp;
                localUpdated = 1;
            }
#endif

            if (localUpdated) updated = 1;

#else
#pragma omp for reduction(+:updated) nowait
            for (vid_t i = 0; i < vertexCount; ++i) {
                if (comp[i] == i) {
                    thread_vector_t haveResult = 0;
                    for (int j = 0; j < threadsCount; ++j)
                        if (localResult[j][i].weight <= MAX_WEIGHT)
                            haveResult |= bit<thread_vector_t>(j);

                    int bestThread = -1;
                    for (int j = 0; j < threadsCount; ++j) if (haveResult & bit<thread_vector_t>(j)) {
                        if (bestThread == -1 || localResult[j][i].weight < localResult[bestThread][i].weight)
                            bestThread = j;
                    }

                    if (bestThread == -1) {
                        bestResult[i].weight = MAX_WEIGHT + 0.1;
                    } else {
                        bestResult[i] = localResult[bestThread][i];
                        comp[i] = localResult[bestThread][i].destComp;
                        //isCoolEdge[localResult[bestThread][i].edgeId] = true;
                        updated = 1;
                    }

                    for (int j = 0; j < threadsCount; ++j) if (haveResult & bit<thread_vector_t>(j)) {
                        localResult[j][i].weight = MAX_WEIGHT + 0.1;
                    }
                } else {
                    //bestResult[i].weight = MAX_WEIGHT + 0.1;
                }
            }
#endif /* USE_COMPRESS */
#elif defined(USE_REDUCTION_TREE)
            int localUpdated = 0;
            const int threadsPerSocket = threadsCount / 2;
            const int socketId = threadId / threadsPerSocket;
            const int threadsFrom = threadsPerSocket * socketId;
            const int threadsTo = threadsPerSocket * (socketId + 1);

            const vid_t reduceInterval = lastUsedVid - prevUsedVid;

            for (int reductionIter = 0; reductionIter < 2; ++reductionIter) {
                const int redThreadId = (reductionIter == 0 ? threadId : (threadId + threadsPerSocket) % threadsCount);
                const vid_t reduceFrom = prevUsedVid + int64_t(reduceInterval) * (redThreadId + 0) / threadsCount;
                const vid_t reduceTo   = prevUsedVid + int64_t(reduceInterval) * (redThreadId + 1) / threadsCount;

                for (int t = threadsFrom; t < threadsTo; ++t) if (t != threadId) {
                    const vid_t tReduceFrom = std::max(reduceFrom, compBound[t][0]);
                    const vid_t tReduceTo   = std::min(reduceTo,   compBound[t][1] + 1);
                    for (vid_t i = tReduceFrom; i < tReduceTo; ++i) {
                        if (localResult[t][i].weight < localResult[threadId][i].weight) {
                            localResult[threadId][i] = localResult[t][i];
                        }
                        localResult[t][i].weight = MAX_WEIGHT + 0.1;
                    }
                }
            }
#pragma omp barrier
            {
                const vid_t reduceFrom = prevUsedVid + int64_t(reduceInterval) * (threadId + 0) / threadsCount;
                const vid_t reduceTo   = prevUsedVid + int64_t(reduceInterval) * (threadId + 1) / threadsCount;
                const int opponentThread = (threadId + threadsPerSocket) % threadsCount;
                for (vid_t v = reduceFrom; v < reduceTo; ++v) {
                    if (localResult[threadId][v].weight < localResult[opponentThread][v].weight)
                        bestResult[v] = localResult[threadId][v];
                    else
                        bestResult[v] = localResult[opponentThread][v];
                    localResult[threadId][v].weight = localResult[opponentThread][v].weight = MAX_WEIGHT + 0.1;
                    if (bestResult[v].weight <= MAX_WEIGHT) {
                        comp[v] = bestResult[v].destComp;
                        localUpdated = 1;
                    }
                }
            }
            if (localUpdated) updated = 1;

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
#else /* REDUCTION_TYPE */
#error reduction type should be set
#endif /* REDUCTION_TYPE */
        }
        times[iterationNumber][threadId][1] = rdtsc.end(threadId);

#pragma omp barrier

        //
        // merge components
        //
        eid_t *mergeData = Answer::data[threadId];
        vid_t mergePos = Answer::pos[threadId][0];
        rdtsc.start(threadId);
        vid_t localGeneratedComps = 0;
        if (!threadId) Eo(lastUsedVid - prevUsedVid);
#ifdef ON_HOME
#pragma omp for reduction(+:tmpTaskResult) nowait
#else
#pragma omp for nowait
#endif /* ON_HOME */
#ifdef USE_COMPRESS
        for (vid_t i = prevUsedVid; i < lastUsedVid; ++i) {
#else
        for (vid_t i = 0; i < vertexCount; ++i) {
#endif /* USE_COMPRESS */
            if (definedIter == 0) {
                const vid_t oc = comp[i];
                if (oc == i) {
                    //bestResult[i].weight = MAX_WEIGHT + 0.1;
                } else {
                    if (comp[oc] == i && i < oc) {
                        comp[i] = i;
#ifdef USE_COMPRESS
                        ++localGeneratedComps;
#endif /* USE_COMPRESS */
                    } else {
#ifdef USE_RESULT_VERTEX
                        const eid_t edgeId = startEdgesIds[i] - 1;
                        //isCoolEdge[edgeId] = true;
                        mergeData[mergePos++] = edgeId;
#else
                        isCoolEdge[best.edgeId] = true;
#endif /* USE_RESULT_VERTEX */
#ifdef ON_HOME
                        tmpTaskResult += edges[edgeId].weight;
#endif /* ON_HOME */
                    }
                }
            } else {
                Result& best = bestResult[i];
                if (best.weight > MAX_WEIGHT) continue;
                vid_t oc = best.destComp;

                if (comp[oc] == i && i < oc) {
                    comp[i] = i;
#ifdef USE_COMPRESS
                    ++localGeneratedComps;
#endif /* USE_COMPRESS */
                } else {
#ifdef ON_HOME
                    tmpTaskResult += best.weight;
#endif /* ON_HOME */
#ifdef USE_RESULT_VERTEX
                    const eid_t edgeId = startEdgesIds[best.from];
                    //isCoolEdge[edgeId] = true;
                    mergeData[mergePos++] = edgeId;
#else
                    isCoolEdge[best.edgeId] = true;
#endif /* USE_RESULT_VERTEX */
                }
                bestResult[i].weight = MAX_DROPPED_WEIGHT;
            }
        }
#ifdef USE_SKIP_LAST_ITER
        diedComponents = localGeneratedComps;
#endif /* USE_SKIP_LAST_ITER */
        Answer::pos[threadId][0] = mergePos;
        times[iterationNumber][threadId][2] = rdtsc.end(threadId);

#ifdef USE_COMPRESS
        generatedComps[threadId] = localGeneratedComps;
#pragma omp barrier
#pragma omp single
        {
            for (int i = 1; i < threadsCount; ++i)
                generatedComps[i] += generatedComps[i - 1];
        }
#pragma omp barrier

        //
        // compress component data
        //
        vid_t currentVid = lastUsedVid + (threadId ? generatedComps[threadId - 1] : 0);
#pragma omp for nowait
        for (vid_t i = prevUsedVid; i < lastUsedVid; ++i) if (comp[i] == i && bestResult[i].weight == MAX_DROPPED_WEIGHT)
            comp[i] = currentVid++;
        assert(currentVid == lastUsedVid + generatedComps[threadId]);
#endif /* USE_COMPRESS */
    }
#ifdef USE_COMPRESS
    prevUsedVid = lastUsedVid;
    lastUsedVid += generatedComps[threadsCount - 1];
#endif /* USE_COMPRESS */

#ifdef USE_SKIP_LAST_ITER
    aliveComponents = diedComponents;
    if (aliveComponents == 1) {
        updated = 0;
    }
#endif /* USE_SKIP_LAST_ITER */

    //
    // pointer jumping
    //
#ifdef USE_COMPRESS
    int changed = 100500;
    if (definedIter == 0) {
        do {
            changed = 0;
#pragma omp parallel 
            {
                const int threadId = omp_get_thread_num();
                rdtsc.start(threadId);
#pragma omp for reduction(+:changed) nowait
                for (vid_t i = 0; i < prevUsedVid; ++i) {
                    const vid_t myComp = comp[i];
                    __builtin_prefetch(comp + comp[i + PREFETCH_PJ_COMP]);
                    if (myComp == i) continue;
                    const vid_t parentComp = comp[myComp];
                    if (myComp != parentComp) {
                        comp[i] = parentComp;
                        changed = 1;
                    }
                }
                times[iterationNumber][threadId][3] += rdtsc.end(threadId);
            }
        } while (changed);
        
    } else {
        while (changed) {
            changed = 0;
#pragma omp parallel 
            {
                const int threadId = omp_get_thread_num();
                rdtsc.start(threadId);
#pragma omp for reduction(+:changed) nowait
                for (vid_t i = vertexCount; i < prevUsedVid; ++i) {
                    const vid_t myComp = comp[i];
                    if (myComp == i) continue;
                    const vid_t parentComp = comp[myComp];
                    if (myComp != parentComp) {
                        comp[i] = parentComp;
                        changed = 1;
                    }
                }
                times[iterationNumber][threadId][3] += rdtsc.end(threadId);
            }
        }

#pragma omp parallel 
        {
            const int threadId = omp_get_thread_num();
            rdtsc.start(threadId);
#pragma omp for nowait
            for (vid_t i = 0; i < vertexCount; ++i) {
                const vid_t myComp = comp[i];
                if (myComp == i) continue;
                const vid_t parentComp = comp[myComp];
                if (myComp != parentComp) {
                    comp[i] = parentComp;
                }
            }
            times[iterationNumber][threadId][3] += rdtsc.end(threadId);
        }
    }

#else /* NOT USE_COMPRESS */
    int changed = 100500;
    while (changed) {
        changed = 0;
#pragma omp parallel 
        {
            const int threadId = omp_get_thread_num();
            rdtsc.start(threadId);
            stickThisThreadToCore(threadId);
#pragma omp for reduction(+:changed) nowait
            for (vid_t i = 0; i < vertexCount; ++i) {
                const vid_t myComp = comp[i];
                if (myComp == i) continue;
                const vid_t parentComp = comp[myComp];
                if (myComp != parentComp) {
                    comp[i] = parentComp;
                    changed = 1;
                }
            }
            times[iterationNumber][threadId][3] += rdtsc.end(threadId);
        }
    }
#endif /* USE_COMPRESS */

#ifdef USE_SKIP_LOOPS
    if (loopsViewEdges) {
        double skipRatio = double(loopsSkipEdges) / loopsViewEdges;
        E(loopsViewEdges);
        E(loopsSkipEdges);
        Eo(skipRatio);
        loopsViewEdges = loopsSkipEdges = 0;
        if (skipRatio < 0.7) useSkipLoops = false;
    }
#endif

    taskResult += tmpTaskResult;
    ++iterationNumber;
    return updated;
}

void doReset() {
    iterationNumber = 0;
    taskResult = 0;
    useMagic = maxWindowSize < vertexCount / 2;
#ifdef USE_ANSWER_IN_VECTOR
    Answer::reset();
#endif

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);

#if 0
#pragma omp for nowait
    for (eid_t e = 0; e < edgesCount; ++e)
        isCoolEdge[e] = false;
#elif 0
    eid_t startOffset = int64_t(edgesCount) * (threadId + 0) / threadsCount;
    eid_t endOffset   = int64_t(edgesCount) * (threadId + 1) / threadsCount;
    memset(isCoolEdge + startOffset, 0, endOffset - startOffset);
#endif

#ifdef USE_COMPRESS
    /*
#pragma omp for nowait
    for (vid_t i = 0; i < vertexCount; ++i) {
        comp[i] = i;
        //bestResult[i].weight = 0;
    }
    */
#pragma omp for nowait
    for (vid_t i = vertexCount; i < lastUsedVid; ++i) { // TODO fix len
        comp[i] = i;
        bestResult[i].weight = 0; // TODO remove?
    }
#else
#pragma omp for nowait
    for (vid_t i = 0; i < vertexCount; ++i) {
        comp[i] = i;
        //bestResult[i].weight = 0;
    }
#endif /* USE_COMPRESS */

        // TODO fix
#ifdef USE_COMPRESS
#pragma omp for nowait
        for (vid_t i = vertexCount; i < lastUsedVid; ++i) // TODO decrase to one?
            localResult[threadId][i] = Result{MAX_WEIGHT + 0.1, 0, 0};
#else
#pragma omp for nowait
        for (vid_t i = vertexCount; i < vertexCount; ++i)
            localResult[threadId][i] = Result{MAX_WEIGHT + 0.1, 0, 0};
#endif /* USE_COMPRESS */

#if 0 // TODO copy and skip loops
        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            const eid_t startEdge = edgesIds[i];
            const eid_t endEdge = edgesIds[i + 1];
            startEdgesIds[i] = edgesIds[i];
            while (startEdgesIds[i] < endEdge && edges[startEdgesIds[i]].dest == i) ++startEdgesIds[i];
        }
#else
        memcpy(startEdgesIds + vertexIds[threadId], edgesIds + vertexIds[threadId], sizeof(eid_t) * (vertexIds[threadId + 1] - vertexIds[threadId]));
#endif
    }

#ifdef USE_COMPRESS
    lastUsedVid = vertexCount;
    prevUsedVid = 0;
#endif /* USE_COMPRESS */
#ifdef USE_SKIP_LOOPS
    useSkipLoops = true;
#endif
}

void doPrepare() {
#pragma omp parallel for
    for (vid_t v = 0; v < vertexCount; ++v)
        sort(edges + edgesIds[v], edges + edgesIds[v + 1], [&](const Edge& a, const Edge& b) -> bool {
                vid_t va = a.dest;
                eid_t da = edgesIds[va + 1] - edgesIds[va];
                vid_t vb = b.dest;
                eid_t db = edgesIds[vb + 1] - edgesIds[vb];
                //return a.weight < b.weight;
                return da < db;
            });
    doReorder();
#ifdef USE_SKIP_LAST_ITER
    aliveComponents = vertexCount;
#endif /* USE_SKIP_LAST_ITER */
#ifdef USE_ANSWER_IN_VECTOR
    Answer::init();
#endif

    vertexIds = new vid_t[threadsCount + 1]; // TODO threads & align
    vertexIds[0] = 0;
    vertexBound = static_cast<pvv*>(malloc(sizeof(pvv) * vertexCount));
    isCoolEdge = static_cast<bool*>(malloc(edgesCount));// new bool[edgesCount](); // TODO NUMA
#pragma omp parallel
    {
        stickThisThreadToCore(omp_get_thread_num());
#pragma omp for
    for (eid_t i = 0; i < edgesCount; ++i) isCoolEdge[i] = false;
    }
    //memset(isCoolEdge, 0, edgesCount);

#ifdef USE_COMPRESS
    // bestResult = new Result[vertexCount * 2]; ALLOC
    bestResult = static_cast<Result*>(malloc(sizeof(Result) * vertexCount * 2));
    // comp = new vid_t[vertexCount * 2]; ALLOC
    comp = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount * 2));
#pragma omp parallel
    {
stickThisThreadToCore(omp_get_thread_num());
#pragma omp for
    for (vid_t i = 0; i < vertexCount; ++i) {
        comp[i] = i;
        //bestResult[i].weight = 0;
    }
    }
#pragma omp parallel for
    for (vid_t i = vertexCount; i < vertexCount * 2; ++i) {
        comp[i] = i;
        //bestResult[i].weight = 0;
    }

#else
    bestResult = new Result[vertexCount];
    comp = new vid_t[vertexCount];
    for (vid_t i = 0; i < vertexCount; ++i)
        comp[i] = i;
#endif /* USE_COMPRESS */

    // startEdgesIds = new eid_t[vertexCount]; ALLOC
    startEdgesIds = static_cast<eid_t*>(malloc(sizeof(eid_t) * (vertexCount + PREFETCH_START_EDGE)));
    for (vid_t i = vertexCount; i < vertexCount + PREFETCH_START_EDGE; ++i) { // TODO allocate on last thread
	    startEdgesIds[i] = edgesCount-1;
    }
#ifdef USE_COMPRESS
    //generatedComps = new vid_t[threadsCount];
    generatedComps = static_cast<int64_t*>(malloc(sizeof(int64_t) * threadsCount));
    lastUsedVid = vertexCount;
    prevUsedVid = 0;
#endif /* USE_COMPRESS */

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        stickThisThreadToCore(threadId);
        const int curThreadsCount = omp_get_num_threads();
        if (!threadId) {
            threadsCount = curThreadsCount;
            localResult = new Result*[threadsCount];
        }
#pragma omp barrier

#if 0
        const int64_t paramA = 10;
        const int64_t paramB = 1;
        const int64_t sumAll = paramA * vertexCount + paramB * edgesCount;
        const int64_t expectedSum = sumAll * (threadId + 1) / threadsCount;
        //eid_t degreeEnd = int64_t(edgesCount) * (threadId + 1) / threadsCount;
        //eid_t degreeSum = 0;
        int64_t currentSum = 0;
        for (vid_t i = 0; i < vertexCount; ++i) {
            const eid_t diff = edgesIds[i + 1] - edgesIds[i];
            /*
            degreeSum += diff;
            if (degreeSum >= degreeEnd) {
                vertexIds[threadId + 1] = i + 1;
                break;
            }
            */
            currentSum += paramA + paramB * diff;
            if (currentSum >= expectedSum) {
                vertexIds[threadId + 1] = i + 1;
                break;
            }
        }
#pragma omp critical
        {
            E(threadId); Eo(vertexIds[threadId + 1]);
        }
#elif 1
        vertexIds[threadId + 1] = int64_t(vertexCount) * (threadId + 1) / threadsCount;
#else
        const eid_t degreeEnd = int64_t(edgesCount) * (threadId + 1) / threadsCount;
        eid_t degreeSum = 0;
        vertexIds[threadId + 1] = -1;
        for (vid_t i = 0; i < vertexCount; ++i) {
            eid_t diff = edgesIds[i + 1] - edgesIds[i];
            degreeSum += diff;
            if (degreeSum >= degreeEnd && componentEnd[i]) {
                vertexIds[threadId + 1] = i + 1;
                break;
            }
        }
        assert(vertexIds[threadId + 1] > 0);
#endif /* vertexes distribution */

#pragma omp barrier

#ifdef USE_COMPRESS
        localResult[threadId] = new Result[vertexCount * 2];
        for (vid_t i = 0; i < vertexCount * 2; ++i)
            localResult[threadId][i] = Result{MAX_WEIGHT + 0.1, 0, 0};
#else
        localResult[threadId] = new Result[vertexCount];
        for (vid_t i = 0; i < vertexCount; ++i)
            localResult[threadId][i] = Result{MAX_WEIGHT + 0.1, 0, 0};
#endif /* USE_COMPRESS */

        for (vid_t i = vertexIds[threadId]; i < vertexIds[threadId + 1]; ++i) {
            const eid_t startEdge = edgesIds[i];
            const eid_t endEdge = edgesIds[i + 1];
            vid_t minv = vertexCount, maxv = 0;
            for (eid_t e = startEdge; e < endEdge; ++e) {
                minv = std::min(minv, edges[e].dest);
                maxv = std::max(maxv, edges[e].dest);
            }
            vertexBound[i] = pvv(minv, maxv);
            maxWindowSize = std::max(maxWindowSize, maxv - minv);
            sort(edges + startEdge, edges + endEdge, EdgeWeightCmp());
            startEdgesIds[i] = edgesIds[i];
            while (startEdgesIds[i] < edgesIds[i + 1] && edges[startEdgesIds[i]].dest == i) ++startEdgesIds[i];
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

#ifndef ON_DISLAB
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

    int64_t calcTime;
    for (int i = 0; i < 1; ++i) {
        memset(times, 0, sizeof(times));
        calcTime = -currentNanoTime();
        doReset();
        if (doAll<0>()) {
            while (doAll<1>());
        }
        //while (doAll());
        calcTime += currentNanoTime();
    }

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


    diedComps.print("Died comps:", "%8d ");
    potentComps.print("Ptnt comps:", "%8d ");
    catchComps.print("Cthd comps:", "%8d ");

    return 0;
}
#else
void init_mst(graph_t *G) {   
#ifndef USE_HYPERTHREADING
    omp_set_num_threads(16);
#else
    omp_set_num_threads(32);
#endif
    convertAll(G);
    warmup();
    doPrepare();
    doReset();
}   

void* MST(graph_t *) {
    memset(times, 0, sizeof(times));
    int64_t prepareTime = -currentNanoTime();
    doReset();
    prepareTime += currentNanoTime();
    if (doAll<0>())
        while (doAll<1>());
    //while (doAll());
    //fprintf(stderr, "prepare time is %.7f\n", double(prepareTime) / 1e9);
    return NULL;
}

void convert_to_output(graph_t *G, void* , forest_t *trees_output) {
#ifdef USE_ANSWER_IN_VECTOR
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        for (vid_t i = 0; i < Answer::pos[threadId][0]; ++i)
            isCoolEdge[Answer::data[threadId][i]] = true;
    }
#endif /* USE_ANSWER_IN_VECTOR */
    map<vid_t, vector<eid_t>> treesInMap;
    for (vid_t i = 0; i < vertexCount; ++i) if (comp[i] == i) {
        treesInMap[i] = vector<eid_t>();
    }
    double sum = 0;
    eid_t lower = 0, higher = 0;
    for (vid_t v = 0; v < vertexCount; ++v) {
        const eid_t startEdge = edgesIds[v];
        const eid_t endEdge = edgesIds[v + 1];
        for (eid_t e = startEdge; e < endEdge; ++e) {
            assert(0 <= e);
            assert(e < edgesCount);
            if (isCoolEdge[e]) {
                vid_t to = comp[edges[e].dest];
                sum += edges[e].weight;
                //treesInMap[to].push_back(startEdge + edges[e].origOffset);
                treesInMap[to].push_back(G->rowsIndices[que[v]] + edges[e].origOffset);
                if (edges[e].weight < 0.1)
                    ++lower;
                else
                    ++higher;
            }
        }
    }
    //Eo(lower); Eo(higher);
    trees_output->numTrees = treesInMap.size();
    trees_output->numEdges = vertexCount - trees_output->numTrees;
    trees_output->p_edge_list = static_cast<edge_id_t*>(malloc(sizeof(edge_id_t) * 2 * trees_output->numTrees));
    trees_output->edge_id     = static_cast<edge_id_t*>(malloc(sizeof(edge_id_t) * trees_output->numEdges));
    vid_t currentTree = 0;
    vid_t currentEdge = 0;
    for (const auto& tree : treesInMap) {
        trees_output->p_edge_list[currentTree++] = currentEdge;
        for (eid_t eid : tree.second)
            trees_output->edge_id[currentEdge++] = eid;
        trees_output->p_edge_list[currentTree++] = currentEdge;
    }
}

void finalize_mst(graph_t *) {
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
}
#endif
