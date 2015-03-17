#include "gen.h"
#include <cstdio>
#include <algorithm>
#include <tuple>
#include <vector>
#include <numeric>
using namespace std;

struct BiEdge {
    vid_t from, to;
    weight_t w;

    bool operator<(const BiEdge& o) const {
        return w < o.w;
    }
} *biedges;
eid_t biEdgesCount;

struct  {
    vid_t *cmp;

    void init(vid_t n) {
        cmp = static_cast<vid_t*>(malloc(sizeof(vid_t) * n));
        iota(cmp, cmp + n, 0);
    }

    vid_t get(vid_t v) {
        return cmp[v] == v ? v : cmp[v] = get(cmp[v]);
    }

    void merge(vid_t v, vid_t u) {
        //if (rand()&1) swap(v, u);
        cmp[u] = v;
    }

} snm;

struct Result {
    vid_t to;
    weight_t weight;
} *result;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }

    int64_t timeRead = -currentNanoTime();
    readAll(argv[1]);
    timeRead += currentNanoTime();

    int64_t timeInit = -currentNanoTime();
    biedges = static_cast<BiEdge*>(malloc(sizeof(BiEdge) * edgesCount));
    for (vid_t v = 0; v < vertexCount; ++v) {
        for (eid_t i = edgesIds[v]; i < edgesIds[v + 1]; ++i) {
            const vid_t u = edges[i].dest;
            if (u <= v) continue;
            const weight_t w = edges[i].weight;
            biedges[biEdgesCount++] = BiEdge{v, u, w};
        }
    }
    snm.init(vertexCount);
    result = static_cast<Result*>(malloc(sizeof(Result) * vertexCount));
    timeInit += currentNanoTime();

    int64_t timeBuild = -currentNanoTime();
    weight_t wresult = 0.0;
    bool changed = 0;
    do {
        changed = 0;
        for (vid_t v = 0; v < vertexCount; ++v)
            result[v].weight = MAX_WEIGHT + 0.1;
        for (eid_t e = 0; e < biEdgesCount; ++e) {
            const vid_t v = biedges[e].from;
            const vid_t cv = snm.get(v);
            const vid_t u = biedges[e].to;
            const vid_t cu = snm.get(u);
            if (cv != cu) {
                //E(cv); Eo(cu);
                changed = true;
                const weight_t w = biedges[e].w;
                if (w < result[cv].weight) {
                    result[cv].weight = w;
                    result[cv].to = cu;
                }
                if (w < result[cu].weight) {
                    result[cu].weight = w;
                    result[cu].to = cv;
                }
            }
        }
        for (vid_t v = 0; v < vertexCount; ++v) if (snm.get(v) == v && result[v].weight <= MAX_WEIGHT) {
            wresult += result[v].weight;
            snm.merge(v, snm.get(result[v].to));
        }
    } while (changed);
    timeBuild += currentNanoTime();
    
    printf("%.10lf\n", double(wresult));
    fprintf(stderr, "Read time: %.5lf\n", double(timeRead) / 1e9);
    fprintf(stderr, "Init time: %.5lf\n", double(timeInit) / 1e9);
    fprintf(stderr, "Bild time: %.5lf\n", double(timeBuild) / 1e9);
    fprintf(stderr, "\n");
    fprintf(stderr, "%.3lf\n%.3lf\n", double(timeRead + timeInit) / 1e9, double(timeBuild) / 1e9);

    return 0;
}
