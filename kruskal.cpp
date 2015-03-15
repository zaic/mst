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
        cmp[v] = u;
    }

} snm;

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
    timeInit += currentNanoTime();

    int64_t timeSort = -currentNanoTime();
    sort(biedges, biedges + biEdgesCount);
    timeSort += currentNanoTime();

    int64_t timeBuild = -currentNanoTime();
    weight_t result = 0.0;
    snm.init(vertexCount);
    for (eid_t e = 0; e < biEdgesCount; ++e) {
        const int cv = snm.get(biedges[e].from);
        const int cu = snm.get(biedges[e].to);
        if (cv != cu) {
            result += biedges[e].w;
            snm.merge(cv, cu);
        }
    }
    timeBuild += currentNanoTime();
    
    printf("%.10lf\n", double(result));
    fprintf(stderr, "Read time: %.5lf\n", double(timeRead) / 1e6);
    fprintf(stderr, "Init time: %.5lf\n", double(timeInit) / 1e6);
    fprintf(stderr, "Sort time: %.5lf\n", double(timeSort) / 1e6);
    fprintf(stderr, "Bild time: %.5lf\n", double(timeBuild) / 1e6);

    return 0;
}
