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
        return tie(w, from, to) < tie(o.w, o.from, o.to);
    }
};

struct  {
    vector<vid_t> cmp;

    void init(vid_t n) {
        cmp = vector<vid_t>(n, 0);
        iota(cmp.begin(), cmp.end(), 0);
    }

    vid_t get(vid_t v) {
        return cmp[v] == v ? v : cmp[v] = get(cmp[v]);
    }

    void merge(vid_t v, vid_t u) {
        if (rand()&1) swap(v, u);
        cmp[v] = u;
    }

} snm;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }
    readAll(argv[1]);

    int64_t prepareTime = -currentNanoTime();
    vector<BiEdge> es;
    for (vid_t v = 0; v < vertexCount; ++v) {
        for (eid_t i = edgesIds[v]; i < edgesIds[v + 1]; ++i) {
            vid_t u = edges[i].dest;
            if (u <= v) continue;
            weight_t w = edges[i].weight;
            es.push_back(BiEdge{v, u, w});
        }
    }
    prepareTime += currentNanoTime();

    int64_t calcTime = -currentNanoTime();
    sort(es.begin(), es.end());

    weight_t result = 0.0;
    snm.init(vertexCount);
    for (const BiEdge& e : es) {
        int cv = snm.get(e.from);
        int cu = snm.get(e.to);
        if (cv != cu) {
            result += e.w;
            snm.merge(cv, cu);
        }
    }
    calcTime += currentNanoTime();
    
    printf("%.10lf\n", double(result));
    fprintf(stderr, "%.3lf\n%.3lf\n", double(prepareTime) / 1e9, double(calcTime) / 1e9);

    return 0;
}
