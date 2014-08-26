#include <cstdio>
#include <cassert>
#include <random>
#include <string>
#include "gen.h"
using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s n m filename\n", argv[0]);
        return 0;
    }
    vid_t n = static_cast<vid_t>(stoll(string(argv[1])));
    eid_t m = static_cast<eid_t>(stoll(string(argv[2])));

    random_device rd;
    uniform_int_distribution<vid_t> random_v(0, n - 1);
    uniform_real_distribution<weight_t> random_w(0, MAX_WEIGHT);

    vector<vector<Edge>> es(n);
    for (eid_t i = 0; i < m; ++i) {
        vid_t from = random_v(rd);
        vid_t to   = random_v(rd);
        weight_t w = random_w(rd);
        es[from].push_back(Edge{to, w});
        es[to].push_back(Edge{from, w});
    }
    m *= 2;

    edgesIds = new eid_t[n + 1];
    edges = new Edge[m];
    edgesIds[0] = 0;
    eid_t curEdge = 0;
    for (vid_t v = 0; v < n; ++v) {
        edgesIds[v + 1] = es[v].size() + edgesIds[v];
        for (auto& e : es[v]) {
            edges[curEdge++] = e;
        }
    }
    assert(curEdge == m);

    FILE *f = fopen(argv[3], "wb");
    fwrite(&n, sizeof(vid_t), 1, f);
    fwrite(&m, sizeof(eid_t), 1, f);
    fwrite(edgesIds, sizeof(eid_t), n + 1, f);
    fwrite(edges, sizeof(Edge), m, f);
    fclose(f);

    return 0;
}
