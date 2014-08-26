#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include "gen.h"
using namespace std;

random_device rd;
uniform_int_distribution<vid_t> random_v;
uniform_real_distribution<weight_t> random_w;

vector<int> dims, value;
vector<vector<Edge>> es;
vid_t n;
eid_t m;

void go(int dim) {
    if (dim == (int)dims.size()) {
        vid_t id = 0;
        vid_t mul = 1;
        vector<vid_t> ngbh;
        for (int i = 0; i < (int)dims.size(); ++i) {
            id += value[i] * mul;
            if (value[i])
                ngbh.push_back(-mul);
            if (value[i] < dims[i] - 1)
                ngbh.push_back(mul);
            mul *= dims[i];
        }
        vid_t v = id;
        for (vid_t i : ngbh) {
            vid_t u = id + i;
            assert(0 <= u && u < n);
            weight_t w = random_w(rd);
            es[v].push_back(Edge{u, w});
            es[u].push_back(Edge{v, w});
            m += 2;
        }

    } else {
        value.push_back(0);
        for (int i = 0; i < dims[dim]; ++i) {
            value.back() = i;
            go(dim + 1);
        }
        value.pop_back();
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s n [c1, c2, ..., cn] filename\n", argv[0]);
        return 1;
    }
    n = static_cast<vid_t>(stoll(string(argv[1])));
    vid_t oldn = n;
    if (argc != 3 + n) {
        fprintf(stderr, "Usage: %s n [c1, c2, ..., cn] filename\n", argv[0]);
        return 1;
    }
    vid_t tmpn = 1;
    for (int i = 0; i < n; ++i) {
        dims.push_back(atoi(argv[2 + i]));
        tmpn *= dims.back();
    }
    n = tmpn;

    random_v = uniform_int_distribution<vid_t>(0, n - 1);
    random_w = uniform_real_distribution<weight_t>(0, MAX_WEIGHT);
    es = vector<vector<Edge>>(n);

    go(0);

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

    FILE *f = fopen(argv[2 + oldn], "wb");
    fwrite(&n, sizeof(vid_t), 1, f);
    fwrite(&m, sizeof(eid_t), 1, f);
    fwrite(edgesIds, sizeof(eid_t), n + 1, f);
    fwrite(edges, sizeof(Edge), m, f);
    fclose(f);

    return 0;
}
