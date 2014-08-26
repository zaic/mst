#include "gen.h"
#include <cassert>
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

vector<vid_t> comp;
vector<BiEdge> es;
weight_t result;

void doFindMin() {
    vector<weight_t> curMin(vertexCount, MAX_WEIGHT + 0.1);
    vector<pvv> curRes(vertexCount, pvv());
    for (const BiEdge& e : es) {
        vid_t cv = comp[e.from];
        vid_t cu = comp[e.to];
        auto relax = [&](vid_t c) {
            if (curMin[c] > e.w) {
                curMin[c] = e.w;
                curRes[c] = pvv(e.from, e.to);
            }
        };
        relax(cv);
        relax(cu);
    }

    for (vid_t i = 0; i < vertexCount; ++i) {
        if (curMin[i] > MAX_WEIGHT) continue;
        const pvv& t = curRes[i];
        vid_t otherComp = (comp[t.first] != i ? comp[t.first] : comp[t.second]);
        assert(otherComp != i);
        if (otherComp < i) {
            result += curMin[i];
            comp[i] = otherComp;
        }
    }
}

void doPointerJumping() {
    bool changed = true;
    while (changed) {
        changed = false;
        for (vid_t i = 0; i < vertexCount; ++i) {
            vid_t myComp = comp[i];
            vid_t parentComp = comp[myComp];
            if (myComp != parentComp) {
                comp[i] = parentComp;
                changed = true;
            }
        }
    }
}

void doReduceEdges() {
    vector<BiEdge> nextEdges;
    for (const BiEdge& e : es) {
        int compv = comp[e.from];
        int compu = comp[e.to];
        if (compv != compu)
            nextEdges.push_back(e);
    }
    es.swap(nextEdges);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }
    readAll(argv[1]);
    for (vid_t v = 0; v < vertexCount; ++v) {
        for (eid_t i = edgesIds[v]; i < edgesIds[v + 1]; ++i) {
            vid_t u = edges[i].dest;
            if (u <= v) continue;
            weight_t w = edges[i].weight;
            es.push_back(BiEdge{v, u, w});
        }
    }
    //sort(es.begin(), es.end());
    comp = vector<vid_t>(vertexCount, 0);
    iota(comp.begin(), comp.end(), 0);
    
    //doReduceEdges(); // remove self-loops
    while (!es.empty()) {
        doFindMin();
        doPointerJumping();
        doReduceEdges();
    }
    
    printf("%.10lf\n", double(result));

    return 0;
}
