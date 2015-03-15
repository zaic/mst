#include "gen.h"
#include <cstdio>
#include <algorithm>
#include <tuple>
#include <vector>
#include <numeric>
#include <queue>
using namespace std;

vid_t *comp;

struct Que {
    weight_t w;
    vid_t to;

    bool operator<(const Que& other) const {
        return w > other.w;
    }
};

void addVertex(priority_queue<Que>& que, vid_t v) {
    for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e) {
        vid_t u = edges[e].dest;
        if (comp[u] == u)
            que.push(Que{edges[e].weight, u});
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }

    int64_t timeRead = -currentNanoTime();
    readAll(argv[1]);
    timeRead += currentNanoTime();

    int64_t timeInit = -currentNanoTime();
    comp = static_cast<vid_t*>(malloc(sizeof(vid_t) * vertexCount));
    iota(comp, comp + vertexCount, 0);
    timeInit += currentNanoTime();

    int64_t timeBuild = -currentNanoTime();
    weight_t result = 0;
    for (vid_t v = 0; v < vertexCount; ++v) if (comp[v] == v) {
        priority_queue<Que> que;
        addVertex(que, v);
        while (!que.empty()) {
            Que best = que.top();
            que.pop();
            if (comp[best.to] == v) continue;
            comp[best.to] = v;
            result += best.w;
            addVertex(que, best.to);
        }
    }
    timeBuild += currentNanoTime();
    
    printf("%.10lf\n", double(result));
    fprintf(stderr, "Read time: %.5lf\n", double(timeRead) / 1e9);
    fprintf(stderr, "Init time: %.5lf\n", double(timeInit) / 1e9);
    fprintf(stderr, "Bild time: %.5lf\n", double(timeBuild) / 1e9);
    fprintf(stderr, "\n");
    fprintf(stderr, "%.3lf\n%.3lf\n", double(timeRead + timeInit) / 1e9, double(0 + timeBuild) / 1e9);

    return 0;
}
