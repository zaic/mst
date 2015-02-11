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
#include "vector.h"
using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input\n", argv[0]);
        return 1;
    }
    readAll(argv[1]);
    fprintf(stderr, "Done\n");

    puts("Graph G {");
    for (vid_t v = 0; v < vertexCount; ++v) {
        for (eid_t e = edgesIds[v]; e < edgesIds[v + 1]; ++e) {
            vid_t u = edges[e].dest;
            if (u > v)
                printf(" %d -- %d;\n", v, u);
        }
    }
    puts("}");

    return 0;
}
