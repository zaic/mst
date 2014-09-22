#include "gen.h"
#include <thread>
#ifdef __clang__
#include "omp.h"
#else
#include <omp.h>
#endif
using namespace std;



void go(int threadNumber) {
    volatile int64_t fpre = 0, fcur = 1;
    stickThisThreadToCore(threadNumber);
    while (true) {
        int64_t fnext = fpre + fcur;
        fpre = fcur;
        fcur = fnext;
    }
}


int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s start_core  end_core\n", argv[0]);
        fprintf(stderr, "         [start_core; end_core)\n");
        return 1;
    }
    const int from = atoi(argv[1]);
    const int to = atoi(argv[2]);
    const int count = to - from;
    thread ts[to];
    for (int i = from; i < to; ++i) {
        ts[i] = thread(go, i);
    }
    ts[from].join();

    return 0;
}
