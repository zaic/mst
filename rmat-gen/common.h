#pragma once

#include <iostream>

using std::cerr;
using std::endl;

#ifdef DEBUG
#   define E(x) { cerr << #x << " = " << (x) << "   "; }
#   define Eo(x) { cerr << #x << " = " << (x) << endl; }
#else
#   define E(x)
#   define Eo(x)
#endif
#define EO(x) Eo(x)

static int getConcurrencyLevel() {
#ifdef NUM_THREADS
    return NUM_THREADS;
#else
    static const int res = sysconf(_SC_NPROCESSORS_ONLN);
    return res;
#endif
}

#include <time.h>

#if defined(CLOCK_MONOTONIC)
#define CLOCK CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define CLOCK CLOCK_REALTIME
#else
#error "Failed to find a timing clock."
#endif

#include <pthread.h>
#include <unistd.h>

static int stick_this_thread_to_core(int core_id) {
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (core_id >= num_cores)
        return 0;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t current_thread = pthread_self();
    return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}
