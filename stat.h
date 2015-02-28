#pragma once

#include "gen.h"

template<typename T>
struct Stat {
    T data[kMaxThreads][kMaxIterations];

    Stat() {
        for (int i = 0; i < kMaxThreads; ++i)
            for (int j = 0; j < kMaxIterations; ++j)
                data[i][j] = T();
    }

    void set(int threadId, const T& value) {
        data[threadId][iterationNumber] = value;
    }

    void add(int threadId, const T& value) {
        data[threadId][iterationNumber] += value;
    }

    void print(const char *msg, const char *fmt) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Stat: %s\n", msg);
        for (int i = 0; i < iterationNumber; ++i) {
            fprintf(stderr, "iteration %2d: ", i);
            for (int j = 0; j < threadsCount; ++j) {
                fprintf(stderr, fmt, data[j][i]);
            }
            fputs("\n", stderr);
        }
        fprintf(stderr, "\n");
    }
};
