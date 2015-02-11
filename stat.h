#pragma once

#include "gen.h"

template<typename T>
struct Stat {
    T data[kMaxThreads][kMaxIterations];

    Stat() {
        for (int i = 0; i < kMaxIterations; ++i)
            for (int j = 0; j < kMaxThreads; ++j)
                data[i][j] = T();
    }

    void set(int threadId, int iterationNumber, const T& value) {
        data[threadId][iterationNumber] = value;
    }

    void add(int threadId, int iterationNumber, const T& value) {
        data[threadId][iterationNumber] += value;
    }

    void print(int iterationNumber, int threadsCount, const char *msg, const char *fmt) {
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
