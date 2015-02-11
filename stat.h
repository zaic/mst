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
};
