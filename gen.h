#pragma once

#include <cinttypes>
#include <utility>

typedef int32_t vid_t;
typedef int32_t eid_t;
typedef double weight_t;
typedef std::pair<vid_t, vid_t> pvv;

struct Edge {
    vid_t dest;
    weight_t weight;
};

struct EdgeDestCmp { bool operator()(const Edge& a, const Edge& b) const; };
struct EdgeWeightCmp { bool operator()(const Edge& a, const Edge& b) const; };

extern vid_t vertexCount;
extern eid_t edgesCount;
extern eid_t *edgesIds;
extern Edge *edges;

extern const weight_t MAX_WEIGHT;

void readAll(char *filename);

int64_t currentNanoTime();
