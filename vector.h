#include <cinttypes>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include "gen.h"

template<typename T, int64_t N, size_t AllocationFactor = 2, size_t InitAllocationFactor = 1>
struct Vector {
    T onStack[N];

    T* onHeap;
    int64_t heapAllocated;
    int64_t vectorSize;



    Vector() {
        onHeap = NULL;
        heapAllocated = 0;
        vectorSize = 0;
    }

    ~Vector() {
        if (onHeap)
            free(onHeap);
    }

    void reserve(int64_t size) {
        size -= N;
        if (size <= 0) return ;
        onHeap = (T*)malloc(sizeof(T) * size);
        heapAllocated = size;
    }

    void clear() {
        if (onHeap) {
            free(onHeap);
            onHeap = NULL;
            heapAllocated = 0;
        }
        vectorSize = 0;
    }

    int64_t size() const {
        return vectorSize;
    }

    bool empty() const {
        return !vectorSize;
    }

    void pushBack(const T& value) {
        if (vectorSize < N) {
            onStack[vectorSize++] = value;
            return ;
        }
        int64_t curSize = vectorSize - N;
        if (curSize == heapAllocated) {
            if (!heapAllocated) {
                heapAllocated = N;
                onHeap = (T*)malloc(sizeof(T) * heapAllocated * InitAllocationFactor);
            } else {
                int64_t newSize = heapAllocated * AllocationFactor;
                T *newHeap = (T*)malloc(sizeof(T) * newSize);
                memcpy(newHeap, onHeap, sizeof(T) * curSize);
                free(onHeap);
                onHeap = newHeap;
                heapAllocated = newSize;
            }
        }
        onHeap[curSize] = value;
        vectorSize++;
    }

    void push_back(const T& value) {
        pushBack(value);
    }

    void removeAt(int64_t index) {
        at(index) = at(--vectorSize);
    }

    T& at(int64_t index) {
        if (index < N) {
            return onStack[index];
        } else {
            return onHeap[index - N];
        }
    }
    
    T& operator[](int64_t index) {
        return at(index);
    }

    void append(const Vector& other) {
        int64_t appendedSize = other.size();
    }
};

#ifdef TEST
struct TestVector {
    TestVector();
};
extern TestVector testVector;
#endif

template<typename T, typename index_t = size_t>
struct LargeVector {
    T *data;
    index_t size;

    index_t requestedSize;

    LargeVector() : data(NULL), size(0) { }
    ~LargeVector() {
        if (data)
            free(data);
    }

    void init(T *dataPtr) {
        size = 0;
        data = dataPtr; //(T*)(malloc(sizeof(T) * 1000));
    }

    void push_back(const T& value) {
        data[size++] = value;
    }

    void push_back(const LargeVector& other) {
        /*
        if (size + other.size > 1000) {
            T *rdata = (T*)malloc(sizeof(T) * requestedSize);
            memcpy(rdata, data, sizeof(T) * size);
            free(data);
            data = rdata;
        }
        */
        memcpy(data + size, other.data, sizeof(T) * other.size);
        size += other.size;
    }

    void swap(LargeVector& other) {
        std::swap(data, other.data);
        std::swap(size, other.size);
    }

    void removeAt(size_t index) {
        size--;
        if (index < size) std::swap(data[index], data[size]);
    }

    void clear() {
        /*
        if (data)
            free(data);
        size = 0;
        */
    }
};
