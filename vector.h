#include <cinttypes>
#include <cstdlib>
#include <cstring>

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

    int64_t size() {
        return vectorSize;
    }

    bool empty() {
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
};

#ifdef TEST
struct TestVector {
    TestVector();
};
extern TestVector testVector;
#endif
