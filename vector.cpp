#include "vector.h"
#include <cassert>

TestVector testVector;

TestVector::TestVector() {
    Vector<int, 10> v0;
    v0.reserve(15);
    for (int i = 0; i < 30; ++i)
        v0.pushBack(i);
    int sum = 29 * 30 / 2;
    for (int i = 0; i < v0.size(); ++i)
        sum -= v0.at(i);
    assert(sum == 0);

    Vector<int, 5> v1;
    for (int i = 0; i < 1000; ++i)
        v1.pushBack(i);
    sum = 999 * 1000 / 2;
    for (int i = 0; i < v1.size(); ++i)
        sum -= v1.at(i);
    assert(sum == 0);
}

int main() {
    return 0;
}
