#pragma once

#define rep(i, a, b) for (auto i{a}; i < (b); ++i)

#include <vector>

template <class T>
void clear_vector(std::vector<T>& v) {
    return std::vector<T>().swap(v);
}
