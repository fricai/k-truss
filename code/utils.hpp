#pragma once

#define rep(i, a, b) for (auto i{a}; i < (b); ++i)

#include <algorithm>
#include <vector>

template <class T>
void clear_vector(std::vector<T>& v) {
    return std::vector<T>().swap(v);
}

template <class T>
void flatten_vector(const std::vector<std::vector<T>>& to_flatten,
                    std::vector<T>& flattened, std::vector<int>& cnt,
                    std::vector<int>& offsets) {
    const int n = (int)to_flatten.size();
    cnt.resize(n);
    offsets.resize(n + 1);

    rep(i, 0, n) cnt[i] = (int)to_flatten[i].size();
    rep(i, 0, n) offsets[i + 1] = offsets[i] + cnt[i];

    flattened.resize(offsets[n]);
    rep(i, 0, n) std::copy_n(to_flatten[i].begin(), cnt[i],
                             flattened.begin() + offsets[i]);
}

template <class T>
void destructive_flatten_vector(std::vector<std::vector<T>>& to_flatten,
                                std::vector<T>& flattened,
                                std::vector<int>& cnt,
                                std::vector<int>& offsets) {
    const int n = (int)to_flatten.size();
    cnt.resize(n);
    offsets.resize(n + 1);

    rep(i, 0, n) cnt[i] = (int)to_flatten[i].size();
    rep(i, 0, n) offsets[i + 1] = offsets[i] + cnt[i];

    flattened.resize(offsets[n]);
    rep(i, 0, n) {
        std::copy_n(to_flatten[i].begin(), cnt[i],
                    flattened.begin() + offsets[i]);
        clear_vector(to_flatten[i]);
    }
}
