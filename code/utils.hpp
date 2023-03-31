#pragma once

#define rep(i, a, b) for (auto i = (a); i < (b); ++i)

#include <algorithm>
#include <vector>

#include "graph.hpp"

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

template <class T>
void flatten_3d_vector(
    const std::vector<std::vector<std::vector<T>>>& to_flatten,
    std::vector<T>& flattened, std::vector<int>& cnt, std::vector<int>& offsets,
    std::vector<std::vector<int>>& cnt2, std::vector<std::vector<int>>& off2) {
    // first dimenion is the number of mpi ranks
    // second dimension is the number of threads on this node

    const auto n = to_flatten.size();
    const auto m = to_flatten[0].size();

    cnt.assign(n, 0);
    offsets.resize(n + 1);
    offsets[0] = 0;

    cnt2.resize(n);
    off2.resize(n);

    rep(i, 0u, n) {
        cnt2[i].resize(m);

        off2[i].resize(m + 1);
        off2[i][0] = 0;
        rep(j, 0u, m) {
            cnt2[i][j] = to_flatten[i][j].size();
            cnt[i] += cnt2[i][j];
            off2[i][j + 1] = off2[i][j] + cnt2[i][j];
        }

        offsets[i + 1] = offsets[i] + cnt[i];
    }

    flattened.resize(offsets[n]);
    rep(i, 0u, n) {
#pragma omp parallel for schedule(dynamic, 1)
        rep(j, 0u, m) {
            std::copy_n(to_flatten[i][j].begin(), cnt2[i][j],
                        flattened.begin() + offsets[i] + off2[i][j]);
        }
    }
}

template <class T>
void flatten_3d_vector(
    const std::vector<std::vector<std::vector<T>>>& to_flatten,
    std::vector<T>& flattened, std::vector<int>& cnt,
    std::vector<int>& offsets) {
    // first dimenion is the number of mpi ranks
    // second dimension is the number of threads on this node

    const auto n = to_flatten.size();
    const auto m = to_flatten[0].size();

    std::vector off2(n, std::vector(m + 1, 0));

    cnt.assign(n, 0);
    offsets.resize(n + 1);
    offsets[0] = 0;

    rep(i, 0u, n) {
        rep(j, 0u, m) {
            const int cur = (int)to_flatten[i][j].size();
            cnt[i] += cur;
            off2[i][j + 1] = off2[i][j] + cur;
        }
        offsets[i + 1] = offsets[i] + cnt[i];
    }

    flattened.resize(offsets[n]);
    rep(i, 0u, n) {
#pragma omp parallel for schedule(dynamic, 1)
        rep(j, 0u, m) {
            std::copy(to_flatten[i][j].begin(), to_flatten[i][j].end(),
                      flattened.begin() + offsets[i] + off2[i][j]);
        }
    }
}
