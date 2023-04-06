#pragma once

#include <vector>

struct dsu {
    std::vector<int> par_or_size;
    explicit dsu(int n) : par_or_size(n, -1) {}

    int head(int u) {
        return par_or_size[u] < 0 ? u : par_or_size[u] = head(par_or_size[u]);
    }

    int size(int u) { return -par_or_size[head(u)]; }

    bool merge(int u, int v) {
        u = head(u), v = head(v);
        if (u == v) return false;
        if (-par_or_size[u] > -par_or_size[v]) std::swap(u, v);
        // subtree of u is smaller
        par_or_size[v] += par_or_size[u];
        par_or_size[u] = v;
        return true;
    }
};

