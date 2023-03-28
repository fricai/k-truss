#include <fstream>
#include <iostream>
#include <cassert>
#include <queue>
#include "jngen/graph.h"
#include <stack>
#include <array>

using namespace std;

#define rep(i, a, b) for (auto i = (a); i < (b); ++i)

using vertex_t = int;
using edge_t = std::array<vertex_t, 2>;

edge_t make_edge(int x, int y) {
    if (x > y) std::swap(x, y);
    return {x, y};
}

void swap(uint64_t* a, uint64_t* b) {
    uint64_t temp = *a;
    *a = *b;
    *b = temp;
}

uint64_t edge(uint64_t a, uint64_t b) {
    if (a > b) swap(&a, &b);
    return (a << 32) | b;
}

void next_step(vector<unordered_set<int>>& V, int k, map<edge_t, int>& KTE) {
    int n = V.size();

    // prefilter
    queue<int> verts_to_delete;
    for (int i = 0; i < n; i++) {
        if ((int)V[i].size() <= k) verts_to_delete.push(i);
    }
    while (!verts_to_delete.empty()) {
        int v = verts_to_delete.front();
        verts_to_delete.pop();

        for (int u : V[v]) {
            KTE[make_edge(v, u)] = k - 1;

            V[u].erase(v);
            if ((int)V[u].size() <= k) verts_to_delete.push(u);
        }
        V[v].clear();
    }

    // initialize
    queue<edge_t> edges_to_delete;
    unordered_map<uint64_t, int> count;

    for (int i = 0; i < n; i++) {
        for (int j : V[i]) {
            // edge (i,j)
            if (j < i) continue;
            int c = 0;
            for (int v : V[j]) {
                if (V[i].find(v) != V[i].end()) c++;
            }
            count[edge(i, j)] = c;
            if (c < k) edges_to_delete.push({i, j});
        }
    }

    // filter
    while (!edges_to_delete.empty()) {
        auto e = edges_to_delete.front();
        edges_to_delete.pop();
        V[e[0]].erase(e[1]);
        V[e[1]].erase(e[0]);

        KTE[make_edge(e[0], e[1])] = k - 1;

        for (int v : V[e[1]]) {
            if (V[e[0]].find(v) == V[e[0]].end()) continue;
            uint64_t e1 = edge(e[0], v);
            uint64_t e2 = edge(e[1], v);
            count[e1] -= 1;
            count[e2] -= 1;
            if (count[e1] == k - 1) edges_to_delete.push({e[0], v});
            if (count[e2] == k - 1) edges_to_delete.push({e[1], v});
        }
    }
}

auto vanilla_KTE(Graph g, int n, int max_k) {
    vector<unordered_set<int>> V(n);
    rep(u, 0, n) {
        auto verts = g.edges(u);
        V[u] = unordered_set<int>(verts.begin(), verts.end());
    }

    // now find all k-trusses
    // A group of size k is equivalent to a truss of size k+2
    // for each k
    map<edge_t, int> KTE;

    rep(k, 1, max_k + 1) next_step(V, k, KTE);

    return KTE;
}

int main(int argc, char** argv) {
    cin.tie(nullptr)->sync_with_stdio(false);

    assert(argc == 7);

    random_device rd;
    mt19937 rng{rd()};

    int n = stoi(argv[1]), m = stoi(argv[2]);
    int k1 = stoi(argv[3]), k2 = stoi(argv[4]);

    auto g = Graph::random(n, m).connected().g().shuffle();

    std::ofstream graph_file(argv[5]);
    std::ofstream truss_file(argv[6]);
    graph_file << n << ' ' << m << '\n';
    rep(u, 0, n) {
        const auto neighbourhood = g.edges(u);

        graph_file << u << ' ';
        graph_file << neighbourhood.size() << ' ';
        for (auto v : neighbourhood) graph_file << v << ' ';
        graph_file << '\n';
    }
    graph_file.close();

    auto kte = vanilla_KTE(g, n, k2);
    for (auto p : kte) {
        int trussness = p.second < k1 ? k1 - 1 : p.second;
        truss_file << p.first[0] << ' ' << p.first[1] << ' ' << trussness << '\n';
    }

    rep(u, 0, n) {
        const auto nbd = g.edges(u);
        for (auto v : nbd) 
            if (v > u and !kte.count(make_edge(u, v)))
                truss_file << u << ' ' << v << ' ' << k2 << '\n';
    }
    truss_file.close();
}
