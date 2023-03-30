#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <queue>
#include <stack>

#include "jngen/graph.h"

using namespace std;

#define rep(i, a, b) for (auto i = (a); i < (b); ++i)

auto find_triangles(Graph _g, int n) {
    vector<vector<int>> g(n);
    rep(u, 0, n) {
        auto verts = _g.edges(u);
        g[u] = {verts.begin(), verts.end()};
        sort(g[u].begin(), g[u].end());
    }

    std::vector<std::array<int, 3>> tris;
    rep(u, 0, n) {
        std::vector<int> cmb;
        for (auto v : g[u]) {
            cmb.clear();
            set_intersection(g[u].begin(), g[u].end(), g[v].begin(), g[v].end(),
                             back_inserter(cmb));
            for (auto w : cmb) tris.push_back({{u, v, w}});
        }
    }

    return tris;
}

int main(int argc, char** argv) {
    cin.tie(nullptr)->sync_with_stdio(false);

    assert(argc == 5);

    const int n = stoi(argv[1]), m = stoi(argv[2]);

    auto g = Graph::random(n, m).connected().g().shuffle();

    std::ofstream graph_file(argv[3]);

    std::cerr << "Graph to " << argv[3] << '\n';

    graph_file << n << ' ' << m << '\n';
    rep(u, 0, n) {
        auto neighbourhood = g.edges(u);

        sort(neighbourhood.begin(), neighbourhood.end());

        graph_file << u << ' ';
        graph_file << neighbourhood.size() << ' ';
        for (auto v : neighbourhood) graph_file << v << ' ';
        graph_file << '\n';
    }

    graph_file.close();

    std::ofstream tris_file(argv[4]);

    std::cerr << "Triangles to " << argv[4] << '\n';

    auto tris = find_triangles(g, n);
    sort(tris.begin(), tris.end());
    for (auto [u, v, w] : tris) tris_file << u << ' ' << v << ' ' << w << '\n';

    tris_file.close();
}
