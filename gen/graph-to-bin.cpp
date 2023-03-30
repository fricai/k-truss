#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#define rep(i, a, b) for (auto i{a}; i < (b); ++i)

using namespace std;

int main(int argc, char* argv[]) {
    assert(argc == 3);
    // header, graph

    using u32 = uint32_t;
    u32 n, m;
    cin >> n >> m;

    vector<vector<u32>> g(n);
    rep(i, 0u, n) {
        u32 id, deg;
        cin >> id >> deg;
        g[id].resize(deg + 2);
        g[id][0] = id;
        g[id][1] = deg;
        rep(i, 0u, deg) cin >> g[id][i + 2];
        assert(is_sorted(g[id].begin() + 2, g[id].end()));
    }

    ofstream hf(argv[1], ios::binary);
    ofstream gf(argv[2], ios::binary);

    gf.write((char*)&n, 4);
    gf.write((char*)&m, 4);
    uint32_t ctr = 8;
    rep(id, 0u, n) {
        gf.write((char*)g[id].data(), 4 * g[id].size());
        hf.write((char*)&ctr, 4);
        ctr += 4 * g[id].size();
    }
}
