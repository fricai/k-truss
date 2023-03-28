#include <iostream>
#include <vector>
#include <fstream>

#define rep(i, a, b) for (auto i{a}; i < (b); ++i)

using namespace std;

int main(int argc, char* argv[]) {
    assert(argc == 3);
    // header, graph

    using u32 = uint32_t;
    u32 n, m;
    cin >> n >> m;


    std::cout << n << ' ' << m << '\n';

    vector<vector<u32>> g(n);
    rep(e, 0, m) {
        int u, v;
        cin >> u >> v;
        g[u].push_back(v);
        g[v].push_back(u);

        std::cout << "Processed edge " << e << '\n';
    }

    rep(u, 0, n) {
        sort(g[u].rbegin(), g[u].rend());
        g[u].push_back(g[u].size());
        g[u].push_back(u);
        reverse(g[u].begin(), g[u].end());
    }

    std::cout << "Data to write complete\n";

    ofstream hf(argv[1], ios::binary);
    ofstream gf(argv[2], ios::binary);

    gf.write((char*)&n, 4);
    gf.write((char*)&m, 4);
    uint32_t ctr = 8;
    rep(id, 0, n) {
        gf.write((char*)g[id].data(), 4 * g[id].size());
        hf.write((char*)&ctr, 4);
        ctr += 4 * g[id].size();
    }
}
