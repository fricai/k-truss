#include <iostream>
#include <random>
#include <set>

using namespace std;

#define rep(i, a, b) for (auto i = (a); i < (b); ++i)

int main(int argc, char** argv) {
    cin.tie(nullptr)->sync_with_stdio(false);

    assert(argc == 3);

    random_device rd;
    mt19937 rng{rd()};

    int n = stoi(argv[1]), m = stoi(argv[2]);

    uniform_int_distribution<> dis{0, n - 1};

    set<pair<int, int>> s;

    vector<vector<int>> g(n);
    rep(i, 0, m) {
        int x, y;
        do {
            x = dis(rng);
            y = dis(rng);
            if (x > y) swap(x, y);
        } while (x == y || s.count({x, y}));
        g[x].push_back(y);
        g[y].push_back(x);
        s.insert({x, y});
    }

    rep(u, 0, n) {
        cout << g[u].size() << ' ';
        for (auto v : g[u]) cout << v << ' ';
        cout << '\n';
    }
}
