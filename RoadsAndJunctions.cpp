#include <bits/stdc++.h>
#define REP(i, n) for (int i = 0; (i) < int(n); ++ (i))
#define REP3(i, m, n) for (int i = (m); (i) < int(n); ++ (i))
#define REP_R(i, n) for (int i = int(n) - 1; (i) >= 0; -- (i))
#define REP3R(i, m, n) for (int i = int(n) - 1; (i) >= int(m); -- (i))
#define ALL(x) begin(x), end(x)
using ll = long long;
using namespace std;
template <class T> using reversed_priority_queue = priority_queue<T, vector<T>, greater<T> >;
template <class T> inline void chmax(T & a, T const & b) { a = max(a, b); }
template <class T> inline void chmin(T & a, T const & b) { a = min(a, b); }
template <typename T> ostream & operator << (ostream & out, vector<T> const & xs) { REP (i, int(xs.size()) - 1) out << xs[i] << ' '; if (not xs.empty()) out << xs.back(); return out; }

struct point_t { int y, x; };

double hypot(point_t const & a, point_t const & b) {
    return hypot(b.y - a.y, b.x - a.x);
}

template <class Function>
vector<pair<int, int> > with_junctions_status(int NC, vector<point_t> const & junctions, vector<bool> const & junction_status, Function cont) {
    assert (junctions.size() == junction_status.size());
    int NJ = junctions.size();
    vector<point_t> next_junctions;
    vector<int> index;
    REP (i, NJ) if (junction_status[i]) {
        next_junctions.push_back(junctions[i]);
        index.push_back(i);
    }
    vector<pair<int, int> > roads = cont(next_junctions);
    auto func = [&](int & i) {
        assert (0 <= i and i < NC + (int)next_junctions.size());
        if (i >= NC) {
            i = index[i - NC];
        }
    };
    for (auto & road : roads) {
        func(road.first);
        func(road.second);
    }
    return roads;
}

double construct_spanning_tree_prim(vector<point_t> const & cities, vector<point_t> const & junctions, vector<pair<int, int> > & edges) {
    assert (edges.empty());
    int NC = cities.size();
    int NJ = junctions.size();
    reversed_priority_queue<tuple<double, int, int> > que;
    vector<bool> used(NC + NJ);
    auto use = [&](int i) {
        used[i] = true;
        point_t p = (i < NC ? cities[i] : junctions[i - NC]);
        REP (j, NC + NJ) if (not used[j]) {
            point_t q = (j < NC ? cities[j] : junctions[j - NC]);
            que.emplace(hypot(p, q), i, j);
        }
    };
    double sum_dist = 0;
    use(0);
    REP (size, NC + NJ - 1) {
        double dist; int i, j;
        while (true) {
            assert (not que.empty());
            tie(dist, i, j) = que.top();
            que.pop();
            if (not used[j]) break;
        }
        edges.emplace_back(i, j);
        sum_dist += dist;
        use(j);
    }
    return sum_dist;
}

pair<vector<point_t>, function<vector<pair<int, int> > (vector<bool> const &)> > solve(int S, vector<point_t> const & cities, double junction_cost, double failure_probability) {
    int NC = cities.size();
    vector<point_t> junctions;

    auto cont = [=](vector<bool> const & junction_status) {
        return with_junctions_status(NC, junctions, junction_status, [&](vector<point_t> const & junctions) {
            vector<pair<int, int> > edges;
            construct_spanning_tree_prim(cities, junctions, edges);

            return edges;
        });
    };
    return make_pair(junctions, cont);
}


vector<point_t> pack_points(vector<int> const & b) {
    assert (b.size() % 2 == 0);
    int n = b.size() / 2;
    vector<point_t> a(n);
    REP (i, n) {
        int x = b[2 * i    ];
        int y = b[2 * i + 1];
        a[i] = { y, x };
    }
    return a;
}
vector<int> unpack_points(vector<point_t> const & a) {
    int n = a.size();
    vector<int> b(2 * n);
    REP (i, n) {
        b[2 * i    ] = a[i].x;
        b[2 * i + 1] = a[i].y;
    }
    return b;
}
vector<int> unpack_pairs(vector<pair<int, int> > const & a) {
    int n = a.size();
    vector<int> b(2 * n);
    REP (i, n) {
        b[2 * i    ] = a[i].first;
        b[2 * i + 1] = a[i].second;
    }
    return b;
}

class RoadsAndJunctions {
    function<vector<pair<int, int> > (vector<bool> const &)> cont;
public:
    vector<int> buildJunctions(int S, vector<int> cities, double junctionCost, double failureProbability) {
        vector<point_t> junctions;
        tie(junctions, cont) = solve(S, pack_points(cities), junctionCost, failureProbability);
        return unpack_points(junctions);
    }
    vector<int> buildRoads(vector<int> junctionStatus) {
        vector<pair<int, int> > roads = cont(vector<bool>(ALL(junctionStatus)));
        return unpack_pairs(roads);
    }
};
