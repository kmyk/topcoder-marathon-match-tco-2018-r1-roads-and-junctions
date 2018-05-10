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

constexpr double ticks_per_sec = 2800000000;
constexpr double ticks_per_sec_inv = 1.0 / ticks_per_sec;
inline double rdtsc() { // in seconds
    uint32_t lo, hi;
    asm volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    return (((uint64_t)hi << 32) | lo) * ticks_per_sec_inv;
}
constexpr int TLE = 10; // sec

struct union_find_tree {
    vector<int> data;
    union_find_tree() = default;
    explicit union_find_tree(size_t n) : data(n, -1) {}
    bool is_root(int i) { return data[i] < 0; }
    int find_root(int i) { return is_root(i) ? i : (data[i] = find_root(data[i])); }
    int tree_size(int i) { return - data[find_root(i)]; }
    int unite_trees(int i, int j) {
        i = find_root(i); j = find_root(j);
        if (i != j) {
            if (tree_size(i) < tree_size(j)) swap(i,j);
            data[i] += data[j];
            data[j] = i;
        }
        return i;
    }
    bool is_same(int i, int j) { return find_root(i) == find_root(j); }
};

struct point_t { int y, x; };
inline bool operator <  (point_t const & a, point_t const & b) { return a.y != b.y ? a.y < b.y : a.x < b.x; };
inline bool operator == (point_t const & a, point_t const & b) { return a.y == b.y and a.x == b.x; };
inline ostream & operator << (ostream & out, point_t const & a) { return out << "(" << a.x << ", " << a.y << ")"; };

double calc_distance(point_t const & a, point_t const & b) {
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
            i = NC + index[i - NC];
        }
    };
    for (auto & road : roads) {
        func(road.first);
        func(road.second);
    }
    return roads;
}

/**
 * @note result edges are sorted with costs
 */
vector<tuple<double, int, int> > construct_spanning_tree_prim(vector<point_t> const & cities) {
    // prepare
    int NC = cities.size();
    reversed_priority_queue<tuple<double, int, int> > que;
    vector<bool> used(NC);
    auto use = [&](int a) {
        used[a] = true;
        REP (b, NC) if (not used[b]) {
            que.emplace(calc_distance(cities[a], cities[b]), a, b);
        }
    };

    // run Prim
    vector<tuple<double, int, int> > edges;
    use(0);
    REP (size, NC - 1) {
        double dist; int a, b;
        while (true) {
            assert (not que.empty());
            tie(dist, a, b) = que.top();
            que.pop();
            if (not used[b]) break;
        }
        edges.emplace_back(dist, a, b);
        use(b);
    }

    sort(ALL(edges));
    return edges;
}

double compute_cost_of_spanning_tree_kruskal(vector<point_t> const & cities, vector<point_t> const & junctions, vector<tuple<double, int, int> > const & city_tree) {
    int NC = cities.size();
    int NJ = junctions.size();

    // prepare que
    vector<tuple<double, int, int> > que;
    REP (a, NC) REP (b, NJ) {
        que.emplace_back(calc_distance(cities[a], junctions[b]), a, NC + b);
    }
    REP (a, NJ) REP (b, a) {
        que.emplace_back(calc_distance(junctions[a], junctions[b]), NC + a, NC + b);
    }
    sort(ALL(que));

    // run Kruskal with two pointers
    double acc = 0;
    union_find_tree uft(NC + NJ);
    int i = 0, j = 0;
    REP (size, NC + NJ - 1) {
        double dist_i = INFINITY;
        for (; i < (int)city_tree.size(); ++ i) {
            int a, b; tie(dist_i, a, b) = city_tree[i];
            if (not uft.is_same(a, b)) break;
        }
        double dist_j = INFINITY;
        for (; j < (int)que.size(); ++ j) {
            int a, b; tie(dist_j, a, b) = que[j];
            if (not uft.is_same(a, b)) break;
        }
        acc += min(dist_i, dist_j);
        int a, b;
        tie(ignore, a, b) = (dist_i <= dist_j ? city_tree[i ++] : que[j ++]);
        uft.unite_trees(a, b);
    }

    return acc;
}

pair<vector<point_t>, function<vector<pair<int, int> > (vector<bool> const &)> > solve(int S, vector<point_t> const & cities, double junction_cost, double failure_probability) {
    double clock_begin = rdtsc();

    // prepare
    int NC = cities.size();
    vector<tuple<double, int, int> > city_tree = construct_spanning_tree_prim(cities);
    double reference_score = 0;
    for (auto const & edge : city_tree) {
        reference_score += get<0>(edge);
    }

    // print debug info
#ifdef LOCAL
    ll seed = -1;
    if (getenv("SEED") != nullptr) {
        seed = atoll(getenv("SEED"));
        cerr << "seed = " << seed << endl;
    }
#endif
    cerr << "S = " << S << endl;
    cerr << "NC = " << NC << endl;
    cerr << "junction cost = " << junction_cost << endl;
    cerr << "failure probability = " << failure_probability << endl;
    cerr << "reference score = " << reference_score << endl;

    // build a junction to find candidates
    map<point_t, double> memo;
    auto compute_spanning_tree_with_junction = [&](int y, int x) {
        point_t key = { y, x };
        if (memo.count(key)) return memo[key];
        vector<point_t> junctions(1, (point_t) { y, x });
        return memo[key] = compute_cost_of_spanning_tree_kruskal(cities, junctions, city_tree);
    };
    REP (i, NC) REP (j, i) REP (k, j) {
        int y = round((cities[i].y + cities[j].y + cities[k].y) / 3.0);
        int x = round((cities[i].x + cities[j].x + cities[k].x) / 3.0);
        compute_spanning_tree_with_junction(y, x);
    }
    vector<pair<double, point_t> > candidates;
    for (auto const & it : memo) {
        double score = it.second;
        if (score < reference_score - junction_cost / 10.0) {  // this border is magical
            candidates.emplace_back(score, it.first);
        }
    }
    sort(ALL(candidates));
    cerr << "size of candidates = " << candidates.size() << endl;

    // hill climbing
    vector<point_t> junctions;
    const int dy[] = { -1, 1, 0, 0 };
    const int dx[] = { 0, 0, 1, -1 };
    for (auto const & candidate : candidates) {
        double score; point_t p; tie(score, p) = candidate;
// cerr << "new " << p << " -> " << score<< endl;
        while (true) {
            bool found = false;
            REP (i, 4) {
                int ny = max(0, min(S, p.y + dy[i]));
                int nx = max(0, min(S, p.x + dx[i]));
                double nscore = compute_spanning_tree_with_junction(ny, nx);
                if (nscore < score) {
                    found = true;
                    score = nscore;
                    p.y = ny;
                    p.x = nx;
// cerr << "move " << p << " -> " << score<< endl;
                    break;
                }
            }
            if (not found) {
                break;
            }
        }
        vector<double> modified_score;
        REP (k, 5 + 1) {  // 5 = 1 + 4 neighborhoods
            double p = pow(failure_probability, k);
            modified_score.push_back((1 - p) * score + p * reference_score + k * junction_cost);
        }
        int k = min_element(ALL(modified_score)) - modified_score.begin();
// if (k == 0) cerr << "rejected " << p << " -> " << reference_score - score << " (k = " << k << ")" << endl;
        if (k == 0) continue;
        bool found = false;
        for (point_t q : junctions) {
            if (calc_distance(p, q) < 30) {
                found = true;
                break;
            }
        }
        if (not found) {
            cerr << "accepted " << p << " -> " << reference_score - score << " (k = " << k << ")" << endl;
            junctions.push_back(p);
            REP (i, k - 1) {
                point_t np = { p.y + dy[i], p.x + dx[i] };
                junctions.push_back(np);
            }
        }
    }
    sort(ALL(junctions));
    junctions.erase(unique(ALL(junctions)), junctions.end());

    // solve roads
    auto cont = [=](vector<bool> const & junction_status) {
        return with_junctions_status(NC, junctions, junction_status, [&](vector<point_t> const & junctions) {
            // use all junctions
            vector<point_t> points;
            copy(ALL(cities), back_inserter(points));
            copy(ALL(junctions), back_inserter(points));
            vector<tuple<double, int, int> > weighted_edges = construct_spanning_tree_prim(points);
            vector<pair<int, int> > edges;
            double score = 0;
            for (auto const & edge : weighted_edges) {
                double cost; int a, b; tie(cost, a, b) = edge;
                score += cost;
                edges.emplace_back(a, b);
            }

            // print debug info
            double elapsed = rdtsc() - clock_begin;
            cerr << "score = " << score << endl;
            cerr << "reference delta = " << reference_score - score << endl;
            cerr << "elapsed = " << elapsed << endl;
#ifdef LOCAL
            cerr << "{\"seed\":" << seed
                 << ",\"S\":" << S
                 << ",\"NC\":" << NC
                 << ",\"junction_cost\":" << junction_cost
                 << ",\"failure_probability\":" << failure_probability
                 << ",\"reference_score\":" << reference_score
                 << ",\"score\":" << score
                 << ",\"NJ\":" << junction_status.size()
                 << ",\"reference_delta\":" << reference_score - score
                 << ",\"elapsed\":" << elapsed
                 << "}" << endl;
#endif
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
