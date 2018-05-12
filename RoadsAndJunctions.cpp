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

class xor_shift_128 {
public:
    typedef uint32_t result_type;
    xor_shift_128(uint32_t seed) {
        set_seed(seed);
    }
    void set_seed(uint32_t seed) {
        a = seed = 1812433253u * (seed ^ (seed >> 30));
        b = seed = 1812433253u * (seed ^ (seed >> 30)) + 1;
        c = seed = 1812433253u * (seed ^ (seed >> 30)) + 2;
        d = seed = 1812433253u * (seed ^ (seed >> 30)) + 3;
    }
    uint32_t operator() () {
        uint32_t t = (a ^ (a << 11));
        a = b; b = c; c = d;
        return d = (d ^ (d >> 19)) ^ (t ^ (t >> 8));
    }
    static constexpr uint32_t max() { return numeric_limits<result_type>::max(); }
    static constexpr uint32_t min() { return numeric_limits<result_type>::min(); }
private:
    uint32_t a, b, c, d;
};

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

const int neighborhood4_y[] = { 0, -1,  0, 1 };
const int neighborhood4_x[] = { 1,  0, -1, 0 };
const int neighborhood8_y[] = { 0, -1, -1, -1, 0, 1, 1, 1 };
const int neighborhood8_x[] = { 1, 1, 0, -1, -1, -1, 0, 1 };

constexpr double eps = 1e-8;
pair<vector<point_t>, function<vector<pair<int, int> > (vector<bool> const &)> > solve(int S, vector<point_t> const & cities, double junction_cost, double failure_probability) {
    double clock_begin = rdtsc();
    xor_shift_128 gen(42);

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

    // build a junction to find candidates for hill climbing
    map<point_t, double> memo;
    auto compute_cost_of_spanning_tree_memo = [&](point_t p) {
        if (memo.count(p)) return memo[p];
        vector<point_t> junctions(1, p);
        return memo[p] = compute_cost_of_spanning_tree_kruskal(cities, junctions, city_tree);
    };
    REP (a, NC) REP (b, a) REP (c, b) {
        int y = round((cities[a].y + cities[b].y + cities[c].y) / 3.0);
        int x = round((cities[a].x + cities[b].x + cities[c].x) / 3.0);
        compute_cost_of_spanning_tree_memo((point_t) { y, x });
    }
    REP (a, NC) {
        REP3 (y, max(0, cities[a].y - 1), min(S, cities[a].y + 1) + 1) {
            REP3 (x, max(0, cities[a].x - 1), min(S, cities[a].x + 1) + 1) {
                compute_cost_of_spanning_tree_memo((point_t) { y, x });
            }
        }
    }
    vector<pair<double, point_t> > candidates;
    for (auto const & it : memo) {
        double score = it.second;
        if (score < reference_score - junction_cost / 10.0) {  // this border is magical
            candidates.emplace_back(score, it.first);
        }
    }
    sort(ALL(candidates));
    cerr << "size of candidates for hill climbing = " << candidates.size() << endl;

    { // hill climbing to filter candidates
        vector<pair<double, point_t> > next_candidates;
        array<int, 8> order;
        iota(ALL(order), 0);
        for (auto const & candidate : candidates) {
            double score; point_t p; tie(score, p) = candidate;
            while (true) {
                bool found = false;
                shuffle(ALL(order), gen);
                for (int i : order) {
                    int ny = max(0, min(S, p.y + neighborhood8_y[i]));
                    int nx = max(0, min(S, p.x + neighborhood8_x[i]));
                    point_t np = { ny, nx };
                    double nscore = compute_cost_of_spanning_tree_memo(np);
                    if (nscore < score) {
                        found = true;
                        score = nscore;
                        p = np;
                        break;
                    }
                }
                if (not found) {
                    break;
                }
            }
            double modified_score_1 = (1 - failure_probability) * score + failure_probability * reference_score + junction_cost;
            if (modified_score_1 > reference_score) {
                continue;
            }
            next_candidates.emplace_back(score, p);
        }
        candidates.clear();
        candidates.swap(next_candidates);
        sort(ALL(candidates));
        candidates.erase(unique(ALL(candidates)), candidates.end());  // while it contains floating points, it's OK since they are memoized value
        cerr << "size of candidates after hill climbing = " << candidates.size() << endl;
    }

    { // check conflicts
        vector<pair<double, point_t> > next_candidates;
        vector<point_t> junctions;
        double score = reference_score;
        for (auto const & candidate : candidates) {
            point_t p = candidate.second;
            junctions.push_back(p);
            double next_score = compute_cost_of_spanning_tree_kruskal(cities, junctions, city_tree);
            if (next_score < score) {
                score = next_score;
                next_candidates.push_back(candidate);
            } else {
                junctions.pop_back();
            }
        }
        candidates.clear();
        candidates.swap(next_candidates);
        cerr << "size of candidates after removing conflicts = " << candidates.size() << endl;
    }

    // construct the result
    vector<point_t> junctions;
    for (auto const & candidate : candidates) {
        point_t p = candidate.second;

        // sort neighborhoods with score
        vector<pair<double, point_t> > neighborhoods;
        constexpr int radius = 3;
        REP3 (ny, max(0, p.y - radius), min(S, p.y + radius) + 1) {
            REP3 (nx, max(0, p.x - radius), min(S, p.x + radius) + 1) {
                point_t np = { ny, nx };
                double score = compute_cost_of_spanning_tree_memo(np);
                neighborhoods.emplace_back(score, np);
            }
        }
        sort(ALL(neighborhoods));

        // choose k
        vector<double> modified_scores;
        REP (k, neighborhoods.size()) {
            double p = pow(failure_probability, k);
            double score = neighborhoods[k].first;  // NOTE: this is an approximation
            double connection_cost = 0.1 * max(0, k - 1);  // NOTE: this is also an approximation
            modified_scores.push_back((1 - p) * score + p * reference_score + k * junction_cost + connection_cost);
        }
        int k = min_element(ALL(modified_scores)) - modified_scores.begin();
        if (k == 0) continue;
        cerr << "use " << p << " (k = " << k << ", expected gain = " << reference_score - modified_scores[k] << ")" << endl;

        // use best k neighborhoods
        REP (i, k) {
            junctions.push_back(neighborhoods[i].second);
        }
    }

    sort(ALL(junctions));
    junctions.erase(unique(ALL(junctions)), junctions.end());
    junctions.erase(remove_if(ALL(junctions), [&](point_t p) { return count(ALL(cities), p); }), junctions.end());
    if ((int)junctions.size() > 2 * NC) {
        junctions.resize(2 * NC);
#ifdef LOCAL
        assert (false);
#endif
    }

#ifdef LOCAL
    double local_elapsed = rdtsc();
    vector<double> samples; {
        int NJ = junctions.size();
        REP (iteration, 10000) {
            vector<point_t> constructed_junctions;
            for (point_t p : junctions) {
                if (not bernoulli_distribution(failure_probability)(gen)) {
                    constructed_junctions.push_back(p);
                }
            }
            double score = NJ * junction_cost + compute_cost_of_spanning_tree_kruskal(cities, constructed_junctions, city_tree);
            samples.push_back(reference_score - score);
        }
        sort(ALL(samples));
    }
    double average_reference_delta = accumulate(ALL(samples), 0.0) / samples.size();
    local_elapsed = rdtsc() - local_elapsed;
#endif

    // solve roads
    auto cont = [=](vector<bool> const & junction_status) {
        return with_junctions_status(NC, junctions, junction_status, [&](vector<point_t> const & junctions) {
            int NJ = junction_status.size();

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
            score += NJ * junction_cost;

            // print debug info
            double elapsed = rdtsc() - clock_begin;
            auto format_float = [&](double x) { return round(x * 1e8) * 1e-8 + 0; };  // NOTE: "+ 0" removes negative zeros
            cerr << "score = " << score << endl;
            cerr << "NJ = " << NJ << endl;
            cerr << "raw reference delta = " << reference_score - score << endl;
#ifdef LOCAL
            cerr << "average score = " << reference_score - average_reference_delta << endl;
            cerr << "average reference delta = " << format_float(average_reference_delta) << endl;
            elapsed -= local_elapsed;
#endif
            cerr << "elapsed = " << elapsed << endl;
#ifdef LOCAL
            cerr << "{\"seed\":" << seed
                 << ",\"S\":" << S
                 << ",\"NC\":" << NC
                 << ",\"junction_cost\":" << junction_cost
                 << ",\"failure_probability\":" << failure_probability
                 << ",\"reference_score\":" << reference_score
                 << ",\"score\":" << score
                 << ",\"NJ\":" << NJ
                 << ",\"raw_reference_delta\":" << reference_score - score
                 << ",\"average_score\":" << score - average_reference_delta
                 << ",\"average_reference_delta\":" << format_float(average_reference_delta)
                 << ",\"elapsed\":" << elapsed
                 ;
            cerr << ",\"delta_samples\":[";
            map<double, int> sample_count;
            for (double sample : samples) {
                sample_count[sample] += 1;
            }
            bool is_first = true;
            for (auto it : sample_count) {
                if (not is_first) cerr << ",";
                is_first = false;
                cerr << "[" << it.first << "," << it.second << "]";
            }
            cerr << "]";
            cerr << "}" << endl;
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
