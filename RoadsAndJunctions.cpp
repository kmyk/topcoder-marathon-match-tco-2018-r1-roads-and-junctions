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



/******************************************************************************
 * general libraries
 ******************************************************************************/

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
    xor_shift_128(uint32_t seed = 42) {
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

constexpr double eps = 1e-8;



/******************************************************************************
 * utilities for this problem
 ******************************************************************************/

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

double compute_cost_of_spanning_tree_kruskal_body(
        vector<point_t> const & cities,
        vector<point_t> const & junctions,
        vector<tuple<double, int, int> > const & city_tree,
        vector<tuple<double, int, int> > const & que) {
    int NC = cities.size();
    int NJ = junctions.size();
    double acc = 0;
    static union_find_tree uft;
    uft.data.assign(NC + NJ, -1);
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

double compute_cost_of_spanning_tree_kruskal(
        vector<point_t> const & cities,
        vector<point_t> const & junctions,
        vector<tuple<double, int, int> > const & city_tree) {
    int NC = cities.size();
    int NJ = junctions.size();
    vector<tuple<double, int, int> > que;
    REP (a, NC) REP (b, NJ) {
        que.emplace_back(calc_distance(cities[a], junctions[b]), a, NC + b);
    }
    REP (a, NJ) REP (b, a) {
        que.emplace_back(calc_distance(junctions[a], junctions[b]), NC + a, NC + b);
    }
    sort(ALL(que));
    return compute_cost_of_spanning_tree_kruskal_body(cities, junctions, city_tree, que);
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

vector<pair<int, int> > drop_first_key(vector<tuple<double, int, int> > const & xs) {
    vector<pair<int, int> > ys;
    for (auto const & x : xs) {
        ys.emplace_back(get<1>(x), get<2>(x));
    }
    return ys;
}

const int neighborhood4_y[] = { 0, -1,  0, 1 };
const int neighborhood4_x[] = { 1,  0, -1, 0 };
const int neighborhood8_y[] = { 0, -1, -1, -1, 0, 1, 1, 1 };
const int neighborhood8_x[] = { 1, 1, 0, -1, -1, -1, 0, 1 };


/******************************************************************************
 * solver class
 ******************************************************************************/

class RoadsAndJunctions {

double clock_begin;
xor_shift_128 gen;
int S;
vector<point_t> cities;
double junction_cost;
double failure_probability;
int NC;
vector<point_t> junctions;
vector<bool> junction_status;
vector<pair<int, int> > roads;



/******************************************************************************
 * entry points
 ******************************************************************************/

public:
vector<int> buildJunctions(int a_S, vector<int> a_cities, double a_junctionCost, double a_failureProbability) {
    clock_begin = rdtsc();
    random_device device;
    gen = (decltype(gen))(device());
    S = a_S;
    cities = pack_points(a_cities);
    junction_cost = a_junctionCost;
    failure_probability = a_failureProbability;
    NC = cities.size();
    solve();
    return unpack_points(junctions);
}

vector<int> buildRoads(vector<int> a_junctionStatus) {
    junction_status = vector<bool>(ALL(a_junctionStatus));
    roads = with_junctions_status(NC, junctions, junction_status, [&](vector<point_t> const & junctions) {
        vector<point_t> points;
        copy(ALL(cities), back_inserter(points));
        copy(ALL(junctions), back_inserter(points));
        vector<tuple<double, int, int> > weighted_edges = construct_spanning_tree_prim(points);
        return drop_first_key(weighted_edges);
    });
    debug_print();
    return unpack_pairs(roads);
}



/******************************************************************************
 * solve() function
 ******************************************************************************/

private:

#ifdef LOCAL
ll seed;
#endif
vector<tuple<double, int, int> > city_tree;
double reference_score;
void solve() {
    // prepare
    city_tree = construct_spanning_tree_prim(cities);
    reference_score = 0;
    for (auto const & edge : city_tree) {
        reference_score += get<0>(edge);
    }

    // print debug info
#ifdef LOCAL
    seed = -1;
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
    build_one_junction_many_times();
    vector<pair<double, point_t> > candidates = make_candidates_from_memo(compute_cost_of_spanning_tree_memo_1);
    cerr << "size of candidates for hill climbing = " << candidates.size() << endl;

    // hill climbing to filter candidates
    candidates = filter_candidates_with_hill_climbing(candidates);
    cerr << "size of candidates after hill climbing = " << candidates.size() << endl;

    { // remove not good points
        auto it = lower_bound(ALL(candidates), make_pair(reference_score + eps, (point_t) { -1, -1 }));
        candidates.erase(it, candidates.end());
        cerr << "size of candidates after simple pruning = " << candidates.size() << endl;
    }

    // find pairs of points which make better score
    vector<pair<int, int> > conflicteds = enumerate_conflicted_pairs(candidates);
    cerr << "size of conflicted pairs = " << conflicteds.size() << endl;
    vector<tuple<double, point_t, point_t> > pair_candidates = list_seed_of_pair_candidates(candidates, conflicteds);
    pair_candidates = find_pair_candidates_with_hill_climbing(pair_candidates);
    cerr << "size of pair candidates = " << pair_candidates.size() << endl;

    // prune unused candidates
    candidates.erase(remove_if(ALL(candidates), [&](pair<double, point_t> const & candidate) {
        point_t p = candidate.second;
        int k = compute_neighborhoods_for_candidate_point(p).first.size();
        return k == 0;
    }), candidates.end());
    pair_candidates.erase(remove_if(ALL(pair_candidates), [&](tuple<double, point_t, point_t> const & candidate) {
        point_t p, q; tie(ignore, p, q) = candidate;
        int k = compute_neighborhoods_for_candidate_point_pair(p, q).first.size();
        return k == 0;
    }), pair_candidates.end());

    // check conflicts with simulated annealing
    tie(candidates, pair_candidates) = select_candidates_with_simulated_annealing(candidates, pair_candidates);
    cerr << "size of candidates after simulated annealing = " << candidates.size() << " + " << pair_candidates.size() << endl;

    // construct the result
    auto a = list_junctions_from_candidates(candidates);
    auto b = list_junctions_from_pair_candidates(pair_candidates);
    copy(ALL(a), back_inserter(junctions));
    copy(ALL(b), back_inserter(junctions));

    // normalize junctions
    sort(ALL(junctions));
    junctions.erase(unique(ALL(junctions)), junctions.end());
    junctions.erase(remove_if(ALL(junctions), [&](point_t p) { return count(ALL(cities), p); }), junctions.end());
    if ((int)junctions.size() > 2 * NC) {
        junctions.resize(2 * NC);
#ifdef LOCAL
        assert (false);
#endif
    }
}



/******************************************************************************
 * parts of solve() function
 ******************************************************************************/

map<point_t, double> compute_cost_of_spanning_tree_memo_1;
double compute_cost_of_spanning_tree_1(point_t p) {
    auto & memo = compute_cost_of_spanning_tree_memo_1;
    if (memo.count(p)) return memo[p];
    vector<point_t> junctions(1, p);
    return memo[p] = compute_cost_of_spanning_tree_kruskal(cities, junctions, city_tree);
}

map<pair<point_t, point_t>, double> compute_cost_of_spanning_tree_memo_2;
double compute_cost_of_spanning_tree_2(point_t p, point_t q) {
    auto & memo = compute_cost_of_spanning_tree_memo_2;
    if (q < p) swap(p, q);
    auto key = make_pair(p, q);
    if (memo.count(key)) return memo[key];
    vector<point_t> junctions({ p, q });
    return memo[key] = compute_cost_of_spanning_tree_kruskal(cities, junctions, city_tree);
}

void build_one_junction_many_times() {
    REP (a, NC) {
        vector<int> order(a);
        iota(ALL(order), 0);
        int size = min<int>(order.size(), 30);
        partial_sort(order.begin(), order.begin() + size, order.end(), [&](int b, int c) {
            return calc_distance(cities[a], cities[b]) < calc_distance(cities[a], cities[c]);
        });
        REP (j, size) REP (k, j) {
            int b = order[j];
            int c = order[k];
            int y = round((cities[a].y + cities[b].y + cities[c].y) / 3.0);
            int x = round((cities[a].x + cities[b].x + cities[c].x) / 3.0);
            compute_cost_of_spanning_tree_1((point_t) { y, x });
        }
    }
    REP (a, NC) {
        REP3 (y, max(0, cities[a].y - 1), min(S, cities[a].y + 1) + 1) {
            REP3 (x, max(0, cities[a].x - 1), min(S, cities[a].x + 1) + 1) {
                compute_cost_of_spanning_tree_1((point_t) { y, x });
            }
        }
    }
}

vector<pair<double, point_t> > make_candidates_from_memo(map<point_t, double> const & memo) {
    vector<pair<double, point_t> > candidates;
    for (auto const & it : memo) {
        double score = it.second;
        if (score < reference_score - eps) {
            candidates.emplace_back(score, it.first);
        }
    }
    sort(ALL(candidates));
    return candidates;
}

vector<pair<double, point_t> > filter_candidates_with_hill_climbing(vector<pair<double, point_t> > const & candidates) {
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
                double nscore = compute_cost_of_spanning_tree_1(np);
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
    sort(ALL(next_candidates));
    next_candidates.erase(unique(ALL(next_candidates)), next_candidates.end());  // while it contains floating points, it's OK since they are memoized value
    return next_candidates;
}

vector<pair<int, int> > enumerate_conflicted_pairs(vector<pair<double, point_t> > const & candidates) {
    vector<pair<int, int> > conflicteds;
    REP (i, candidates.size()) {
        REP (j, i) {
            double score_i; point_t p; tie(score_i, p) = candidates[i];
            double score_j; point_t q; tie(score_j, q) = candidates[j];
            double score = compute_cost_of_spanning_tree_2(p, q);
            double expected_delta = (reference_score - score_i) + (reference_score - score_j);
            double actual_delta = reference_score - score;
            if (actual_delta < expected_delta - eps) {
                conflicteds.emplace_back(i, j);
            }
        }
    }
    return conflicteds;
}

vector<tuple<double, point_t, point_t> > list_seed_of_pair_candidates(vector<pair<double, point_t> > const & candidates, vector<pair<int, int> > const & conflicteds) {
    vector<tuple<double, point_t, point_t> > pair_candidates;
    for (auto conflicted : conflicteds) {
        int i, j; tie(i, j) = conflicted;
        point_t p = candidates[i].second;
        point_t q = candidates[j].second;
        double score = compute_cost_of_spanning_tree_2(p, q);
        pair_candidates.emplace_back(score, p, q);
    }
    for (auto candidate : candidates) {
        point_t p = candidate.second;
        vector<int> order(NC);
        iota(ALL(order), 0);
        int size = min(NC, 30);
        partial_sort(order.begin(), order.begin() + size, order.end(), [&](int a, int b) {
            return calc_distance(p, cities[a]) < calc_distance(p, cities[b]);
        });
        order.resize(size);
        for (int a : order) {
            point_t q = cities[a];
            q.y += (q.y < p.y ? +1 : p.y < q.y ? -1 : 0);
            q.x += (q.x < p.x ? +1 : p.x < q.x ? -1 : 0);
            double score = compute_cost_of_spanning_tree_2(p, q);
            pair_candidates.emplace_back(score, p, q);
        }
    }
    return pair_candidates;
}

vector<tuple<double, point_t, point_t> > find_pair_candidates_with_hill_climbing(vector<tuple<double, point_t, point_t> > const & prev_pair_candidates) {
    vector<tuple<double, point_t, point_t> > pair_candidates;
    array<int, 8> order;
    iota(ALL(order), 0);
    for (auto candidate : prev_pair_candidates) {
        double score; point_t p, q; tie(score, p, q) = candidate;
        while (true) {
            shuffle(ALL(order), gen);
            bool found = false;
            for (int i : order) {
                point_t np = p;
                point_t nq = q;
                if (i >= 4) {
                    i -= 4;
                    swap(np, nq);
                }
                np.y += neighborhood4_y[i];
                np.x += neighborhood4_x[i];
                double next_score = compute_cost_of_spanning_tree_2(np, nq);
                if (next_score < score - eps) {
                    p = np;
                    q = nq;
                    score = next_score;
                    found = true;
                }
            }
            if (not found) break;
        }
        double score_1 = compute_cost_of_spanning_tree_1(p);
        double score_2 = compute_cost_of_spanning_tree_1(q);
        if (reference_score - score > (reference_score - score_1) + (reference_score - score_2) + eps + 0.1 * junction_cost) {
            pair_candidates.emplace_back(score, p, q);
        }
    }
    sort(ALL(pair_candidates));
    pair_candidates.erase(unique(ALL(pair_candidates)), pair_candidates.end());
    return pair_candidates;
}

vector<vector<pair<double, int> > > prepare_near_cities(vector<pair<double, point_t> > const & candidates, vector<tuple<double, point_t, point_t> > const & pair_candidates) {
    vector<vector<pair<double, int> > > near_cities;
    const int size = min(NC, 20);
    auto add = [&](point_t p) {
        near_cities.emplace_back();
        auto & it = near_cities.back();
        REP (a, NC) {
            it.emplace_back(calc_distance(p, cities[a]), a);
        }
        sort(ALL(it));
        it.resize(size);
        it.shrink_to_fit();
    };
    for (auto candidate : candidates) {
        add(candidate.second);
    }
    for (auto candidate : pair_candidates) {
        add(get<1>(candidate));
        add(get<2>(candidate));
    }
    return near_cities;
}

vector<tuple<double, int, int> > make_queue_from_near_cities(
        vector<point_t> const & junctions,
        vector<int> const & junction_indices,
        vector<vector<pair<double, int> > > const & near_cities) {
    int NC = cities.size();
    int NJ = junctions.size();
    vector<tuple<double, int, int> > que;
    REP (b, NJ) {
        for (auto const & it : near_cities[junction_indices[b]]) {
            double dist; int a; tie(dist, a) = it;
            que.emplace_back(dist, a, NC + b);
        }
    }
    REP (a, NJ) REP (b, a) {
        que.emplace_back(calc_distance(junctions[a], junctions[b]), NC + a, NC + b);
    }
    sort(ALL(que));
    return que;
}

pair<vector<pair<double, point_t> >, vector<tuple<double, point_t, point_t> > > select_candidates_with_simulated_annealing(vector<pair<double, point_t> > const & candidates, vector<tuple<double, point_t, point_t> > const & pair_candidates) {
    // prepare
    int NJ1 = candidates.size();
    int NJ2 = pair_candidates.size();
    int NJ = NJ1 + NJ2;
    if (NJ == 0) return make_pair(vector<pair<double, point_t> >(), vector<tuple<double, point_t, point_t> >());

#ifdef DEBUG
    map<vector<bool>, double> memo;
#else
    unordered_map<vector<bool>, double> memo;
#endif
    vector<vector<pair<double, int> > > near_cities = prepare_near_cities(candidates, pair_candidates);
    auto compute_cost_of_spanning_tree_mask = [&](vector<bool> const & mask) {
        if (memo.count(mask)) return memo[mask];
        vector<point_t> junctions;
        vector<int> junction_indices;
        REP (i, NJ1) if (mask[i]) {
            junctions.push_back(candidates[i].second);
            junction_indices.push_back(i);
        }
        REP (i, NJ2) if (mask[NJ1 + i]) {
            junctions.push_back(get<1>(pair_candidates[i]));
            junctions.push_back(get<2>(pair_candidates[i]));
            junction_indices.push_back(NJ1 + 2 * i);
            junction_indices.push_back(NJ1 + 2 * i + 1);
        }
        vector<tuple<double, int, int> > que = make_queue_from_near_cities(junctions, junction_indices, near_cities);
        return memo[mask] = compute_cost_of_spanning_tree_kruskal_body(cities, junctions, city_tree, que);
    };
    vector<bool> mask(NJ);
    double score = reference_score;
    vector<bool> result = mask;
    double highscore = score;

    // use idea of tabu search
    vector<vector<bool> > tabu_list(100);
    int tabu_head = 0;

    // run iterations
    int iteration = 0;
    double sa_clock_begin = rdtsc();
    double sa_clock_end = clock_begin + TLE * 0.95;
    double temperature = 1;
    for (; ; ++ iteration) {
        if (temperature < 0.1 or iteration % 100 == 0) {
            double t = rdtsc();
            if (t > sa_clock_end) break;
            temperature = 1 - (t - sa_clock_begin) / (sa_clock_end - sa_clock_begin);
        }
        int i = uniform_int_distribution<int>(0, NJ - 1)(gen);
        mask[i] = not mask[i];
        int j = -1;
        if (bernoulli_distribution(0.7)(gen)) {
            j = uniform_int_distribution<int>(0, NJ - 1)(gen);
            if (j == i) {
                j = -1;
            } else {
                mask[j] = not mask[j];
            }
        }
        int k = -1;
        if (j != -1 and bernoulli_distribution(0.7)(gen)) {
            k = uniform_int_distribution<int>(0, NJ - 1)(gen);
            if (k == i or k == j) {
                k = -1;
            } else {
                mask[k] = not mask[k];
            }
        }
        double next_score = compute_cost_of_spanning_tree_mask(mask);
        next_score += count(mask.begin(), mask.begin() + NJ1, true) * junction_cost;
        next_score += count(mask.begin() + NJ1, mask.end(),   true) * junction_cost * (2 + 0.5);  // pair candidates have much cost since both of the pair need to be successfully built
        double delta = score - next_score;
        bool accept = next_score < score + eps or bernoulli_distribution(max(0.0, min(1.0, 0.03 * exp(delta) / temperature)))(gen);
        if (accept and find(ALL(tabu_list), mask) != tabu_list.end()) accept = false;
        if (accept) {
            score = next_score;
            if (score < highscore) {
                highscore = score;
                result = mask;
                cerr << "simulated annealing: ";
                REP (i, NJ1) cerr << int(mask[i]);
                cerr << '.';
                REP3 (i, NJ1, NJ) cerr << int(mask[i]);
                cerr << " " << reference_score - highscore << endl;
            }
            tabu_list[tabu_head ++] = mask;
            if (tabu_head == (int)tabu_list.size()) tabu_head = 0;
        } else {
            mask[i] = not mask[i];
            if (j != -1) mask[j] = not mask[j];
            if (k != -1) mask[k] = not mask[k];
        }
    }
    cerr << "simulated annealing: iteration = " << iteration << endl;

    // make the result
    mask = result;
    vector<pair<double, point_t> > next_candidates;
    REP (i, NJ1) if (mask[i]) {
        next_candidates.push_back(candidates[i]);
    }
    vector<tuple<double, point_t, point_t> > next_pair_candidates;
    REP (i, NJ2) if (mask[NJ1 + i]) {
        next_pair_candidates.push_back(pair_candidates[i]);
    }
    return make_pair(next_candidates, next_pair_candidates);
}

pair<vector<point_t>, double> compute_neighborhoods_for_candidate_point(point_t p) {
    // sort neighborhoods with score
    vector<pair<double, point_t> > neighborhoods;
    constexpr int radius = 3;
    REP3 (ny, max(0, p.y - radius), min(S, p.y + radius) + 1) {
        REP3 (nx, max(0, p.x - radius), min(S, p.x + radius) + 1) {
            point_t np = { ny, nx };
            double score = compute_cost_of_spanning_tree_1(np);
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

    // make result
    vector<point_t> result;
    REP (i, k) {
        result.push_back(neighborhoods[i].second);
    }
    return make_pair(result, modified_scores[k]);
}

vector<point_t> list_junctions_from_candidates(vector<pair<double, point_t> > const & candidates) {
    vector<point_t> junctions;
    for (auto const & candidate : candidates) {
        point_t p = candidate.second;
        vector<point_t> qs; double modified_score; tie(qs, modified_score) = compute_neighborhoods_for_candidate_point(p);
        int k = qs.size();
        cerr << "use " << p << " (k = " << k << ", expected gain = " << reference_score - modified_score << ")" << endl;
#ifdef LOCAL
        assert (k != 0);
#endif
        copy(ALL(qs), back_inserter(junctions));
    }
    return junctions;
}

pair<vector<point_t>, double> compute_neighborhoods_for_candidate_point_pair(point_t p, point_t q) {
    // sort neighborhoods with score
    vector<tuple<double, point_t, point_t> > neighborhoods;
    constexpr int radius = 2;
    REP3 (npy, max(0, p.y - radius), min(S, p.y + radius) + 1) {
        REP3 (npx, max(0, p.x - radius), min(S, p.x + radius) + 1) {
            point_t np = { npy, npx };
            REP3 (nqy, max(0, q.y - radius), min(S, q.y + radius) + 1) {
                REP3 (nqx, max(0, q.x - radius), min(S, q.x + radius) + 1) {
                    point_t nq = { nqy, nqx };
                    double score = compute_cost_of_spanning_tree_2(np, nq);
                    neighborhoods.emplace_back(score, np, nq);
                }
            }
        }
    }
    sort(ALL(neighborhoods));

    // choose k
    vector<double> modified_scores;
    REP (k, neighborhoods.size()) {
        double p1 = pow(failure_probability, k);
        double p = 1 - pow(1 - p1, 2);
        double score = get<0>(neighborhoods[k]);
        double connection_cost = 0.2 * max(0, k - 1);
        modified_scores.push_back((1 - p) * score + p * reference_score + 2 * k * junction_cost + connection_cost);
    }
    int k = min_element(ALL(modified_scores)) - modified_scores.begin();

    // use best k neighborhoods
    vector<point_t> result;
    set<point_t> used_p;
    set<point_t> used_q;
    for (auto neighborhood : neighborhoods) {
        point_t np, nq; tie(ignore, np, nq) = neighborhood;
        if ((int)used_p.size() < k and not used_p.count(np)) {
            used_p.insert(np);
            result.push_back(np);
        }
        if ((int)used_q.size() < k and not used_q.count(nq)) {
            used_q.insert(nq);
            result.push_back(nq);
        }
    }
    return make_pair(result, modified_scores[k]);
}

vector<point_t> list_junctions_from_pair_candidates(vector<tuple<double, point_t, point_t> > const & pair_candidates) {
    vector<point_t> junctions;
    for (auto const & candidate : pair_candidates) {
        point_t p, q; tie(ignore, p, q) = candidate;
        vector<point_t> rs; double modified_score; tie(rs, modified_score) = compute_neighborhoods_for_candidate_point_pair(p, q);
        int k = rs.size();
        cerr << "use " << p << " " << q << " (k = " << k << ", expected gain = " << reference_score - modified_score << ")" << endl;
#ifdef LOCAL
        assert (k != 0);
#endif
        copy(ALL(rs), back_inserter(junctions));
    }
    return junctions;
}



/******************************************************************************
 * debug_print
 ******************************************************************************/

void debug_print() {
    double elapsed = rdtsc() - clock_begin;
    int NJ = junctions.size();
    double score = NJ * junction_cost;
    for (auto road : roads) {
        int a, b; tie(a, b) = road;
        point_t p = a < NC ? cities[a] : junctions[a - NC];
        point_t q = b < NC ? cities[b] : junctions[b - NC];
        score += calc_distance(p, q);
    }

#ifdef LOCAL
    // measure average score
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
#endif

    auto format_float = [&](double x) { return round(x * 1e8) * 1e-8 + 0; };  // NOTE: "+ 0" removes negative zeros
    cerr << "score = " << score << endl;
    cerr << "NJ = " << NJ << endl;
    cerr << "raw reference delta = " << reference_score - score << endl;
#ifdef LOCAL
    cerr << "average score = " << reference_score - average_reference_delta << endl;
    cerr << "average reference delta = " << format_float(average_reference_delta) << endl;
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
/*
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
*/
    cerr << "}" << endl;
#endif
}

};
