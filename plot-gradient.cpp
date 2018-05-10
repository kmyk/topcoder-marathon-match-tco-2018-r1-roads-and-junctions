#include <bits/stdc++.h>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#define REP(i, n) for (int i = 0; (i) < int(n); ++ (i))
using namespace std;
using namespace boost;

int main() {
    // input
    int S; cin >> S;
    int NC; cin >> NC; NC /= 2;
    vector<pair<int, int> > cities(NC);
    REP (i, NC) {
        cin >> cities[i].first;
        cin >> cities[i].second;
    }
    double junction_cost; cin >> junction_cost;
    double failure_probability; cin >> failure_probability;

    // construct a graph
    typedef adjacency_matrix< undirectedS, no_property, property < edge_weight_t, double > > Graph;
    typedef graph_traits< Graph >::edge_descriptor Edge;
    Graph g(NC + 1);
    auto distance = [](pair<int, int> const & a, pair<int, int> const & b) {
        return hypot(b.first - a.first, b.second - a.second);
    };
    REP (i, NC) REP (j, i) {
        add_edge(i, j, distance(cities[i], cities[j]), g);
    }

    auto kruskal = [&]() {
        vector< Edge > spanning_tree;
        kruskal_minimum_spanning_tree(g, back_inserter(spanning_tree));
        property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
        double acc = 0;
        for (auto e : spanning_tree) {
            acc += weight[e];
        }
        return acc;
    };

    // reference value
    add_edge(0, NC, 0, g);
    cout << kruskal() << endl;
    remove_edge(0, NC, g);

    // generate data
    int N = (getenv("N") ? atoi(getenv("N")) : 100);
    vector<int> values;
    for (int z = 0; z < S; z += max(1, S / N)) {
        values.push_back(z);
    }
    values.push_back(S);
    for (int y : values) {
        for (int x : values) {
            REP (i, NC) {
                add_edge(i, NC, distance(cities[i], make_pair(x, y)), g);
            }
            double value = kruskal();
            cout << x << ' ' << y << ' ' << value << endl;
            REP (i, NC) {
                remove_edge(i, NC, g);
            }
        }
    }
    return 0;
}
