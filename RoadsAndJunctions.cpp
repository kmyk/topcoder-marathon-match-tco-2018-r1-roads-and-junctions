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

class RoadsAndJunctions {
public:
    int NC;
    vector<int> buildJunctions(int S, vector<int> cities, double junctionCost, double failureProbability) {
        // store number of cities for building the roads
        NC = cities.size() / 2;
        // build one junction at (0,0)
        return vector<int>({0, 0});
    }
    vector<int> buildRoads(vector <int> junctionStatus) {
        // build a road from the single junction to each city
        // (we assume that it was built, but don't check it)
        vector<int> ret(2 * NC);
        for (int i = 0; i < NC; ++i) {
            ret[2 * i] = NC;
            ret[2 * i + 1] = i;
        }
        return ret;
    }
};
