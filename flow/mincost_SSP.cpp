// successive shortest paths algorithm for min cost flow
// runs in O(f * E log V), which is good for small flows f.


// fast priority queue for integers
// source: https://judge.yosupo.jp/submission/8774
template <class Key, class T>
struct Radix_Heap{
    static_assert(is_integral<Key>::value, "only integers supported");
    static int lg(Key a) {
        return a ? __lg(a) : -1;
    }

    void emplace(Key key, T val) {
        assert(key >= last);
        v[lg(key ^ last) + 1].emplace_back(key, val);
        ++sz;
    }
    void pull() {
        if (ptr < v[0].size()) return;
        v[0].clear(), ptr = 0;
        int i = 1;
        while (v[i].empty()) ++i;
        last = min_element(begin(v[i]), end(v[i]), [](auto a, auto b) {
            return a.first < b.first;
        })->first;
        for (auto e : v[i]) v[lg(e.first ^ last) + 1].push_back(e);
        v[i].clear();
    }
    pair<Key, T> top(){
        pull();
        return v[0][ptr];
    }
    void pop(){
        pull();
        --sz;
        ++ptr;
    }
    bool empty() const {
        return sz == 0;
    }
    void clear(){
        sz = 0;
        ptr = 0;
        last = 0;
        for(auto &e:v) e.clear();
    }

    array<vector<pair<Key, T>>, sizeof(Key) * 8 + 1> v;
    Key last = 0;
    size_t sz = 0, ptr = 0;
};


typedef int flow_t;
typedef long long cost_t;
struct mcFlow{
    struct Edge{
        cost_t c;
        flow_t f;
        int to, rev;
        Edge(int _to, cost_t _c, flow_t _f, int _rev) : c(_c), f(_f), to(_to), rev(_rev) {}
    };
    const cost_t INFCOST = numeric_limits<cost_t>::max()/2;
    const cost_t INFFLOW = numeric_limits<flow_t>::max()/2;

    mcFlow(int _N, int _S, int _T):N(_N), S(_S), T(_T), G(_N) {}
    void add_edge(int a, int b, cost_t cost, flow_t cap, flow_t rcap=0){
        if(a==b){
            assert(cost >= 0 || (cap == 0 && rcap == 0));
            return;
        }
        assert(a>=0&&a<N && b>=0&&b<N);
        G[a].emplace_back(b, cost, cap, G[b].size());
        G[b].emplace_back(a, -cost, rcap, G[a].size()-1);
    }

    template<bool require_max_flow>
    pair<flow_t, cost_t> minCostFlow(vector<cost_t> phi = vector<cost_t>(N)){
        assert((int)phi.size() == N);
        vector<cost_t> dist(N);
        vector<int> state(N);
        vector<Edge*> from(N, 0);

        Radix_Heap<cost_t, int> q;
        cost_t retCost=0;
        cost_t bestCost=0;
        flow_t retFlow=0;
        do {
            fill(dist.begin(), dist.end(), INFCOST);
            fill(state.begin(), state.end(), 0);
            fill(from.begin(), from.end(), (Edge*)0);
            q.clear();
            dist[S]=0; state[S]=1;
            q.emplace(0, S);
            while(!q.empty()){
                int cur;
                cur = q.top().second;q.pop();
                if(state[cur] == 2) continue;
                state[cur]=2;
                for(Edge &e:G[cur]){
                    if(e.f==0) continue;
                    cost_t newDist = dist[cur] + phi[cur] - phi[e.to] + e.c;
                    if(newDist < dist[e.to]){
                        dist[e.to]=newDist; from[e.to]=&e;
                        q.emplace(newDist, e.to);
                        state[e.to]=1;
                    }
                }
            }
            if(from[T]==0) break;
            flow_t augment=INFFLOW;
            for(Edge*e = from[T];e;e=from[G[e->to][e->rev].to]){
                augment=min(augment, e->f);
            }
            for(Edge*e = from[T];e;e=from[G[e->to][e->rev].to]){
                retCost+=e->c*augment;
                e->f-=augment;
                G[e->to][e->rev].f+=augment;
            }
            retFlow+=augment;
            for(int i=0;i<N;++i){
                phi[i]+=dist[i];
            }
            if(!require_max_flow && bestCost <= retCost) break;
            bestCost = retCost;
        } while(from[T]);
        return make_pair(retFlow, bestCost);
    }
    flow_t getFlow(Edge const&e){
        return G[e.to][e.rev].f;
    }

    int N, S, T;
    vector<vector<Edge> > G;
};
