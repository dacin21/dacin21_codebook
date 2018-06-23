#include <bits/stdc++.h>
using namespace std;
// Push relabel matching in O(V^0.5 * E)
template<typename flow_t = int>
struct PushRelabelMatching {
    struct Edge {
        int to, rev;
    };
    int n, m, hi, tim;
    vector<vector<Edge> > g;
    vector<Edge*> cur;
    vector<vector<int> > hs;
    vector<int> H, match;
    queue<int> q;
    PushRelabelMatching(int N, int M) : n(N+1), m(n+M+1), g(m), cur(n), hs(m+1), H(m), match(m, -1){}
    int add_edge(int s, int t) {
        t+=n;
        assert(0<=s && s<n && n<=t && t<m);
        Edge a = {t, (int)g[t].size()};
        Edge b = {s, (int)g[s].size()};
        g[s].push_back(a);
        g[t].push_back(b);
        return b.rev;
    }
    void enq(int u){
        hs[H[u]].push_back(u);
    }
    flow_t max_matching() {
        int v = g.size();
        for(int i=0;i<n;++i) hs[H[i]=1].push_back(i);
        for(int i=0;i<n;++i) cur[i] = g[i].data();
        if(n)
        for (hi = 1, tim=0;hi<m;) {
            int u = hs[hi].back();
            hs[hi].pop_back();
            for(;;)
                if (cur[u] == g[u].data() + g[u].size()) {
                    H[u] = 1e9;
                    for(auto &e:g[u])
                        if (H[u] > H[e.to]+1){
                            H[u] = H[e.to]+1;
                            cur[u] = &e;
                        }
                    if(H[u]>=m)break;
                    if(tim++==v){
                        tim = hi =0;
                        fill(match.begin(), match.begin()+n, -1);
                        fill(H.begin(), H.end(), m);
                        for(int i=n;i<m;++i)
                            if(~match[i]) match[match[i]] = i;
                            else {q.push(i); H[i]=0;}
                        for(auto &e:hs) e.clear();
                        while(!q.empty()){
                            int a=q.front(); q.pop();
                            for(auto const&e:g[a])
                            if(H[e.to] == m){
                                H[e.to] = H[a]+1;
                                if(~match[e.to]){
                                    H[match[e.to]] = H[a]+2;
                                    q.emplace(match[e.to]);
                                } else enq(e.to);
                            }
                        }
                        break;
                    }
                } else if (H[u] == H[cur[u]->to]+1){
                    if(~match[cur[u]->to]) enq(match[cur[u]->to]);
                    match[cur[u]->to] = u;
                    H[cur[u]->to]+= 2;
                    break;
                } else ++cur[u];
            if(hi) --hi; if(hi) --hi;
            while (hi<m && hs[hi].empty()) ++hi;
        }
        int ans = 0;
        for(int i=n;i<m;++i)
            if(~match[i]){
                ++ans;
                match[match[i]] = i;
            }
        return ans;
    }
};
namespace FIO{
    char buf[32*1042|1];
    int bc=0, be=0;
    char gc(){
        if(bc>=be){
            be = fread(buf, 1, sizeof(buf)-1, stdin);
            buf[be] = bc = 0;
        }
        return buf[bc++];
    }
    void readint(){}
    void readuint(){}
    template<typename T, typename ...S>
    void readuint(T &a, S& ...b){
        a=0;
        int x=gc();
        while(x<'0' || x>'9') x=gc();
        while(x>='0' && x<='9'){
            a = a*10+x-'0'; x=gc();
        }
        readuint(b...);
    }
    template<typename T, typename ...S>
    void readint(T &a, S& ...b){
        a=0;
        int x=gc(), s=1;;
        while(x!='-' && (x<'0' || x>'9')) x=gc();
        if(x=='-'){ s=-s; x=gc(); }
        while(x>='0' && x<='9'){
            a = a*10+x-'0'; x=gc();
        }
        if(s<0) a=-a;
        readint(b...);
    }
}
using FIO::readint;
using FIO::readuint;
signed matching(){
    int N, M, E;
    readuint(N, M, E);
    //scanf("%d%d%d", &N, &M, &E);
    PushRelabelMatching<> fl(N+1, M+1);
    int u, v;
    for (int i = 0; i < E; ++i) {
        readuint(u, v);
        //scanf("%d%d", &u, &v);
        fl.add_edge(u, v);
    }
    int ans = fl.max_matching();
    printf("%d\n",ans);
    return 0;
}

signed main() {
    #ifdef LOCAL_RUN
    freopen("in.txt", "r", stdin);
    #endif // LOCAL_RUN
    return matching();
}

