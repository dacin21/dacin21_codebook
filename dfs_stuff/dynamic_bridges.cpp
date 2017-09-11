/*
 *  Off-line dynamic bridges
 *  add/remove edges and report number of bridges
 *  runs in O(N + Q log Q)
 *  CF: (0.75 s for Q ~1e5)
 */
#include <bits/stdc++.h>
using namespace std;

struct Bicomp{
	int n;
    vector<vector<pair<int, int> > > g;
    vector<int> comp, comp2, add;
    int compcnt, compcnt2;
    vector<array<int, 3> > bridges;
    vector<char> used, used2, vis;
    stack<int> s;
    Bicomp(int _n):n(_n), g(n), comp(n), comp2(n, -1), add(n, 0), compcnt(0), compcnt2(0), used(n, 0), used2(n, 0), vis(n, 0){}
    void add_edge(int a, int b){
        g[a].push_back({b, 1});
        g[b].push_back({a, 1});
        used[a]=used[b]=1;
    }
    void add_bridge_edge(int a, int b, int w){
        g[a].push_back({b, w});
        g[b].push_back({a, w});
    }
    int rec(int c){
        vis[c]=1;
        s.push(c);
        for(auto const&e:g[c]){
            if(vis[e.first]==2) continue;
            else if(vis[e.first]==1){
                --add[c];
                ++add[e.first];
            } else {
                int x = rec(e.first);
                add[c]-=x;
                if(x==1){
                    do{
                        x = s.top();
                        comp[x] = compcnt;
                        s.pop();
                    } while(x!=e.first);
                    ++compcnt;
                    bridges.push_back({c, e.first, e.second});
                }
            }
        }
        vis[c]=2;
        return -add[c];
    }
    pair<int, int> rec2(int c){
        vis[c]=1;
        pair<int, int> ret={-1, 0};
        if(used2[c]) ret.first = comp2[c] =  compcnt2++;
        for(auto const&e:g[c]){
            if(!vis[e.first]){
                auto a = rec2(e.first);
                a.second+=e.second;
                if(a.first>=0){
                    if(ret.first<0) ret=a;
                    else {
                        if(comp2[c]<0){
                            comp2[c] = compcnt2++;
                            bridges.push_back({comp2[c], ret.first, ret.second});
                            ret = {comp2[c], 0};
                        }
                        bridges.push_back({comp2[c], a.first, a.second});
                    }
                }
            }
        }
        return ret;
    }
    pair<size_t, size_t> calc(){
        size_t br = 0, br2 = 0;
        // compute bridges
        for(int i=0;i<n;++i){
            if(vis[i]==0){
                assert(rec(i)==0);
                while(!s.empty()){
                    comp[s.top()] = compcnt;
                    s.pop();
                }
                ++compcnt;
            }
        }
        // compress forest
        fill(vis.data(), vis.data()+compcnt, 0);
        for(int i=0;i<compcnt;++i) g[i].clear();
        for(int i=0;i<n;++i) used2[comp[i]]|=used[i];
        for(auto &e:bridges){
            br+=e[2];
            add_bridge_edge(comp[e[0]], comp[e[1]], e[2]);
        }
        bridges.clear();
        for(int i=0;i<compcnt;++i)
            if(!vis[i]) rec2(i);
        for(auto &e:bridges) br2+=e[2];
        return {br2, br};
    }
	void reset(int _n){
	    n=_n;
		assert(n<=(int) g.size());
		for(int i=0;i<n;++i) g[i].clear();
		fill(comp.begin(), comp.begin()+n, -1);
		fill(comp2.begin(), comp2.begin()+n, -1);
		compcnt = compcnt2 = 0;
		bridges.clear();
		fill(add.begin(), add.begin()+n, 0);
		fill(used.begin(), used.begin()+n, 0);
		fill(used2.begin(), used2.begin()+n, 0);
		fill(vis.begin(), vis.begin()+n, 0);
	}
};

struct Bidynacon{
    int N;
    vector<int> ans;
    vector<array<int, 3> > ups;
    vector<int> qs;
    long long ops = 0, nods=0;
	Bicomp com;
    Bidynacon(int _N):N(_N), com(N){}
    void add_edge(int a, int b){
        if(a>b) swap(a, b);
        assert(a>=0 && b<N);
        ups.push_back({1, a, b});
    }
    void rem_edge(int a, int b){
        if(a>b) swap(a, b);
        assert(a>=0 && b<N);
        ups.push_back({-1, a, b});
    }
    void add_query(){
        qs.push_back(ups.size());
    }
    void rec(int l, int r, int n, int k, vector<array<int, 4> > const&eds, vector<array<int, 3> > const&br){
        // remove pruning if queries are dense
        if(accumulate(ans.begin()+l, ans.begin()+r+1, 0)==0) return;
        com.reset(n);
        vector<array<int, 4> > eds2;
        eds2.reserve(eds.size());
        for(auto const&e:eds){
            if(e[1]<l || e[0]>=r) continue;
            if(e[2]==e[3]) continue;
            if(e[0]<l && e[1]>=r){
                com.add_edge(e[2], e[3]);
            } else {
                eds2.push_back(e);
                com.used[e[2]] = 1;
                com.used[e[3]] = 1;
            }
        }
        for(auto const&e:br) com.add_bridge_edge(e[0], e[1], e[2]);
        auto e = com.calc();
        n = com.compcnt2;
        k+=e.second-e.first;
        auto br2 = com.bridges;
        for(auto &e:eds2){
            e[2] = com.comp2[com.comp[e[2]]];
            e[3] = com.comp2[com.comp[e[3]]];
        }
        if(eds2.empty()){
            for(int i=l;i<=r;++i){
                ans[i] = k+e.first;
            }
            return;
        }
        rec(l, l+(r-l)/2, n, k, eds2, br2);
        rec(l+(r-l)/2+1, r, n, k, eds2, br2);
    }
    vector<int> calc(){
        ops=nods=0;
        vector<array<int, 4> > edges;
        map<pair<int, int>, vector<int> > cur;
        int Q = ups.size();
        for(int i=0;i<Q;++i){
            auto const&e=ups[i];
            if(e[0]==1){
                cur[{e[1], e[2]}].push_back(i);
            } else {
                auto &f = cur[{e[1], e[2]}];
                assert(!f.empty());
                edges.push_back({f.back(), i, e[1], e[2]});
                f.pop_back();
            }
        }
        for(auto &e:cur){
            for(auto &f:e.second){
                edges.push_back({f, 1+Q, e.first.first, e.first.second});
            }
        }
        ans.resize(Q+1);
        for(auto e:qs) ans[e]=1;
        rec(0, Q, N, 0, edges, vector<array<int, 3>>());
        vector<int> ret(qs.size());
        for(size_t i=0;i<qs.size();++i){
            ret[i] = ans[qs[i]];
        }
        return ret;
    }
};