/*
 *  Highest label push-relabel algorithm in O(V^2 E^0.5)
 *  Uses various heuristics
 *  Code by me (dacin21)
 *
 *  Special thanks to min_25 and Chilli for showing me a faster version of gap-relabeling.
 *
 */

#include <bits/stdc++.h>
using namespace std;

template<typename cap_t, typename excess_t, bool global_relabeling = true, bool min_cut_only = false, bool shuffle_edges = false>
class Push_Relabel{
public:
    struct Edge{
        int to, rev;
        cap_t f;
    };

    Push_Relabel(int n_):n(n_), m(0){}

    void add_edge(int u, int v, cap_t c, cap_t c_rev = 0){
        edge_pool.emplace_back(u, v, c, c_rev);
        ++m;
    }
    excess_t max_flow(int s_, int t_){
        s = s_; t = t_;
        run_pr();
        return excess[t]-1;
    }

private:
    void compile_g(){
        g_pos.resize(n+1);
        if(shuffle_edges) random_shuffle(edge_pool.begin(), edge_pool.end());
        for(auto &e:edge_pool){
            ++g_pos[get<0>(e)];
            ++g_pos[get<1>(e)];
        }
        for(int i=0;i<n;++i){
            g_pos[i+1]+=g_pos[i];
        }
        g.resize(g_pos.back());
        for(auto &e:edge_pool){
            int u, v; cap_t c, c_rev;
            tie(u, v, c, c_rev) = e;
            const int i = --g_pos[u], j = --g_pos[v];
            g[i] = Edge{v, j, c};
            g[j] = Edge{u, i, c_rev};
        }
    }
    void global_relabel(){
        q.reserve(n);
        fill(h.begin(), h.end(), n);
        fill(h_cnt.begin(), h_cnt.end(), 0);
        h_cnt[n] = 1;
        q.push_back(t);
        h[t] = 0;
        for(auto &e:buck) e.clear();
        for(auto &e:buck_all) e.clear();
        for(auto it = q.begin();it<q.end();++it){
            const int u = *it;
            if(u != t && excess[u]){
                hi = h[u];
                buck[h[u]].push_back(u);
            }
            if(u != t) buck_all[h[u]].push_back(u);
            ++h_cnt[h[u]];
            for(int i = g_pos[u],i_end = g_pos[u+1];i < i_end;++i){
                Edge const&e = g[i];
                if(g[e.rev].f && h[e.to] == n){
                    h[e.to] = h[u]+1;
                    q.push_back(e.to);
                }
            }
        }
        hi_all = h[q.back()];
        assert(h[s] == n);
        q.clear();
    }
    void push(int u, Edge &e, excess_t f){
        if(!excess[e.to]){
            buck[h[e.to]].push_back(e.to);
        }
        Edge&back = g[e.rev];
        e.f-=f;
        back.f+=f;
        excess[e.to]+=f;
        excess[u]-=f;
    }
    void init_pr(){
        compile_g();
        cur.assign(n, 0);
        copy(g_pos.begin(), prev(g_pos.end()), cur.begin());
        h.resize(n);
        excess.assign(n, 0);
        buck.resize(2*n);
        buck_all.resize(n+1);
        h_cnt.assign(2*n, 0);
        h[s] = n;
        h_cnt[n] = 1;
        h_cnt[0] = n-1;
        excess[t] = 1;
    }
    void run_pr(){
        init_pr();
        for(int i = g_pos[s],i_end = g_pos[s+1];i < i_end;++i){
            push(s, g[i], g[i].f);
        }
        hi = hi_all = 0;
        if(global_relabeling) global_relabel();
        if(!buck[hi].empty())
        for(;hi>=0;){
            int u = buck[hi].back(); buck[hi].pop_back();
            int u_cur = cur[u];
            //discharge
            if(!min_cut_only || h[u] < n)
            while(excess[u] > 0){
                if(__builtin_expect(u_cur == g_pos[u+1], false)){
                    int new_h = 1e9;
                    for(int i = g_pos[u],i_end = g_pos[u+1];i < i_end;++i){
                        auto const&e = g[i];
                        if(e.f && new_h > h[e.to]+1){
                            new_h = h[e.to]+1;
                            u_cur = i;
                        }
                    }
                    ++h_cnt[new_h];
                    h[u] = new_h;
                    if(__builtin_expect(!--h_cnt[hi] && hi < n, false)){
                        // gap relabel
                        for(int j = hi;j <= hi_all;++j){
                            for(auto &f:buck_all[j]) if(!min_cut_only || h[f] < n){
                                --h_cnt[h[f]];
                                h[f] = n+1;
                            }
                            buck_all[j].clear();
                        }
                    }
                    hi = h[u];
                } else {
                    Edge &e_cur = g[u_cur];
                    if(e_cur.f && h[u] == h[e_cur.to]+1){
                        push(u, e_cur, min<excess_t>(excess[u], e_cur.f));
                    } else ++u_cur;
                }
            }
            if(h[u] < n) {
                hi_all = max(hi_all, h[u]);
                buck_all[h[u]].push_back(u);
            }
            cur[u] = u_cur;
            while(hi>=0 && buck[hi].empty()) --hi;
        }
    }

    int n, m, s, t, hi, hi_all;
    vector<tuple<int, int, cap_t, cap_t> > edge_pool;
    vector<int> g_pos;
    vector<Edge> g;
    vector<int> q, cur, h, h_cnt;
    vector<excess_t> excess;
    vector<vector<int> > buck, buck_all;
};
