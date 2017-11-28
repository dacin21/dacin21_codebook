/*
 *  directed strongly connected components in O(N+M)
 *
 */

#include <bits/stdc++.h>
using namespace std;
struct StrongComp{
    int N;
    vector<vector<int> > G, invG;
    vector<int> component;
    int number_of_components;
    vector<bool> vis;
    stack<int> s;
    StrongComp(int _N):N(_N), G(N), invG(N){}
    void add_edge(int a, int b){
        assert(a>=0 && a<N);
        assert(b>=0 && b<N);
        G[a].push_back(b);
        invG[b].push_back(a);
    }
    void dfs(int cur){
        vis[cur]=true;
        for(int to:G[cur])
            if(!vis[to])
                dfs(to);
        s.push(cur);
    }
    void invDfs(int cur){
        vis[cur]=false;
        component[cur] = number_of_components;
        for(int to:invG[cur])
            if(vis[to])
                invDfs(to);
    }
    void calc_components(){
        vis.clear(); vis.resize(N, false);
        component.clear();component.resize(N);
        number_of_components = 0;
        for(int i=0;i<N;++i)
            if(!vis[i])
                dfs(i);
        for(;!s.empty();s.pop()){
            if(vis[s.top()]){
                invDfs(s.top());
                ++number_of_components;
        }   }
    }
};
