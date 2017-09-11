// bridges and articulation point in O(N+M)
#include <bits/stdc++.h>
using namespace std;
template<typename T1, typename T2>
void xmin(T1 &a, T2 const&b){ if(a>b) a=b;}
struct Artpoints{
	int n, m;
	vector<vector<pair<int, int> > > g;
	vector<char> isap, isbr;
	vector<array<int, 3> > brid;
	Artpoints(int _n):n(_n), g(n+1){}
	void add_edge(int a, int b){
		g[a].emplace_back(b, m);
		g[b].emplace_back(a, m);
		++m;
	}
	void calc(){
		isap.assign(n+1, 0);
		isbr.assign(m, 0);
		vector<int> d(n+1, -1), l(n+1), p(n+1);
		vector<pair<int, int>*>cur(n+1);
		for(int i=0;i<=n;++i) cur[i] = g[i].data();
		stack<int> s;
		for(int i=0;i<=n;++i){
			if(d[i]<0){
				p[i] = -1;
				d[i] = 0;
				int rchild = 0;
				s.push(i);
				while(!s.empty()){
					int a = s.top();
					if(cur[a]>g[a].data()){
					    xmin(l[a], l[cur[a][-1].first]);
					    if(l[cur[a][-1].first]>=d[a]) isap[a]=1;
					}
					for(auto&x = cur[a]; x<g[a].data()+g[a].size() && d[x->first]>=0;++x){
						if(x->second != p[a]){
							xmin(l[a], d[x->first]);
						}
					}
					if(cur[a]<g[a].data()+g[a].size()){
						if(a==i) ++rchild;
						auto&x = cur[a];
						l[x->first] = d[x->first] = d[a]+1;
						p[x->first] = x->second;
						s.push(x->first);
						++cur[a];
					} else {
						if(p[a]>=0 && l[a]==d[a]){
							isbr[p[a]] = 1;
						}
						s.pop();
					}
				}
				isap[i] = (rchild>1);
			}
		}
		for(int i=0;i<=n;++i){
			for(auto const&e:g[i]){
				if(e.first<i && isbr[e.second]){
					bridges.push_back({e.first, i, e.second});
				}
			}
		}
	}
};