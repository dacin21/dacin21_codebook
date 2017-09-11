template<typename T = int, typename comp = less<T>>
vector<int> offline_RMQ(vector<pair<int, int> > const&qs, vector<T> const&v){
    int n=v.size(), qq=qs.size();
    vector<int> p(n, -1), ans(qq);
    auto f = [&](int x){return ~p[x] ? p[x]=f(p[x]):x;};
    vector<vector<pair<int, int> > > q(n);
    for(int i=0;i<qq;++i) q[qs[i].first].emplace_back(qs[i].second-1, i);
    stack<int> s;
    for(int i=0;i<n;++i){
        for(;!s.empty()&& comp()(v[i], v[s.top()]);s.pop()) p[s.top()]=i;
        s.push(i);
        for(auto &e:q[i]) ans[e.second] = f(e.first);
    }
    return ans;
}
