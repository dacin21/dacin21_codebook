/*
 *  computes an optimal matrix chain-multiplication in O(N^2)
 *  Algorithm due to F. FRANCES YAO, SPEED-UP IN DYNAMIC PROGRAMMING
 */

 #include <bits/stdc++.h>
 using namespace std;


long long matrix_chain_multiplication(vector<int> const&v){
    int N=v.size();
    vector<pair<int, int> > m(N);
    for(int i=0;i<N;++i)m[i] = make_pair(v[i], i);
    rotate(m.begin(), min_element(m.begin(), m.end()), m.end());
    for(int i=0;i<N;++i)m[i].second=i;
    m.push_back(m[0]);
    vector<pair<int, int> > br;
    map<pair<int, int>, pair<int, int> > childs;
    // find bridges
    stack<pair<int, int> > s;
    for(int i=1;i<=N;++i){
        s.push(m[i-1]);
        while(s.top()>m[i]){
            pair<int, int> t = s.top(); s.pop();
            br.emplace_back(s.top().second, t.second);
            br.emplace_back(t.second, m[i].second);
            childs[make_pair(s.top().second, m[i].second)] = make_pair(br.size()-1, br.size()-2);
        }
    }
    vector<vector<long long> > DP(br.size(), vector<long long>(N+1, -1));
    for(int ind=0;ind<br.size();++ind){
        auto e = br[ind];
        int i=e.first, j=e.second;
        if(m[i]>m[j]) swap(i, j);
        if(!childs.count(e)){ //leaf
            DP[ind][i] = 0;//line
            for(int k=0;k<N;++k)
                if(m[k]<m[i]) DP[ind][k] = m[i].first*(long long)m[j].first*m[k].first;
        } else {
            int s1 =childs[e].first, s2 = childs[e].second;
            DP[ind][i] = DP[s1][i]+DP[s2][i];
            for(int k=0;k<N;++k){
                if(m[k]<m[i]) DP[ind][k] = min(DP[ind][i]+m[i].first*(long long)m[j].first*m[k].first,DP[s1][k]+DP[s2][k]);
            }
        }
    }
    long long ans = DP[childs[make_pair(0, 0)].first][0]+DP[childs[make_pair(0, 0)].second][0];
    return ans;
}