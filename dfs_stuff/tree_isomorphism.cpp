/*
 *  Tree isomorphism by sorting in O(N * log N)
 */
void calc_iso(vector<vector<vector<int> > > const&G, vector<pair<int, int> > const&C, vector<int> &ret){
    int N=G.size();
    vector<vector<int> > pre(2*N);
    vector<vector<vector<int> > > sig(2*N);
    for(int i=0;i<N;++i){
        pre[i].resize(G[i].size(), -1);
        sig[i].resize(G[i].size());
        if(C[i].first!=C[i].second){
            pre[i+N].resize(G[i].size(), -1);
            sig[i+N].resize(G[i].size());
        }
    }
    vector<vector<pair<int, int> > > lay(1);
    for(int i=0;i<2*N;++i){
        int j= (i<N?i:i-N);
        if(i>=N && C[j].first==C[j].second) continue;
        int d=0;
        queue<int> q;
        q.emplace(i<N?C[j].first:C[j].second);
        q.emplace(-1);
        while(q.size()>1){
            int a=q.front();q.pop();
            if(a==-1){
                ++d;
                if(d==(int)lay.size()) lay.emplace_back();
                q.emplace(a);
            } else {
                lay[d].emplace_back(i, a);
                for(int e:G[j][a]){
                    if(e!=pre[i][a]){
                        pre[i][e] = a;
                        q.emplace(e);
        }   }   }   }
    }
    for(int i=(int)lay.size()-1;i>=0;--i){
        sort(lay[i].begin(), lay[i].end(), [&sig, N](pair<int, int> const&a, pair<int, int> const&b){return sig[a.first][a.second] < sig[b.first][b.second];});
        int cur=0;
        vector<int> lastSig;
        for(int j=0;j<(int)lay[i].size();++j){
            vector<int> &si = sig[lay[i][j].first][lay[i][j].second];
            if(j&& si!= lastSig){
                ++cur;
            }
            lastSig.swap(si);
            vector<int>(1, cur).swap(si);
            if(pre[lay[i][j].first][lay[i][j].second]!=-1){
                sig[lay[i][j].first][pre[lay[i][j].first][lay[i][j].second]].emplace_back(si[0]);
            }
        }
    }
    for(int i=0;i<N;++i){
        if(C[i].first==C[i].second){
            ret[i] = sig[i][C[i].first][0];
        } else{
            ret[i] = min(sig[i][C[i].first][0], sig[i+N][C[i].second][0]);
        }
    }
}
pair<int, int> get_center(vector<vector<int> > const&G){
    queue<int> q;
    vector<int> pre(G.size(), -1);
    q.emplace(0);
    int a;
    while(!q.empty()){
        a=q.front();q.pop();
        for(int const&e:G[a]){
            if(e!=pre[a]){
                q.emplace(e);
                pre[e]=a;
    }   }   }
    q.emplace(a);
    fill(pre.begin(), pre.end(), -1);
    while(!q.empty()){
        a=q.front();q.pop();
        for(int const&e:G[a]){
            if(e!=pre[a]){
                q.emplace(e);
                pre[e]=a;
    }   }   }
    int len=0;
    for(int i=a;i!=-1;i=pre[i])++len;
    for(int i=1;i<(len+1)/2;++i)a=pre[a];
    if(len%2) return make_pair(a, a);
    return make_pair(a, pre[a]);
}
vector<int> label(vector<vector<vector<int> > > const&G){
    int N=G.size();
    vector<int> labels(N);
    vector<pair<int, int> > centers(N);
    for(int i=0;i<N;++i){
        centers[i] = get_center(G[i]);
    }
    calc_iso(G, centers, labels);
    return labels;
}
