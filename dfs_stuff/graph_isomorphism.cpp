/*
 *  Graph isomorphism by backtracking
 *  Runs quite well on graphs with small degree (up to 50 nodes for random 4-regular graphs).
 */

bool graph_iso(vector<vector<int> > const&g1, vector<vector<int> > const& g2, vector<int> &code, vector<int> const&ord, int d){
    if(d == ord.size()) return true;
    int v = ord[d];

    // check consistency
    auto check = [&](){
        int v2 = code[v];
        // duplicate
        if(count(code.begin(), code.end(), v2)>1) return false;
        if(g1[v].size() != g2[v2].size()) return false;
        for(auto &e:g1[v]){
            int w2 = code[e];
            if(w2!=-1){
                if(find(g2[v2].begin(), g2[v2].end(), w2) == g2[v2].end()){
                    return false;
                }
            }
        }
        return true;
    };
    int n = g1.size();
    for(int i=0;i<n;++i){
        code[v] = i;
        if(check()){
            if(graph_iso(g1, g2, code, ord, d+1)){
                return true;
            }
        }
    }
    code[v] = -1;
    return false;
}

vector<int> get_isomo(vector<vector<int> > const&g1, vector<vector<int> > const&g2){
    int n = g1.size();
    vector<int> ord, code(n, -1);
    queue<int> s;
    s.emplace(uniform_int_distribution<int>(0, n-1)(rng));
    vector<char> vis(n);
    while(!s.empty()){
        int u=s.front();
        s.pop();
        if(vis[u]) continue;
        vis[u]=1;
        ord.push_back(u);
        for(auto &e:g1[u]) s.emplace(e);
    }
    if(ord.size()!=n) {
        cerr << "not connected?\n";
        return vector<int>();
    }
    // random order, probably worse than the BFS-order from above.
    //iota(ord.begin(), ord.end(), 0);
    //shuffle(ord.begin(), ord.end(), rng);
    bool did = graph_iso(g1, g2, code, ord, 0);
    if(!did){
        return vector<int>();
    } else {
        return code;
    }
}
