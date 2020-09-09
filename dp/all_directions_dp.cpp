// All directions DP in O(n)
template<typename Traits>
class All_Directions_DP{
public:
    using Val = typename Traits::Val;

    All_Directions_DP(int n_) : n(n_), d(n), up(n), down(n), vertex(n){
        pre.reserve(n);
        suff.reserve(n);
    }
    void run(vector<vector<int> > g_){
        g = move(g_);
        dfs_1<true>(0, -1);
        dfs_2<true>(0, -1);
    }

    Val const& get_edge_val(int a, int b){
        return d[a] < d[b] ? down[a] : up[b];
    }
    Val const& get_vertex_val(int u){
        return vertex[u];
    }

private:
    template<bool is_root>
    void dfs_1(int u, int p){
        if(!is_root){
            up[u] = Traits::init(u, p);
            d[u] = d[p]+1;
            rotate(g[u].begin(), find(g[u].begin(), g[u].end(), p), g[u].end());
        } else {
            d[u] = 0;
        }
        for(auto &e:g[u]) if(e != p){
            dfs_1<false>(e, u);
            if(!is_root){
                up[u] = Traits::merge(up[u], up[e]);
            }
        }
        if(!is_root){
            Traits::post_edge(up[u], u, p);
        }
    }
    template<bool is_root>
    void dfs_2(int u, int p){
        const int k = g[u].size();
        pre.assign(k+1, Traits::neutral());
        suff.assign(k+1, Traits::neutral());
        if(!is_root){
            pre[1] = suff[0] = down[u];
        }
        for(int i=!is_root; i<k; ++i){
            auto const e = g[u][i];
            pre[i+1] = suff[i] = up[e];
        }
        for(int i=1; i<=k; ++i){
            pre[i] = Traits::merge(pre[i-1], pre[i]);
        }
        vertex[u] = pre.back();
        Traits::post_vertex(vertex[u], u);
        for(int i=k-1; i>=0; --i){
            suff[i] = Traits::merge(suff[i], suff[i+1]);
        }
        for(int i=!is_root; i<k; ++i){
            auto const e = g[u][i];
            down[e] = Traits::merge(Traits::init(u, e), Traits::merge(pre[i], suff[i+1]));
            Traits::post_edge(down[e], u, e);
        }
        for(auto const&e:g[u]) if(e != p){
            dfs_2<false>(e, u);
        }
    }

    int n;
    vector<vector<int> > g;
    vector<int> d;
    vector<Val> up, down, pre, suff, vertex;
};

// Example: Sum of distances for each vertex
struct AD_Traits{
    struct Val{
        int siz;
        ll total;
    };
    ///using Val = int;

    static Val neutral(){
        return Val{0, 0};
    }
    // edge a -> b
    static Val init(int a, int b){
        return Val{1, 0};
    }
    static Val merge(Val const&a, Val const&b){
        return Val{a.siz+b.siz, a.total+b.total};
    }
    static void post_edge(Val &v, int a, int b){
        v.total += v.siz;
    }
    static void post_vertex(Val &v, int a){   
    }
};
