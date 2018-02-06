// Heavy light decomposition
// Tree path updates and queries in O(log(n)^2)
template<class Segtree_Data>
class Heavy_Light{
private:
    struct Segtree{
        int N, height;
        vector<typename Segtree_Data::node_t> data;
        vector<typename Segtree_Data::update_t> lazy;
        Segtree(vector<typename Segtree_Data::node_t> const&base):N(base.size()), height(__builtin_clz(1)-__builtin_clz(N)), data(2*N, Segtree_Data::node_ne()), lazy(2*N, Segtree_Data::update_ne()){
            copy(base.begin(), base.end(), data.begin()+N);
            for(int i=N-1;i>=0;--i)
                data[i]=Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
        }
        Segtree(int n):N(n), height(__builtin_clz(1)-__builtin_clz(N)), data(2*N, Segtree_Data::node_init()), lazy(2*N, Segtree_Data::update_ne()){
            for(int i=N-1;i>=0;--i)
                data[i]=Segtree_Data::merge_nodes(data[i<<1], data[i<<1|1]);
        }
        void localUpdate(int pos, typename Segtree_Data::update_t const&val){
            Segtree_Data::update_node(data[pos], val);
            if(pos<N) lazy[pos] = Segtree_Data::merge_lazy(lazy[pos], val);
        }
        void push(int pos){
            for(int s=height;s>0;--s){
                int i=pos>>s;
                if(lazy[i]!=Segtree_Data::update_ne()){
                    localUpdate(i<<1, lazy[i]);
                    localUpdate(i<<1|1, lazy[i]);
                    lazy[i]=Segtree_Data::update_ne();
                }
            }
        }
        void rePath(int pos){
            while(pos>>=1) Segtree_Data::update_node(data[pos] = Segtree_Data::merge_nodes(data[pos<<1], data[pos<<1|1]), lazy[pos]);
        }
        void update(int l, int r, typename Segtree_Data::update_t const&val){
            int l2=l+=N, r2=r+=N;
            push(l2); push(r2-1);
            for(;l<r;l>>=1, r>>=1){
                if(l&1) localUpdate(l++, val);
                if(r&1) localUpdate(--r, val);
            }
            rePath(l2);rePath(r2-1);
        }
        typename Segtree_Data::node_t query(int l, int r){
            push(l+N); push(r+N-1);
            typename Segtree_Data::node_t retL=Segtree_Data::node_ne(), retR=Segtree_Data::node_ne();
            for(l+=N, r+=N;l<r;l>>=1, r>>=1){
                if(l&1)retL=Segtree_Data::merge_nodes(retL, data[l++]);
                if(r&1)retR=Segtree_Data::merge_nodes(data[--r], retR);
            }
            return Segtree_Data::merge_nodes(retL, retR);
        }
    };
    int n;
    vector<vector<int> > const&g;
    vector<int> p, hp, d, heavy, pos;
    Segtree st;
    int dfs(int c){
        int s = 1, mas = 0;
        for(int e:g[c])
            if(e!=p[c]){
                d[e] = d[c]+1; p[e] = c;
                int es = dfs(e); s +=es;
                if(es>mas){
                    mas = es; heavy[c] = e;
                }
            }
        return s;
    }
    void init(){
        dfs(0);
        for(int i=0, lpos=0;i<n;++i)
            if(p[i]==-1 || i!=heavy[p[i]])
                for(int j=i;~j;j=heavy[j]){
                    hp[j] = i;
                    pos[j] = lpos++;
                }
    }
    template<class OP>
    void path(int u, int v, OP const&op){
        for(;hp[u]!=hp[v];v=p[hp[v]]){
            if(d[hp[u]]>d[hp[v]]) swap(u, v);
            op(pos[hp[v]], pos[v]+1);
        }
        if(d[u]>d[v]) swap(u, v);
        op(pos[u], pos[v]+1); // node aggregate
        //op(pos[u]+1, pos[v]+1); // edge aggregate
    }
public:
    Heavy_Light(vector<vector<int> > const&G):n(G.size()), g(G), p(n, -1), hp(n, -1), d(n, 0), heavy(n, -1), pos(n, -1), st(n){init();}
	Heavy_Light(vector<vector<int> > const&G, vector<typename Segtree_Data::node_t> const&_base):n(G.size()), g(G), p(n, -1), hp(n, -1), d(n, 0), heavy(n, -1), pos(n, -1), st(n){
        init();
        vector<typename Segtree_Data::node_t> base(n);
        for(int i=0;i<n;++i) base[pos[i]] = _base[i];
        st = Segtree(base);
    }

    void update(int u, int v, typename Segtree_Data::update_t const&val){
        path(u, v, [this, &val](int l, int r){st.update(l, r, val);});
    }
    typename Segtree_Data::node_t query(int u, int v){
        typename Segtree_Data::node_t ans = Segtree_Data::node_ne();
        path(u, v, [this, &ans](int l, int r){ans = Segtree_Data::merge_nodes(ans, st.query(l, r));});
        return ans;
    }
    int lca(int u, int v){
        for(;hp[u]!=hp[v];v=p[hp[v]]){
            if(d[hp[u]]>d[hp[v]]) swap(u, v);
        }
        if(d[u]>d[v]) swap(u, v);
        return u;
    }
};
struct Segtree_Data{
    typedef int node_t;
    typedef int update_t;
    static constexpr node_t node_ne() {return numeric_limits<int>::min();}
    static constexpr node_t node_init() {return numeric_limits<int>::min();}
    static constexpr update_t update_ne(){return numeric_limits<int>::max();}
    static node_t merge_nodes(node_t const&left, node_t const&right){
        return node_t(max(left,right));
    }
    static void update_node(node_t &node, update_t const&update){
        if(node == node_ne() || update<node) node = update;
    }
    static update_t merge_lazy(update_t const&old, update_t const&update){
        return min(old, update);
    }
};
