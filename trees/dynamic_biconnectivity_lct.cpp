// off-line dynamic bridges
// runs in O(n + q log(n))
// uses link-cut trees and is faster than divide&conquer

template<typename S>
inline void xmin(S&a, S const&b){
    a=min(a, b);
}
template<typename S>
inline void xmax(S&a, S const&b){
    a=max(a, b);
}
template<typename S>
inline void xmin_2(S&a, S const&b, S const& inf){
    if((a>b && b!=inf) || a==inf) a=b;
}
template<typename S>
inline void xmax_2(S&a, S const&b, S const& neg_inf){
    if((a<b && b!=neg_inf) || a==neg_inf) a=b;
}

// Link-cut tree with path aggregates
#define lct_assert(x) assert(x)
//#define lct_assert(x) do{}while(0);
class Link_Cut_Tree{
public:
    struct Agg{
        static constexpr int no_cover = numeric_limits<int>::min();
        static constexpr int no_deletion = numeric_limits<int>::max();
        int edge_cnt, covered_cnt_subtree;
        int cover_lazy, cover_me, cover_subtree;
        int deletion_me, deletion_min;
        Agg(bool is_edge):edge_cnt(is_edge), covered_cnt_subtree(0), cover_lazy(no_cover), cover_me(no_cover), cover_subtree(no_cover), deletion_me(no_deletion), deletion_min(no_deletion){}
        int get_covered_cnt(){
            return cover_lazy == no_cover ? covered_cnt_subtree : edge_cnt;
        }
    };
    struct Node{
        static constexpr int neg_inf = -1e9;
        Node*p=0, *pp=0, *c[2]={0, 0};
        Agg agg;
        int id;
        bool ev=0, is_edge;
        Node(int _id=-1, bool _is_edge=false):agg(_is_edge), id(_id), is_edge(_is_edge){}
        void push(){
            if(ev){
                ev=0;
                for(int i:{0, 1})
                    if(c[i]) c[i]->ev^=1;
                swap(c[0], c[1]);
            }
            if(agg.cover_lazy != Agg::no_cover){
                xmax(agg.cover_me, agg.cover_lazy);
                xmin_2(agg.cover_subtree, agg.cover_lazy, Agg::no_cover);
                agg.covered_cnt_subtree = agg.edge_cnt;
                for(int i:{0, 1})
                    if(c[i]) xmax(c[i]->agg.cover_lazy, agg.cover_lazy);
                agg.cover_lazy = Agg::no_cover;
            }
        }
        void recalc(){
            push();
            assert(agg.cover_lazy == Agg::no_cover);
            agg.edge_cnt = is_edge;
            agg.covered_cnt_subtree = is_edge && (agg.cover_me != Agg::no_cover);
            agg.cover_subtree = agg.cover_me;
            agg.deletion_min = agg.deletion_me;
            for(int i:{0, 1})
                if(c[i]){
                    c[i]->p=this;
                    xmin_2(agg.cover_subtree, c[i]->agg.cover_subtree, Agg::no_cover);
                    xmin_2(agg.cover_subtree, c[i]->agg.cover_lazy, Agg::no_cover);
                    agg.edge_cnt+=c[i]->agg.edge_cnt;
                    agg.covered_cnt_subtree+=c[i]->agg.get_covered_cnt();
                    xmin(agg.deletion_min, c[i]->agg.deletion_min);
                }
        }
        void unlink(int i){
            lct_assert(c[i]);
            c[i]->p=0;
            c[i]=0;
            recalc();
        }
        int up(){return p?p->c[1]==this:-1;}
        void rot(){
            swap(pp, p->pp);
            int i = up(), j=p->up();
            p->c[i]=c[!i];
            c[!i]=p; p=p->p;
            c[!i]->recalc(); //recalc();
            if(p) {
                p->c[j]= this;
                //p->recalc();
            }
        }
        Node* splay(){
            for(push();p;){
                if(p->p) p->p->push();
                p->push(); push();
                if(p->up()==up()){
                    p->rot();
                    p->push();
                    push();
                }
                rot();
            }
            recalc();
            return this;
        }
        Node*first(){
            push();
            return c[0]?c[0]->first():splay();
        }
        // self-driving removal
        void remove_cover(int cover){
            if(agg.cover_lazy <= cover){
                agg.cover_lazy = Agg::no_cover;
            }
            if(agg.cover_subtree == Agg::no_cover || agg.cover_subtree > cover) return;
            if(agg.cover_me <= cover){
                agg.cover_me = Agg::no_cover;
            }
            for(int i:{0, 1})
                if(c[i]) c[i]->remove_cover(cover);
            recalc();
        }
    };
    vector<Node> g;
    // edges get counted twice, covered edge once
    int edges, covered_edges;
    Link_Cut_Tree(size_t n, size_t m){
        size_t N = n+m;
        g.reserve(N);
        for(size_t i=0;i<N;++i)
            g.emplace_back(i, i>=n);
        edges = covered_edges = 0;
    }
    bool connected(int u, int v){
        Node*x = access(&g[u])->first();
        Node*y = access(&g[v])->first();
        return x==y;
    }
    void link(int u, int p){
        lct_assert(!connected(u, p));
        make_root(&g[u]);
        g[u].pp=&g[p];
        ++edges;
    }
    void cut_up(int u){
        Node*x = &g[u];
        access(x);
        x->unlink(0);
        --edges;
    }
    void cut(int u, int v){
        lct_assert(connected(u, v));
        Node*x = &g[u],*y=&g[v];
        make_root(x); make_root(y);
        x->splay();
        lct_assert(x->pp?y==x->pp:(y==x->c[0] || (x->c[0] && y==x->c[0]->c[0])));
        if(x->pp) x->pp=0;
        else x->unlink(0);
        --edges;
    }
private:
    void make_root(Node*u){
        access(u);
        u->ev^=1;
        access(u);
    }
    Node*access(Node*x){
        Node*u = x;
        u->splay();
        while(Node*pp=u->pp){
            pp->splay();
            assert(pp->agg.cover_lazy == Agg::no_cover);
            u->pp=0;
            if(pp->c[1]){
                swap(pp->c[1]->p,pp->c[1]->pp);
            }
            pp->c[1]=u;
            pp->recalc();
            u=pp;
        }
        x->splay();
        if(x->c[1]){
            x->c[1]->pp=x;
            x->unlink(1);
        }
        return u;
    }
public:
    int lca(int u, int v){
        access(&g[u]);
        return access(&g[v])->id;
    }
    int path_min_cover(int u, int v){
        make_root(&g[u]);
        access(&g[v]);
        return g[v].agg.cover_subtree;
    }
    void cover(int u, int v, int cover_id){
        make_root(&g[u]);
        Node*n_v = &g[v];
        access(n_v);
        covered_edges-=n_v->agg.get_covered_cnt();
        n_v->agg.cover_lazy = cover_id;
        covered_edges+=n_v->agg.get_covered_cnt();
    }
    void un_cover(int u, int v, int cover_id){
        make_root(&g[u]);
        Node*n_v = &g[v];
        access(n_v);
        covered_edges-=n_v->agg.get_covered_cnt();
        n_v->remove_cover(cover_id);
        covered_edges+=n_v->agg.get_covered_cnt();
    }
    void set_deletion(int u, int tim){
        g[u].agg.deletion_me = tim;
        g[u].recalc();
        access(&g[u]);
    }
    int path_min_deletion(int u, int v){
        make_root(&g[u]);
        access(&g[v]);
        return g[v].agg.deletion_min;
    }
};

class Dynamic_Biconnectivity{
    struct Query{
        int u, v; // edge endpoints
        int deletion_time; // time of other event on edge
    };
    int n, q;
    Link_Cut_Tree lct;
    vector<map<int, int> > present_edges;
    vector<Query> queries;
    vector<char> is_tree_edge;
public:
    Dynamic_Biconnectivity(int _n, int _q):n(_n), q(_q), lct(n, q), present_edges(n){
        queries.reserve(q);
    }
    void add_edge(int a, int b){
        assert(0<=a && a<n && 0<=b && b<n && a!=b);
        if(a>b) swap(a, b);
        int tim = queries.size();
        present_edges[a][b] = tim;
        queries.push_back(Query{a, b, q+1});
    }
    void remove_edge(int a, int b){
        assert(0<=a && a<n && 0<=b && b<n && a!=b);
        if(a>b) swap(a, b);
        int tim = queries.size();
        auto it = present_edges[a].find(b);
        assert(it != present_edges[a].end());
        int insertion_time = it->second;
        present_edges[a].erase(it);
        queries[insertion_time].deletion_time = tim;
        queries.push_back(Query{a, b, insertion_time});
    }
    vector<int> compute(){
        is_tree_edge.assign(q, 0);
        vector<int> ret; ret.reserve(q);
        for(int i=0;i<q;++i){
            auto const& query = queries[i];
            if(query.deletion_time<i){ // deletion
                if(is_tree_edge[query.deletion_time]){
                    remove_edge_from_tree(query.deletion_time);
                } else {
                    lct.un_cover(query.u, query.v, i);
                }
            } else { // insertion
                if(lct.connected(query.u, query.v)){
                    int min_cover = lct.path_min_deletion(query.u, query.v);
                    const Query* covering_query = &query;
                    if(min_cover < query.deletion_time){
                        int min_cover_insertion_id = queries[min_cover].deletion_time;
                        covering_query = &queries[min_cover_insertion_id];
                        remove_edge_from_tree(min_cover_insertion_id);
                        insert_edge_into_tree(i);
                    }
                    lct.cover(covering_query->u, covering_query->v, covering_query->deletion_time);
                } else {
                    insert_edge_into_tree(i);
                }
            }
            ret.push_back(lct.edges/2 - lct.covered_edges);
        }
        return ret;
    }
private:
    void insert_edge_into_tree(int id){
        lct.set_deletion(n+id, queries[id].deletion_time);
        lct.link(queries[id].u, n+id);
        lct.link(queries[id].v, n+id);
        is_tree_edge[id] = true;
    }
    void remove_edge_from_tree(int id){
        lct.cut(queries[id].u, n+id);
        lct.cut(queries[id].v, n+id);
        lct.covered_edges-=lct.g[n+id].agg.get_covered_cnt();
        is_tree_edge[id] = false;
    }
};
constexpr int Link_Cut_Tree::Agg::no_cover;
constexpr int Link_Cut_Tree::Agg::no_deletion;
