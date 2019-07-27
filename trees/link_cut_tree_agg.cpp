// Link-cut tree with path aggregates
#define lct_assert(x) assert(x)
struct LinkCutTree{
    struct Node{
        Node*p=0, *pp=0, *c[2]={0, 0};
        int id;
        bool ev=0;
        pair<int, int> agg, sub_agg;
        Node(int _id=-1):id(_id), agg(-1e9, -1), sub_agg(-1e9, -1){}
        void recalc(){
            sub_agg = agg;
            for(int i:{0, 1})
                if(c[i]){
                    c[i]->p=this;
                    sub_agg = max(sub_agg, c[i]->sub_agg);
                }
        }
        void push(){
            if(ev){
                ev=0;
                swap(c[0], c[1]);
                for(int i:{0, 1})
                    if(c[i]) c[i]->ev^=1;
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
            if(p) p->c[j]= this;
            c[!i]->recalc(); recalc();
            if(p) p->recalc();
        }
        Node* splay(){
            for(push();p;){
                if(p->p) p->p->push();
                p->push(); push();
                if(p->up()==up())
                    p->rot();
                rot();
            }
            return this;
        }
        Node*first(){
            push();
            return c[0]?c[0]->first():splay();
        }
    };
    vector<Node> g;
    LinkCutTree(size_t N){
        g.reserve(N);
        for(size_t i=0;i<N;++i)
            g.emplace_back(i);
    }
    // Note: u and v will NOT be at the root after this.
    bool connected(int u, int v){
        access(&g[u]);
        Node*x = g[u].first();
        access(&g[v]);
        Node*y = g[v].first();
        return x==y;
    }
    void link(int u, int p){
        lct_assert(!connected(u, p));
        make_root(&g[u]);
        access(&g[p]);
        g[u].pp=&g[p];
    }
    void cut_up(int u){
        Node*x = &g[u];
        access(x);
        x->unlink(0);
    }
    void cut(int u, int v){
        lct_assert(connected(u, v));
        Node*x = &g[u],*y=&g[v];
        make_root(x); make_root(y);
        x->splay();
        lct_assert(x->pp?y==x->pp:(y==x->c[0] || (x->c[0] && y==x->c[0]->c[0])));
        if(x->pp) x->pp=0;
        else x->unlink(0);
    }
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
    int lca(int u, int v){
        access(&g[u]);
        return access(&g[v])->id;
    }
    pair<int, int> pathagg(int u, int v){
        make_root(&g[u]);
        access(&g[v]);
        return g[v].sub_agg;
    }
	void set_agg(int u, pair<int, int> const& new_agg){
		access(&g[u]);
		g[u].agg = new_agg;
		g[u].recalc();
	}
};