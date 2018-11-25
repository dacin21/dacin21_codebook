/*
 *	Fully dynamic connectiviy
 *	add or remove edges in O(log n^2)
 *	query connectivity in O(log n)
 *	code by dacin21, got cleaned up in september 2017
 *	0.33 seconds for n=m=1e5 on spoj
 */

#ifdef LOCAL_RUN
#define asser(x) do{if(1){assert(x);}}while(0)
#define asser2(x) do{if(1){assert(x);}}while(0)
#else
#define asser(x) do{if(0){assert(x);}}while(0)
#define asser2(x) do{if(0){assert(x);}}while(0)
#endif

struct Treap{
    struct Node{
		static mt19937 rng;
		static bool rng_init;
        Node*l, *r, *p;
        size_t y;
        unsigned int size;
        //int from, to;
        char mark, sub_mark; // used to find edges on current level
        Node():l(0), r(0), p(0), y(rng()), size(1), mark(0), sub_mark(0){if(!rng_init) rng = mt19937(std::chrono::duration_cast<std::chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count()); rng_init=1;}
        //Node(int _from, int _to):Node(){from=_from; to=_to;}
        Node* recalc(){
            sub_mark = mark;
            size = 1;
            if(l){
                sub_mark|=l->sub_mark;
                size+=l->size;
            }
            if(r){
                sub_mark|=r->sub_mark;
                size+=r->size;
            }
            return this;
        }
        Node* set_l(Node*_l){
            l = _l;
            if(l) l->p = this;
            return recalc();
        }
        Node* set_r(Node*_r){
            r = _r;
            if(r) r->p = this;
            return recalc();
        }
        Node* set_ch(Node*_l, Node*_r){
            l = _l; r = _r;
            if(l) l->p = this;
            if(r) r->p = this;
            return recalc();
        }
    };
    static Node* root(Node*a){
        if(!a) return a;
        while(a->p) a = a->p;
        return a;
    }
    static unsigned int size(Node*n){
        return n?n->size:0;
    }
    // splits tree into <c, =c, >c
    static pair<Node*, Node*> split(Node* c){
        Node*l = c->l, *r = c->r;
        if(l) l->p = 0;
        if(r) r->p = 0;
        c->set_ch(0, 0);
        while(c->p){
            Node* p = c->p;
            c->p = 0;
            if(p->l == c){
                r = p->set_l(r);
            } else {
                l = p->set_r(l);
            }
            c = p;
        }
        return make_pair(l, r);
    }
    // splits tree into <c, >=c
    static pair<Node*, Node*> lower_split(Node*c){
        Node*l = c->l, *r = c;
        if(l) l->p = 0;
        c->set_l(0);
        while(c->p){
            Node* p = c->p;
            c->p = 0;
            if(p->l == c){
                r = p->set_l(r);
            } else {
                l = p->set_r(l);
            }
            c = p;
        }
        return make_pair(l, r);
    }
    static Node* join(Node*l, Node*r){
        if(!l) return r;
        if(!r) return l;
        unsigned int depth;
        // bit-mask is used to reduce recursion depth, idea from anta (codeforces)
        unsigned long long path_mask = 0;
        Node*m = 0;
        for(depth=1;;++depth){
            // with high probability, this does not get called
            if(depth+1>= sizeof(path_mask) * CHAR_BIT){
                m = join(l->r, r->l);
                ++depth;
                path_mask = path_mask<<2|1<<(l->y>=r->y);
                break;
            }
            if(l->y < r->y){
                path_mask = path_mask<<1;
                Node*c = l->r;
                if(!c){
                    m = r;
                    r = r->p;
                    break;
                }
                l = c;
            } else {
                path_mask = path_mask<<1|1;
                Node*c = r->l;
                if(!c){
                    m = l;
                    l = l->p;
                    break;
                }
                r = c;
            }
        }
        // m: middle tree, l: next left ancestor, r: next right ancestor
        while(depth--){
            if(path_mask&1){
                m = r->set_l(m);
                r = r->p;
            } else {
                m = l->set_r(m);
                l = l->p;
            }
            path_mask>>=1;
        }
        return m;
    }
    // move a to the front
    static Node* evert(Node*a){
        pair<Node*, Node*> t = lower_split(a);
        return join(t.second, t.first);
    }
    static Node* push_front(Node*t, Node*l){
        t = root(t);
        asser(l->size == 1);
        if(!t) return l;
        for(;;){
            if(l->y < t->y){
                Node* tp = t->p;
                l->set_r(t);
                t = tp;
                break;
            }
            if(!t->l){
                break;
            }
            t = t->l;
        }
        while(t){
            l = t->set_l(l);
            t = t->p;
        }
        return l;
    }
    static Node* push_back(Node*t, Node*r){
        t = root(t);
        asser(r->size == 1);
        if(!t) return r;
        for(;;){
            if(r->y < t->y){
                Node* tp = t->p;
                r->set_l(t);
                t = tp;
                break;
            }
            if(!t->r){
                break;
            }
            t = t->r;
        }
        while(t){
            r = t->set_r(r);
            t = t->p;
        }
        return r;
    }
    static Node* begin(Node*a){
        a = root(a);
        while(a->l) a = a->l;
        return a;
    }
    static Node* end(Node*a){
        a = root(a);
        while(a->r) a = a->r;
        return a;
    }
    static unsigned int lower_cnt(Node*a){
        unsigned int ret = size(a->l);
        for(Node*p = a->p;p;p=a->p){
            if(a==p->r){
                ret+=size(p->l)+1;
            }
            a = p;
        }
        return ret;
    }
	static Node* get(Node*a, unsigned int ind){
	    assert(ind<size(a));
		for(;;){
			assert(a);
			if(ind == size(a->l)) return a;
			if(ind < size(a->l)) a = a->l;
			else {
				ind-=size(a->l)+1;
				a = a->r;
			}
		}
	}
    static void re_path(Node*a){
        for(;a;a=a->p){
            a->recalc();
        }
    }
};
mt19937 Treap::Node::rng = mt19937(83466);
bool Treap::Node::rng_init = false;

struct Euler_Tour_Tree{
    vector<tuple<int, int, bool> > const&graph;
    vector<Treap::Node> edges;
    vector<Treap::Node> first_edge;
    Euler_Tour_Tree(int n, int m, vector<tuple<int, int, bool> >&_graph):graph(_graph), edges(2*m), first_edge(n){
        for(int i=0;i<n;++i){
            first_edge[i] = Treap::Node();
        }
    }
    bool is_single(int a){
        return first_edge[a].size == 0 && !first_edge[a].p;
    }
    bool connected(int a, int b){
        return Treap::root(&first_edge[a]) == Treap::root(&first_edge[b]);
    }
    void reroot(int a){
        Treap::evert(&first_edge[a]);
    }
    void link(int edge_index, char mark){
        int a = get<0>(graph[edge_index]), b = get<1>(graph[edge_index]);
        Treap::Node* e_ab = new (edges.data()+2*edge_index) Treap::Node();
        Treap::Node* e_ba = new (edges.data()+2*edge_index+1) Treap::Node();
        e_ab->mark = mark;
        e_ab->recalc();
        Treap::evert(&first_edge[a]);
        Treap::push_back(&first_edge[a], e_ab);
        Treap::evert(&first_edge[b]);
        Treap::push_back(&first_edge[b], e_ba);
        Treap::join(Treap::root(&first_edge[a]), Treap::root(&first_edge[b]));
    }
    void cut(int edge_index){
        Treap::Node* e_ab = &edges[2*edge_index], *e_ba = &edges[2*edge_index|1];
        pair<Treap::Node*, Treap::Node*> ta = Treap::split(e_ab);
        unsigned int tar_size = Treap::size(ta.second);
        pair<Treap::Node*, Treap::Node*> tb = Treap::split(e_ba);
        // ensure ta is to the left of tb
        if(ta.second != e_ba && tar_size==Treap::size(ta.second)){
            swap(e_ab, e_ba);
            swap(ta, tb);
        }
        // order is ta.first, e_ab, ta.second or tb.first, e_ba, tb.second
        Treap::join(Treap::root(tb.second), Treap::root(ta.first));
    }
    unsigned int size(int a){
        if(is_single(a)) return 1;
        return (Treap::size(Treap::root(&first_edge[a]))+2)/3;
    }
    void set_mark(int a, char mark){
        first_edge[a].mark = mark;
        Treap::re_path(&first_edge[a]);
    }
    void set_edge_mark(int edge_index, char mark){
        edges[2*edge_index].mark = mark;
        Treap::re_path(&edges[2*edge_index]);
    }
    void set_edge_mark(Treap::Node*c, char mark){
        asser(c);
        if((c-edges.data())%2) return set_edge_mark((c-edges.data())/2, mark);
        c->mark = mark;
        Treap::re_path(c);
    }
    // calls op on marked nodes until op returns true or there are no such nodes
    template<char mark_mask, class OP>
    bool call_on_nodes(int a, OP &&op){
        Treap::Node*c = &first_edge[a];
        while(c->p && !(c->sub_mark&mark_mask)) c = c->p;
        while(c->sub_mark&mark_mask){
            if(c->mark&mark_mask){
                if(op(c)) return true;
            }
            // escape subtree if mark got unset
            while(c->p && !(c->sub_mark&mark_mask)) c = c->p;
            // find successor with mark_mask set in sub_mark
            if(c->r && c->r->sub_mark&mark_mask){
                c = c->r;
                while(c->l && c->l->sub_mark&mark_mask) c = c->l;
            } else {
                while(c->p && c == c->p->r) c = c->p;
                if(!c->p){
                    while(c->l && c->l->sub_mark&mark_mask) c = c->l;
                } else {
                    c = c->p;
                }
            }
        }
        asser(!(Treap::root(&first_edge[a])->sub_mark & mark_mask));
        return false;
    }
    unsigned int index_edge(Treap::Node* n){
        return (n-edges.data())/2;
    }
    unsigned int index_vertex(Treap::Node*n){
        return n - first_edge.data();
    }
};

struct Layer_Structure{
    int logn, n, m;
    vector<Euler_Tour_Tree> levels;
    vector<tuple<int, int, bool> > edges;
    vector<unsigned int> free_edges;
    // stores all edge and their index
    vector<map<int, int> > graph;
    Layer_Structure(int _n, int _m):n(_n), m(_m){
        for(logn=1;1ll<<logn<=n;++logn);
        levels.reserve(logn);
        levels.emplace_back(n, m, edges);
        edges.reserve(m);
        graph.resize(n*logn);
    }
    int calc_graph_pos(int level, int vertex){
        return level*n+vertex;
    }
    bool connected(int a, int b){
        return levels[0].connected(a, b);
    }
    void insert_edge(int edge_index, int level){
        int a = get<0>(edges[edge_index]), b = get<1>(edges[edge_index]);
        bool add_to_forest = get<2>(edges[edge_index]);
        asser(a!=b);
        int x = calc_graph_pos(level, a);
        int y = calc_graph_pos(level, b);
        asser(graph[x].find(b) == graph[x].end());
        asser(graph[y].find(a) == graph[y].end());
        if(add_to_forest){
            levels[level].link(edge_index, 2);
        }
        if(graph[x].empty()){
            levels[level].set_mark(a, 1);
        }
        if(graph[y].empty()){
            levels[level].set_mark(b, 1);
        }
        graph[x][b] = edge_index;
        graph[y][a] = edge_index;
    }
    void remove_edge(unsigned int edge_index, int level){
        int a = get<0>(edges[edge_index]), b = get<1>(edges[edge_index]);
        int x = calc_graph_pos(level, a), y = calc_graph_pos(level, b);
        asser(graph[x].find(b)!=graph[x].end());
        asser(graph[y].find(a)!=graph[y].end());
        asser(graph[x][b] == graph[y][a]);
        graph[x].erase(b);
        if(graph[x].empty()) levels[level].set_mark(a, 0);
        graph[y].erase(a);
        if(graph[y].empty()) levels[level].set_mark(b, 0);
    }
    void link(int a,int b){
        bool is_tree = !levels[0].connected(a, b);
        unsigned int edge_index;
        if(!free_edges.empty()){
            edge_index = free_edges.back();
            free_edges.pop_back();
        } else {
            edge_index = edges.size();
            edges.emplace_back();
        }
        edges[edge_index] = make_tuple(a, b, is_tree);
        insert_edge(edge_index, 0);
    }
    void push_edge(unsigned int edge_index, int level){
        asser(level+1<logn);
        if(level+1== (int)levels.size()) levels.emplace_back(n, m, edges);

        remove_edge(edge_index, level);
        insert_edge(edge_index, level+1);
    }
    void cut(int a, int b){
        int level = (int)levels.size()-1;
        while(graph[calc_graph_pos(level, a)].find(b) == graph[calc_graph_pos(level, a)].end()){
            asser(level>0);
            --level;
        }
        unsigned int edge_index = graph[calc_graph_pos(level, a)][b];
        asser2(levels[level].connected(a, b));
        bool was_in_forest = get<2>(edges[edge_index]);
        remove_edge(edge_index, level);
        free_edges.push_back(edge_index);
        if(!was_in_forest) return;
        asser2(level+1 == (int)levels.size() || !levels[level+1].connected(a, b));
        for(int i=level;i>=0;--i){
            levels[i].cut(edge_index);
            int small = a, big = b;
            if(levels[i].size(small)>levels[i].size(b)){
                swap(small, big);
            }
            // push tree edges first
            levels[i].call_on_nodes<2>(small, [this, i](Treap::Node*c){
                auto const& edge = edges[levels[i].index_edge(c)];
                int x = calc_graph_pos(i, get<0>(edge));
                if(graph[x].find(get<1>(edge))!=graph[x].end()){
                    push_edge(levels[i].index_edge(c), i);
                }
                levels[i].set_edge_mark(c, 0);
                if(graph[x].empty()){
                    levels[i].set_mark(get<0>(edge), 0);
                }
                return false;
            });
            // search for non-tree edges
            if(levels[i].call_on_nodes<1>(small, [this, a, b, i, small, edge_index](Treap::Node*c){
                int from = levels[i].index_vertex(c);
                int x = calc_graph_pos(i, from);
                asser(!graph[x].empty());
                for(auto it = graph[x].begin();it!=graph[x].end();){
                    // edge in same tree
                    if(levels[i].connected(small, it->first)){
                        unsigned int pushed_edge_index = it->second;
                        ++it;
                        push_edge(pushed_edge_index, i);
                    // replacement edge
                    } else {
                        get<2>(edges[it->second]) = true;
                        levels[i].link(it->second, 2);
                        for(int j=i-1;j>=0;--j){
                            levels[j].cut(edge_index);
                            levels[j].link(it->second, 0);
                        }
                        return true;
                    }
                }
                asser(graph[x].empty());
                levels[i].set_mark(from, 0);
                return false;
            })){
                break;
            }
        }
    }
};

namespace FIO{
    char buf[32*1042|1];
    int bc=0, be=0;
    char gc(){
        if(bc>=be){
            be = fread(buf, 1, sizeof(buf)-1, stdin);
            buf[be] = bc = 0;
        }
        return buf[bc++];
    }
    void readint(){}
    void readuint(){}
    template<typename T, typename ...S>
    void readuint(T &a, S& ...b){
        a=0;
        int x=gc();
        while(x<'0' || x>'9') x=gc();
        while(x>='0' && x<='9'){
            a = a*10+x-'0'; x=gc();
        }
        readuint(b...);
    }
	template<typename T, typename ...S>
    void readint(T &a, S& ...b){
        a=0;
        int x=gc(), s=1;;
        while(x!='-' && (x<'0' || x>'9')) x=gc();
		if(x=='-'){ s=-s; x=gc(); }
        while(x>='0' && x<='9'){
            a = a*10+x-'0'; x=gc();
        }
		if(s<0) a=-a;
        readint(b...);
    }
}
using FIO::readuint;

int dynacon2(){
  cin.tie(0);
  ios_base::sync_with_stdio(false);
  int N, M;
  readuint(N, M);
  Layer_Structure l(N+1, M+1);
  char c;
  int a, b;
  for(;M>0;--M){
    c = FIO::gc();
    while(c<'a' || c > 'z') c = FIO::gc();
    readuint(a, b);
    if(c=='a') l.link(a, b);
    else if(c == 'c') cout << (l.connected(a, b) ? "YES\n" : "NO\n");
    else l.cut(a, b);
  }
  return 0;
}
signed main(){
    return dynacon2();
}
