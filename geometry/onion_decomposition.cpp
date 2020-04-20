struct Point{
    Point& operator+=(Point const&o) {
        x+=o.x;
        y+=o.y;
        return *this;
    }
    Point& operator-=(Point const&o) {
        x-=o.x;
        y-=o.y;
        return *this;
    }
    Point operator-() const {
        return Point{-x, -y};
    }
    friend Point operator+(Point a, Point const&b){
        a+=b;
        return a;
    }
    friend Point operator-(Point a, Point const&b){
        a-=b;
        return a;
    }
    friend ll dot(Point const&a, Point const&b){
        return a.x*b.y + a.y*b.y;
    }
    friend ll cross(Point const&a, Point const&b){
        return a.x*b.y - a.y*b.x;
    }
    friend bool operator<(Point const&a, Point const&b){
        return tie(a.x, a.y) < tie(b.x, b.y);
    }

    ll x=0, y=0;
};

int ccw(Point const&a, Point const&b){
    ll x = cross(a, b);
    return (x>0)-(x<0);
}
int ccw(Point const&a, Point const&b, Point const&c){
    return ccw(b-a, c-a);
}
// Decremental convex hull in O(n log n)
// From "Applications of a semi-dynamic convex hull algorithm" by J. Hershberger and S. Suri
struct Upper_Hull{
    struct Link{
        Point p;
        Link *prev = nullptr, *next = nullptr;
        int id;
    };
    struct Node{
        Link *chain, *chain_back;
        Link *tangent;
    };
    void fix_chain(int u, Link*l, Link*r){
        for(;l->next || r->next;){
            if(!r->next || (l->next && ccw(l->next->p - l->p, r->next->p - r->p) <= 0)){
                if(ccw(l->p, l->next->p, r->p) <= 0){
                    l = l->next;
                } else {
                    break;
                }
            } else {
                if(ccw(l->p, r->p, r->next->p) > 0){
                    r = r->next;
                } else {
                    break;
                }
            }
        }
        tree[u].tangent = l;
        tree[u].chain = tree[2*u].chain;
        tree[u].chain_back = tree[2*u+1].chain_back;
        tree[2*u].chain = l->next;
        tree[2*u+1].chain_back = r->prev;
        if(l->next){
            l->next->prev = nullptr;
        } else {
            tree[2*u].chain_back = nullptr;
        }
        if(r->prev){
            r->prev->next = nullptr;
        } else {
            tree[2*u+1].chain = nullptr;
        }
        l->next = r;
        r->prev = l;
    }
    void build(int u, int a, int b){
        if(b-a == 1){
            tree[u].chain = tree[u].chain_back = &lists[a];
            tree[u].tangent = nullptr;
            return;
        }
        const int m = a + (b-a)/2;
        build(2*u, a, m);
        build(2*u+1, m, b);
        auto l = tree[2*u].chain, r = tree[2*u+1].chain;
        fix_chain(u, l, r);
    }

    void rob(int u, int v){
            tree[u].chain = tree[v].chain;
            tree[v].chain = nullptr;
            tree[u].chain_back = tree[v].chain_back;
            tree[v].chain_back = nullptr;
    }

    void remove(int u, int a, int b, int const&i){
        if(i < a || i >= b) return;
        // we should never hit a leaf
        assert(b-a > 1);
        const int m = a + (b-a)/2;
        // one child -> that child contains i
        if(!tree[u].tangent){
            int v = i<m ? 2*u : 2*u+1;
            tree[v].chain = tree[u].chain;
            tree[v].chain_back = tree[u].chain_back;
            if(i < m){
                remove(2*u, a, m, i);
            } else {
                remove(2*u+1, m, b, i);
            }
            rob(u, v);
            return;
        }
        // restore hull of children
        auto l = tree[u].tangent, r = l->next;
        l->next = tree[2*u].chain;
        if(tree[2*u].chain){
            tree[2*u].chain->prev = l;
        } else {
            tree[2*u].chain_back = l;
        }
        tree[2*u].chain = tree[u].chain;
        r->prev = tree[2*u+1].chain_back;
        if(tree[2*u+1].chain_back){
            tree[2*u+1].chain_back->next = r;
        } else {
            tree[2*u+1].chain = r;
        }
        tree[2*u+1].chain_back = tree[u].chain_back;
        // delete i
        const int v = i<m ? 2*u : 2*u+1;
        // only i
        if(tree[v].chain == tree[v].chain_back && tree[v].chain->id == i){
            tree[v].chain = tree[v].chain_back = nullptr;
            rob(u, v^1);
            tree[u].tangent = nullptr;
            return;
        }
        if(i < m){
            if(l->id == i){
                l = l->prev;
            }
            remove(2*u, a, m, i);
            if(!l){
                l = tree[2*u].chain;
            }
            // more points on right half might get visible
            while(r->prev && ccw(l->p, r->prev->p, r->p) <= 0){
                r = r->prev;
            }
        } else {
            if(r->id == i){
                r = r->prev;
            }
            remove(2*u+1, m, b, i);
            if(!r){
                r = tree[2*u+1].chain;
            }
        }
        fix_chain(u, l, r);
    }
    void remove(int i){
        // i is the only point
        if(tree[1].chain == tree[1].chain_back){
            tree[1].chain = tree[1].chain_back = nullptr;
            return;
        }
        remove(1, 0, n, i);
    }
    Upper_Hull(vector<Point> const&v) : n(v.size()), tree(4*n), lists(n){
        assert(is_sorted(v.begin(), v.end()));
        for(int i=0; i<n; ++i){
            lists[i].p = v[i];
            lists[i].id = i;
        }
        build(1, 0, n);
    }
    vector<int> get_hull(){
        vector<int> ret;
        for(Link* u = tree[1].chain; u; u=u->next){
            ret.push_back(u->id);
        }
        return ret;
    }
    vector<Point> get_hull_points(){
        vector<Point> ret;
        for(Link* u = tree[1].chain; u; u=u->next){
            ret.push_back(u->p);
        }
        return ret;
    }

    int n;
    vector<Node> tree;
    vector<Link> lists;
};


// test code from https://codeforces.com/blog/entry/75929
signed main(){
    int N;
    scanf("%d",&N);
    vector<int> layer(N);
    vector<int> ans(N);
    vector<Point> ps;
    map<Point,int> id;
    for(int i=0;i<N;i++){
        int X,Y;
        scanf("%d %d",&X,&Y);
        ps.push_back({X,Y});
        id[{X,Y}]=i;
    }

    sort(ps.begin(),ps.end());
    Upper_Hull left(ps);
    reverse(ps.begin(),ps.end());
    for(auto& p:ps){
        p=-p;
    }
    Upper_Hull right(ps);
    for(auto& p:ps){
        p=-p;
    }
    reverse(ps.begin(),ps.end());
    for(int l=1,cnt=0;cnt<N;l++){
        set<int> hull;
        for(int i:left.get_hull()){
            hull.insert(i);
        }
        for(int i:right.get_hull()){
            hull.insert(N-1-i);
        }
        for(int i:hull){
            assert(!layer[i]);
            cnt++;
            layer[i]=l;
            left.remove(i);
            right.remove(N-1-i);
        }
    }
    for(int i=0;i<N;i++){
        ans[id[ps[i]]]=layer[i];
    }
    for(int i=0;i<N;i++){
        printf("%d\n",ans[i]);
    }
    return 0;
}
