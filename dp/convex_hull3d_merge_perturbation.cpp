// convex hull trick in 3D
// requires distinct monotone x-coordinates
// O(log(n)^2) per query or update (now for real)
// uses symbolic perturbation to deal with degenerate faces

ll searches = 0;
ll builds = 0;

mt19937 RNG(132132);
struct Persnode{
    array<Persnode*, 2> ch;
    array<ll, 3> d, d_next;
    ll y;
    Persnode():ch{0, 0}, d{0,0,0}, d_next{0,0,0}, y(RNG()){}
    Persnode(array<ll, 3> const&v):ch{0, 0}, d(v), d_next{0,0,0}, y(RNG()){}
    Persnode(Persnode* const&n):ch(n->ch), d(n->d), d_next(n->d_next), y(n->y){}
};
constexpr size_t pool_size = 10*1000000;
size_t pool_back = 0;
Persnode pool[pool_size];
Persnode* get_node(){
    #ifdef LOCAL_RUN
    assert(pool_back < pool_size);
    #endif // LOCAL_RUN
    return pool + (pool_back++);
}
void reset_pool(){
    pool_back = 0;
}

struct Perstree{
    static  pair<Persnode*, Persnode* > split(Persnode* const&t, ll const&x, array<ll, 3> *&link_next, array<ll, 3> *&next_d){
        if(!t) return pair<Persnode*, Persnode* > (0, 0);
        Persnode* ret = new(get_node()) Persnode(t);
        if(x<=t->d[0]){
            auto tmp = split(t->ch[0], x, link_next, next_d);
            ret->ch[0] = tmp.second;
            if(tmp.second == nullptr){
                next_d = &(ret->d);
            }
            return {tmp.first, ret};
        } else {
            auto tmp = split(t->ch[1], x, link_next, next_d);
            ret->ch[1] = tmp.first;
            if(tmp.first == nullptr){
                link_next = &(ret->d_next);
            }
            return{ret, tmp.second};
        }
    }
    static Persnode* merge(Persnode* const&t, Persnode* const&o){
        if(!t) return o;
        if(!o) return t;
        if(t->y>=o->y){
            Persnode* ret = new(get_node()) Persnode(t);
            ret->ch[1] = merge(ret->ch[1], o);
            return ret;
        } else {
            Persnode* ret = new(get_node()) Persnode(o);
            ret->ch[0] = merge(t, ret->ch[0]);
            return ret;
        }
    }
    Persnode* root;
    Perstree():root(0){}
    Perstree(Persnode* r):root(r){}
    Perstree insert_or_delete(array<ll, 3> const&p){
        array<ll, 3> * link = nullptr, *nextp = nullptr, *trash = nullptr;
        auto tmp = split(root, p[0], link, trash);
        auto tmp2 = split(tmp.second, p[0]+1, trash, nextp);
        if(!tmp2.first){
            if(link) *link= p;
            auto mid = new(get_node()) Persnode(p);
            if(nextp) mid->d_next = *nextp;
            return Perstree(merge(merge(tmp.first, mid),tmp2.second));
        }
        if(link && nextp){
            *link = *nextp;
        }
        return Perstree(merge(tmp.first, tmp2.second));
    }
    array<ll, 3> mascal(array<ll, 3> const&dir)const{
        assert(root);
        array<ll, 3> res = root->d;
        Persnode* cur = root;
        while(cur){
            if((cur->d[0]-res[0])*dir[0] + (cur->d[1]-res[1])*dir[1] + (cur->d[2]-res[2])*dir[2]>=0){
                res = cur->d;
            }
            ++searches;
            Persnode* r = cur->ch[1];
            if(!r) {
                cur = cur->ch[0];
            } else {
                if((cur->d[0]-cur->d_next[0])*dir[0] + (cur->d[1]-cur->d_next[1])*dir[1] + (cur->d[2]-cur->d_next[2])*dir[2]>=0){
                    cur = cur->ch[0];
                } else {
                    cur = cur->ch[1];
                }
            }
        }
        return res;
    }
    void print()const{
        stack<pair<Persnode*, int> > s;
        s.emplace(root, 0);
        while(!s.empty()){
            auto e = s.top();
            s.pop();
            if(e.first==0)continue;
            if(e.second){
                cerr << "  " << e.first->d[0] << "/" << e.first->d[1] << "/" << e.first->d[2] << " ";
                cerr << "|" << e.first->d_next[0] << "/" << e.first->d_next[1] << "/" << e.first->d_next[2] << " ";
            } else {
                s.emplace(e.first->ch[1], 0);
                s.emplace(e.first, 1);
                s.emplace(e.first->ch[0], 0);
            }
        }
        cerr << "\n";
    }
};


// might get away with int64_t if X*Y + X*Z fits.
using middle_t = __int128;

ostream& operator<<(ostream&o, __int128 const&x){
    if(x < 0) return o << '-' << -x;
    if(x < 1e18){
        return o << (int64_t) x;
    } else {
        return o << x/10 << (int)(x%10);
    }
}

struct Chull3D{
    using ll = long long;

    struct Frac{
        middle_t a, b;
        Frac(middle_t const&x=0, middle_t const&y=1):a(x), b(y){if(b<0){a=-a; b=-b;}}
        Frac& reduce(){
            middle_t g = __gcd(a, b);
            a/=g, b/=g;
            if(a<0){ a=-a, b=-b; }
            return *this;
        }
        Frac& operator+(Frac const&o)const{
            middle_t g = __gcd(b, o.b);
            return Frac(o.b/g*a + b/g*o.a, b/g*o.b).reduce();
        }
        Frac& operator-(Frac const&o)const{
            middle_t g = __gcd(b, o.b);
            return Frac(o.b/g*a - b/g*o.a, b/g*o.b).reduce();
        }
        Frac& operator*(Frac const&o)const{
            return Frac(a*o.a, b*o.b).reduce();
        }
        Frac& operator/(Frac const&o)const{
            return Frac(a*o.b, b*o.a).reduce();
        }
        Frac operator-()const{
            return Frac(-a, b);
        }
        bool operator==(Frac const&o)const{
            return ((a<0) == (o.a<0)) && (a*o.b == b*o.a);
        }
        bool operator<(Frac const&o)const{
            //return a/(long double)b < o.a/(long double)o.b;
            return a*(__int128)o.b < o.a*(__int128)b;
            //return a*(long double)o.b < o.a*(long double)b;
        }
        friend ostream& operator<<(ostream&o, Frac const&f){
            return o << (f.a/(long double)f.b);
        }
        explicit operator double() const {
            return (double)a / (double)b;
        }
    };
    struct Point{
        ll x, y, z;
        Point *prev, *next;
        Point():prev(0), next(0){}
        Point(ll _x, ll _y, ll _z, Point*p, Point*n):x(_x), y(_y), z(_z), prev(p), next(n){}
        void act() {
            if (prev->next != this) {
                //cerr << "in " << prev << "/" << this << "/" << next << "\n";
                prev->next = next->prev = this; // insert
            }
            else {
                //cerr << "out " << prev << "/" << this << "/" << next << "\n";
                prev->next = next; next->prev = prev;
            } // delete
        }
        tuple<ll, ll, ll> tup(){
            return {x, y, z};
        }
        array<ll,3> arr(){
            return {x, y, z};
        }
        friend ostream& operator<<(ostream&o, Point const&p){
            return o << p.x << "/" << p.y << "/" << p.z;
        }
    };
    Point cross(Point const&a, Point const&b, Point const&c){
        Point x(b.x-a.x, b.y-a.y, b.z-a.z, 0, 0);
        Point y(c.x-a.x, c.y-a.y, c.z-a.z, 0, 0);
        return Point(x.y*y.z-x.z*y.y, x.z*y.x-x.x*y.z, x.x*y.y-x.y*y.x, 0, 0);
    }
    static ll inf;
    static Point nil;
    Point *NIL = &nil;
    struct Event{
        static int sign(array<ll, 3> const&v){
            for(int i=0; i<3; ++i){
                if(v[i]) return (v[i]>0) ? 1 : -1;
            }
            return 0;
        }
        // s=1 --> INF,   s=-1 --> -INF
        static Event get_inf(int s){
            Event ret;
            ret.num = s;
            ret.den = 0;
            ret.sgn = 1;
            return ret;
        }
        Event(){}
        Event(Point *a, Point *b, Point *c) {
            auto *p = a, *q = b, *r = c;
        /**
            degenerate version:
                (q->x-p->x)*(r->y-p->y) - (r->x-p->x)*(q->y-p->y)
            perturbation order:
                1 * (q->x-p->x)*(r->y-p->y) - (r->x-p->x)*(q->y-p->y)
                d(p->y) * (-(q->x-p->x) + (r->x-p->x))
                d(q->y) * -(r->x-p->x)
                d(r->y) * (q->x-p->x)
            numerator is similar with z instead of y
         */
            num = (q->x-p->x)*(middle_t)(r->z-p->z) - (r->x-p->x)*(middle_t)(q->z-p->z);
            numerator = {{
                (-(q->x-p->x) + (r->x-p->x)),
                -(r->x-p->x),
                (q->x-p->x)
            }};
            den = (q->x-p->x)*(middle_t)(r->y-p->y) - (r->x-p->x)*(middle_t)(q->y-p->y);
            denominator = {{
                (-(q->x-p->x) + (r->x-p->x)),
                -(r->x-p->x),
                (q->x-p->x)
            }};
            numerator_weights = {p->x, q->x, r->x};
            denominator_weights = {p->x, q->x, r->x};
            if(den){
                sgn = (den>0) - (den<0);
            } else {
                sgn = sign(denominator);
            }
        }
        // a/x <=> b/y
        static int frac_comp(middle_t a, middle_t b, middle_t x, middle_t y) {
            auto tmp = a*(__int128)b - x*(__int128)y;
            return (tmp>0) - (tmp<0);
        }
        static int frac_comp_long(middle_t a, middle_t b, middle_t x, middle_t y) {
            // split to go beyond 128 bits
            constexpr ll BLOCK = (1ll<<32);
            auto tmp = (a>>32)*(__int128)b - (x>>32)*(__int128)y;
            if(tmp > (__int128)1<<94) return 1;
            if(tmp < -((__int128)1<<94)) return -1;
            tmp = tmp*BLOCK + (a&(BLOCK-1))*(__int128)b - (x&(BLOCK-1))*(__int128)y;
            return (tmp>0) - (tmp<0);
        }
        // only non-degenerate term
        int single_comp(Event const&o, middle_t weight, middle_t weight_o) const {
            return frac_comp_long(num, weight, o.num, weight_o);
        }
        // only perturbed numerator
        int numer_comp(Event const&o, middle_t weight, middle_t weight_o) const {
            int ret;
            if((ret = single_comp(o, weight, weight_o))) return ret;
            auto const &N = numerator, &N_o = o.numerator;
            auto const &W = numerator_weights, &W_o = o.numerator_weights;
            int i = 0, j = 0;
            while(i<3 & j<3){
                if(W[i] > W_o[j]){
                    if((ret = frac_comp(0, 0, N_o[j], weight_o))) return ret;
                    ++j;
                } else if(W[i] < W_o[j]){
                    if((ret = frac_comp(N[i], weight, 0, 0))) return ret;
                    ++i;
                } else {
                    if((ret = frac_comp(N[i], weight, N_o[j], weight_o))) return ret;
                    ++i;
                    ++j;
                }
            }
            while(i<3){
                if((ret = frac_comp(N[i], weight, 0, 0))) return ret;
                ++i;
            }
            while(j<3){
                if((ret = frac_comp(0, 0, N_o[j], weight_o))) return ret;
                ++j;
            }
            /*while(i<3 || j<3){
                if(i == 3 || (j!=3 && W[i] > W_o[j])){
                    if((ret = frac_comp(0, 0, N_o[j], weight_o))) return ret;
                    ++j;
                } else if(j == 3 || W[i] < W_o[j]){
                    if((ret = frac_comp(N[i], weight, 0, 0))) return ret;
                    ++i;
                } else {
                    if((ret = frac_comp(N[i], weight, N_o[j], weight_o))) return ret;
                    ++i;
                    ++j;
                }
            }*/
            return 0;
        }
        int full_comp(Event const&o) const {
            int ret;
            if((ret = numer_comp(o, o.den, den))) return ret;
            auto const &D = denominator, &D_o = o. denominator;
            auto const &W = denominator_weights, &W_o = o.denominator_weights;
            int i = 0, j = 0;
            while(i<3 || j<3){
                if(i == 3 || (j!=3 && W[i] > W_o[j])){
                    if((ret = numer_comp(o, D_o[j], 0))) return ret;
                    ++j;
                } else if(j == 3 || W[i] < W_o[j]){
                    if((ret = numer_comp(o, 0, D[i]))) return ret;
                    ++i;
                } else {
                    if((ret = numer_comp(o, D_o[j], D[i]))) return ret;
                    ++i;
                    ++j;
                }
            }
            return 0;
        }
        int comp(Event const&o) const {
            int ret = full_comp(o);
            return sgn * o.sgn * ret;
        }
        int numerator_sign() const {
            if(num) return (num>0) - (num<0);
            return sign(numerator);
        }
        #define COMP_BASED_OP(x)\
        friend bool operator x (Event const&a, Event const&b){\
            return a.comp(b) x 0;\
        }
        COMP_BASED_OP(<)
        COMP_BASED_OP(<=)
        COMP_BASED_OP(>)
        COMP_BASED_OP(>=)
        COMP_BASED_OP(==)
        COMP_BASED_OP(!=)
        #undef COMP_BASED_OP

        explicit operator Frac() const {
            if(den){
                return Frac(sgn*num, sgn*den);
            }
            // now its one of: -inf, inf
            // distinct x-coordinates -> numerator is non-zero
            assert(numerator_sign());
            return Frac(numerator_sign()*sgn, 0);
        }

        //Point *p=nullptr, *q=nullptr, *r=nullptr;
        middle_t num{}, den{};
        array<ll, 3> numerator{}, denominator{};
        array<ll, 3> numerator_weights{}, denominator_weights{};
        int sgn;
    };
    Event INF = Event::get_inf(1);
    Event NEG_INF = Event::get_inf(-1);
    Event make_event(Point *a, Point *b, Point *c){
        if(a == NIL || b == NIL || c == NIL){
            return INF;
        }
        return Event(a, b, c);
    }
    Frac time(Point *p, Point *q, Point *r) { // when turn changes
        if (p == NIL || q == NIL || r == NIL) return (Frac)INF;
        return (Frac) Event(p, q, r);
    }
    Point *sort(Point P[], int n) { // mergesort
        Point *a, *b, *c, head;
        if (n == 1) { P[0].next = NIL; return P; }
        a = sort(P, n/2);
        b = sort(P+n/2, n-n/2);
        c = &head;
        do
            if (a->tup() < b->tup()) { c = c->next = a; a = a->next; }
            else { c = c->next = b; b = b->next; }
        while (c != NIL);
        return head.next;
    }
    void hull(Point *list, int n, Point **A, Point **B) { // the algorithm
        Point *u, *v, *mid;
        Event t[6], oldt, newt;
        int i, j, k, l, minl;
        if (n == 1) { A[0] = list->prev = list->next = NIL; return; }

        for (u = list, i = 0; i < n/2-1; u = u->next, i++) ;
        mid = v = u->next;
        hull(list, n/2, B, A); // recurse on left and right sides
        hull(mid, n-n/2, B+n/2*2, A+n/2*2);

        for (;;) // find initial bridge
            if (make_event(u, v, v->next).sgn <= 0) v = v->next;
            else if (make_event(u->prev, u, v).sgn <= 0) u = u->prev;
            else break;

        // merge by tracking bridge uv over time
        for (i = k = 0, j = n/2*2, oldt = NEG_INF; ; oldt = newt) {
            ++builds;
            // second doesn't matter apparently, just varies degenerate faces????
            t[0] = make_event(B[i]->prev, B[i], B[i]->next);
            t[1] = make_event(B[j]->prev, B[j], B[j]->next);
            t[2] = make_event(u->prev, u, v);
            t[3] = make_event(u, v, v->next);
            t[4] = make_event(u, u->next, v);
            t[5] = make_event(u, v->prev, v);
            for (newt = INF, l = 0; l < 6; l++)
            if (t[l] > oldt && t[l] < newt) { minl = l; newt = t[l]; }
            //for(int l=0;l<6;++l){
            //    cerr << setw(16) << (Frac)t[l];
            //}   cerr << " : " << setw(0) << minl << "\n";
            if (newt == INF) break;
            switch (minl) {
                case 0: if (B[i]->tup() < u->tup()) A[k++] = B[i]; B[i++]->act(); break;
                case 1: if (B[j]->tup() > v->tup()) A[k++] = B[j]; B[j++]->act(); break;
                case 2: A[k++] = u; u = u->prev; break;
                case 3: A[k++] = v; v = v->next; break;
                case 4: A[k++] = u = u->next; break;
                case 5: A[k++] = v = v->prev; break;
            }
        }
        A[k] = NIL;
        u->next = v; v->prev = u; // now go back in time to update pointers
        for (k--; k >= 0; k--) {
            if (A[k]->tup() <= u->tup() || A[k]->tup() >= v->tup()) {
                A[k]->act();
                if (A[k] == u) u = u->prev; else if (A[k] == v) v = v->next;
            } else {
                u->next = A[k];
                A[k]->prev = u;
                v->prev = A[k];
                A[k]->next = v;
                if (A[k]->tup() < mid->tup()) u = A[k];
                else v = A[k];
            }
        }
        //cerr << *u << "  " << *v << "\n";
        //print_times(A);
    }
    size_t n;
    vector<array<ll, 3> > points;
    vector<Point> P;
    vector<Point*> A, B;
    vector<pair<Frac, Perstree>> hulls;
    Point* list;
    int old_pool_back;
    Chull3D(vector<array<ll, 3> > const&p):n(p.size()), points(p){
        std::sort(points.begin(), points.end());
        for(auto const&e:points)
            P.emplace_back(e[0], e[1], e[2], (Point*)0, (Point*)0);
        list = sort(P.data(), n);
        A.resize(2*n, 0); B.resize(2*n, 0);
        hull(list, n, A.data(), B.data());
        build_hull();
    }
    void print_times(Point **A){
        Frac oldt = (Frac)NEG_INF;
        int i;
        stringstream ss;
        cerr << "start\n";
        for(i=0; A[i]!=NIL; A[i++]->act()){
            //auto now = time(A[i]->prev, A[i], A[i]->next);
            //if(now < oldt) cerr << oldt << "/" << now << "\n";
            oldt = time(A[i]->prev, A[i], A[i]->next);
            //cerr << A[i]->arr()[0] << "/" << A[i]->arr()[1] << "/" << A[i]->arr()[2] << "\n";
            ss << oldt << " ";
        }
        ss << "\n";
        cerr << "mid\n";
        for(;i>0;A[--i]->act());
        cerr << "done\n";
        cerr << ss.str();
    }
    void build_hull(){
        old_pool_back = pool_back;
        Perstree hul;
        for(Point*l = list;l!=NIL; l = l->next){
            hul = hul.insert_or_delete(l->arr());
        }
        Frac oldt = (Frac)NEG_INF;
        int i;
        for(i=0; A[i]!=NIL; A[i++]->act()){
            hulls.emplace_back(oldt, hul);
            hul = hul.insert_or_delete(A[i]->arr());
            auto now = time(A[i]->prev, A[i], A[i]->next);
            if(now < oldt){
                cerr << oldt << "/" << now << "\n";
                auto t = time(A[i]->prev, A[i], A[i]->next);
            }
            oldt = time(A[i]->prev, A[i], A[i]->next);
            //cerr << A[i]->arr()[0] << "/" << A[i]->arr()[1] << "/" << A[i]->arr()[2] << "\n";
            //cerr << oldt << " ";
        }
        //cerr << "\n";
        for(;i>0;A[--i]->act());
        hulls.emplace_back(oldt, hul);
        //for(auto &e:hulls) debug << (double)e.first << "  ";
        //debug << "\n";
        //for(auto const&e:hulls){
        //    cerr << e.first << ":  ";
        //    e.second.print();
        //} cerr << "\n";
    }
    void print_faces(){
        cout << "faces:\n";
        for (int i = 0; A[i] != NIL; A[i++]->act()){ // output
            cout << A[i]->prev-P.data() << " " << A[i]-P.data() << " " << A[i]->next-P.data() << " : ";
            cout << (A[i]->prev->next != A[i] ? cross(*A[i]->prev, *A[i], *A[i]->next):cross(*A[i]->next, *A[i], *A[i]->prev)) << "     ";
            //cout << cross(*A[i]->prev, *A[i], *A[i]->next);
            cout << "\n";
        }
    }
    array<ll, 3> extreme_point(array<ll, 3> const&dir)const{
        Frac t(-dir[1], dir[2]);
        auto it = lower_bound(hulls.begin(), hulls.end(), t,[](pair<Frac, Perstree>const&a, Frac const&b){return a.first<b;});
        if(it == hulls.begin()) return it->second.mascal(dir);
        if(it == hulls.end()) return prev(it)->second.mascal(dir);
        array<ll, 3> r1 = it->second.mascal(dir);
        array<ll, 3> r2 = prev(it)->second.mascal(dir);
        return ((r2[0]-r1[0])*dir[0] + (r2[1]-r1[1])*dir[1] + (r2[2]-r1[2])*dir[2])<0 ? r1:r2;
    }
    size_t size(){
        return n;
    }
    void merge(Chull3D &o){
        //cerr << "merge: " << size() << "/" << o.size() << "\n";
        n+=o.n;
        points.insert(points.end(), o.points.begin(), o.points.end());
        //debug << named(points) << "\n";
        vector<Point> pres(P.size()+o.P.size());
        copy(P.begin(), P.end(), pres.begin());
        copy(o.P.begin(), o.P.end(), pres.begin()+P.size());
        A.swap(B);
        A.resize(2*pres.size());
        o.A.swap(o.B);
        for(size_t i=0;i<P.size();++i){
            if(pres[i].next!=NIL) pres[i].next = pres[i].next-P.data()+pres.data();
            if(pres[i].prev!=NIL) pres[i].prev = pres[i].prev-P.data()+pres.data();
        }
        for(size_t i=P.size();i<P.size()+o.P.size();++i){
            if(pres[i].next!=NIL) pres[i].next = pres[i].next-o.P.data()+pres.data()+P.size();
            if(pres[i].prev!=NIL) pres[i].prev = pres[i].prev-o.P.data()+pres.data()+P.size();
        }
        for(int i=0;B[i]!=NIL;++i){
            B[i] = B[i]-P.data()+pres.data();
        }
        for(int i=0;o.B[i]!=NIL;++i){
            o.B[i] = o.B[i]-o.P.data()+pres.data()+P.size();
        }
        list = list-P.data()+pres.data();
        o.list = o.list-o.P.data()+pres.data()+P.size();
        Point *u, *v, *mid;
        Event t[6], oldt, newt;
        int i, j, k, l, minl;

        for (u = list; u->next!=NIL; u = u->next);
        mid = v = o.list;

        for (;;) // find initial bridge
            if (make_event(u, v, v->next).sgn <= 0) v = v->next;
            else if (make_event(u->prev, u, v).sgn <= 0) u = u->prev;
            else break;

        // merge by tracking bridge uv over time
        for (i = k = 0, j = 0, oldt = NEG_INF; ; oldt = newt) {
            ++builds;
            /*t[0] = make_event(B[i]->prev, B[i], B[i]->next);
            t[1] = make_event(o.B[j]->prev, o.B[j], o.B[j]->next);
            t[2] = make_event(u->prev, u, v);
            t[3] = make_event(u, v, v->next);
            t[4] = make_event(u, u->next, v);
            t[5] = make_event(u, v->prev, v);
            for (newt = INF, l = 0; l < 6; l++)
            if (t[l] > oldt && t[l] < newt) { minl = l; newt = t[l]; }*/


            //int minll;
            newt = INF;
            auto consider = [&](int l, Point*a, Point*b, Point*c){
                //if(a == NIL || b == NIL || c == NIL) return;
                auto tmp = make_event(a, b, c);
                if(tmp > oldt && tmp < newt){
                    newt = tmp;
                    minl = l;
                }
            };
            consider(0, B[i]->prev, B[i], B[i]->next);
            consider(1, o.B[j]->prev, o.B[j], o.B[j]->next);
            consider(2, u->prev, u, v);
            consider(3, u, v, v->next);
            consider(4, u, u->next, v);
            consider(5, u, v->prev, v);
            if (newt == INF) break;

            //for(int l=0;l<6;++l){
            //    cerr << setw(16) << (Frac)t[l];
            //}   cerr << " : " << setw(0) << minl << "\n";
            if (newt == INF) break;
            switch (minl) {
                case 0: if (B[i]->tup() < u->tup()) A[k++] = B[i]; B[i++]->act(); break;
                case 1: if (o.B[j]->tup() > v->tup()) A[k++] = o.B[j]; o.B[j++]->act(); break;
                case 2: A[k++] = u; u = u->prev; break;
                case 3: A[k++] = v; v = v->next; break;
                case 4: A[k++] = u = u->next; break;
                case 5: A[k++] = v = v->prev; break;
            }
        }
        A[k] = NIL;
        u->next = v; v->prev = u; // now go back in time to update pointers
        for (k--; k >= 0; k--) {
            if (A[k]->tup() <= u->tup() || A[k]->tup() >= v->tup()) {
                A[k]->act();
                if (A[k] == u) u = u->prev; else if (A[k] == v) v = v->next;
            } else {
                u->next = A[k];
                A[k]->prev = u;
                v->prev = A[k];
                A[k]->next = v;
                if (A[k]->tup() < mid->tup()) u = A[k];
                else v = A[k];
            }
        }
        //cerr << *u << "  " << *v << "\n";

        //cerr << "did"  << "\n";
        o.hulls.clear();
        hulls.clear();
        // free up memory
        ///cerr << pool_back << " -> " << old_pool_back << "\n";
        pool_back = old_pool_back;
        build_hull();
        o.B.clear();
        o.A.clear();
        o.points.clear();
        o.P.clear();
        P.swap(pres);
        o.list  = 0;
        o.n = 0;
        //cerr << "done"  << "\n";

    }
};
ll Chull3D::inf = 1e17;
Chull3D::Point Chull3D::nil =  Chull3D::Point(inf, inf, inf, 0, 0);
ll operator*(array<ll, 3> const&a, array<ll, 3> const&b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
template<typename T>
void pvec(vector<T> const&v, string const&s){
    string c = (s + string(5, ' ')).substr(0, 5);
    cerr << c;
    for(auto const&e:v) cerr << setw(5) << e << " ";
    cerr << setw(0) << "\n";
}


void test(){
    int n;
    cin >> n;
    vector<array<ll, 3> > v(n);
    for(auto &e:v){
        for(int i=0; i<3; ++i){
            cin >> e[i];
        }
    }
    Chull3D h(v);
}

constexpr ll inf = 1e18;
void solve(){
    debug << fixed << setprecision(10);
    //return test();

    mt19937 rng(53252);
    auto get_rand = [&](int64_t l, int64_t r){
        return uniform_int_distribution<int64_t>(l, r)(rng);
    }; (void) get_rand;
    /// SOLVE HERE
    int n;
    //n = 2e5;
    //n = 1e5;
    cin >> n;
    vector<ll> A(n), B(n);
    for(auto &e:A) cin >> e;
    //for(auto &e:A) e = 1e12*get_rand(1, 1);
    //for(int i=0; i<n; i+=2) --A[i];
    for(auto &e:B) cin >> e;
    //for(auto &e:B) e = get_rand(1, 100);

    for(int i=0; i<n; ++i){
        auto e = A[i];
        A.push_back(e);
    }
    {
        int x = max_element(A.begin(), A.end()) - A.begin();
        rotate(A.begin(), A.begin()+x, A.end());
        rotate(B.begin(), B.begin()+x, B.end());
    }
    A.push_back(A[0]);
    B.push_back(B[0]);

    vector<ll> pre_b(n+2, 0);
    vector<ll> pre_bk(n+2, 0);
    for(int i=0; i<=n; ++i){
        pre_b[i+1] = pre_b[i] + B[i];
        pre_bk[i+1] = pre_bk[i] + i*(ll)B[i];
    }

    vector<int> s(1);
    vector<ll> dp(n+1);
    dp[0] = 0;

    /**const ll cost = 2*(i*pre_bk[j] + j*pre_bk[i] + i*(pre_b[j]*j) - j*(pre_b[i]*i));
    const ll cost = 2*(-i*(pre_bk[j] - pre_b[j]*j) + j*(pre_bk[i] - pre_b[i]*i));
    const ll gain = (A[i] + A[j]) * (ll)(i-j);
    const ll gain = i*A[j] - j*A[i];

    const ll gain = dp[j] + i*(2*(pre_bk[j] - j*pre_b[j]) + A[j]) - (2*(pre_bk[i] - pre_b[i]*i) + A[i]) * j;*/

    vector<array<ll, 3>> tmp; tmp.reserve(n);
    vector<Chull3D> hulls;
    const int l = 1+n/1000;
    for(int i=1;i<=n;++i){
        if((int)tmp.size()>=l){
            hulls.emplace_back(tmp);
            tmp.clear();
            while(hulls.size()>1 && hulls.rbegin()[1].size()<=hulls.rbegin()[0].size()){
                hulls.rbegin()[1].merge(hulls.rbegin()[0]);
                hulls.pop_back();
            }
        }
        array<ll, 3> pq{- (2*(pre_bk[i] - pre_b[i]*i) + A[i]), -i, -1};
        ll val = i*A[0], tmp2;
        for(auto const&e:tmp){
            if(val < (tmp2 = e*pq)){
                val = tmp2;
            }
        }
        for(auto const&e:hulls){
            tmp2 = e.extreme_point(pq)*pq;
            if(val < tmp2){
                val = tmp2;
            }
        }
        dp[i] = val;
        tmp.push_back({i, -(2*(pre_bk[i] - i*pre_b[i]) + A[i]), -dp[i]});
    }

    dp[n] -= 2*n*pre_bk[n];
    dp[n] += n*A[n];
    for(int i=0; i<n; ++i){
        dp[n] += 2*i*(ll)i*B[i];
    }
    //debug << named(dp) << "\n";
    fl ans = dp[n] * 0.5 / n;


    /*{
        vector<int> opt(n+1);
        vector<ll> dp(n+1);
        dp[0] = 0;
        for(int i=1; i<=n; ++i){
            ll me = -inf;
            for(int j=0; j<i; ++j){
                const ll cost = 2*((j+i)*(pre_bk[i]-pre_bk[j]) - i*(ll)j*(pre_b[i]-pre_b[j]));
                //debug << j << " " << i << " : " << cost << "\n";
                const ll gain = (A[i] + A[j]) * (i-j);
                const ll cand = dp[j]+gain-cost;
                if(cand > me){
                    me = cand;
                    opt[i] = j;
                }
            }
            dp[i] = me;
        }
        for(int i=0; i<n; ++i){
            dp[n] += 2*i*(ll)i*B[i];
        }
        fl ans_ = dp[n] * 0.5 / n;
        //debug << named(opt) << "\n";
        debug << " ==> " << ans_ << "\n";
        assert(abs(ans - ans_)/(abs(ans)+abs(ans_)) < 1e-8);
        //debug << named(dp) << "\n";
    }*/

    cout << fixed << setprecision(20) << ans << "\n";
    cerr << "usage: " << pool_back << "\n";

}
