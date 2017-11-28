// convex hull trick in 3D
// requires distinct x-coordinates
// updates need monotone x-coordinates
// O(log(n)^2) per query or update

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

struct Perstree{
    static mt19937 RNG;
    struct Persnode{
        array<shared_ptr<Persnode>, 2> ch;
        array<ll, 3> d;
        ll y;
        Persnode():ch{0, 0}, d{0,0,0}, y(RNG()){}
        Persnode(array<ll, 3> const&v):ch{0, 0}, d(v), y(RNG()){}
        Persnode(shared_ptr<Persnode> const&n):ch(n->ch), d(n->d), y(n->y){}
    };
    static  pair<shared_ptr<Persnode>, shared_ptr<Persnode> > split(shared_ptr<Persnode> const&t, ll const&x){
        if(!t) return pair<shared_ptr<Persnode>, shared_ptr<Persnode> > (0, 0);
        shared_ptr<Persnode> ret = make_shared<Persnode>(t);
        if(x<=t->d[0]){
            auto tmp = split(t->ch[0], x);
            ret->ch[0] = tmp.second;
            return {tmp.first, ret};
        } else {
            auto tmp = split(t->ch[1], x);
            ret->ch[1] = tmp.first;
            return{ret, tmp.second};
        }
    }
    static shared_ptr<Persnode> merge(shared_ptr<Persnode> const&t, shared_ptr<Persnode> const&o){
        if(!t) return o;
        if(!o) return t;
        if(t->y>=o->y){
            shared_ptr<Persnode> ret = make_shared<Persnode>(t);
            ret->ch[1] = merge(ret->ch[1], o);
            return ret;
        } else {
            shared_ptr<Persnode> ret = make_shared<Persnode>(o);
            ret->ch[0] = merge(t, ret->ch[0]);
            return ret;
        }
    }
    shared_ptr<Persnode> root;
    Perstree():root(0){}
    Perstree(shared_ptr<Persnode> r):root(r){}
    Perstree insert_or_delete(array<ll, 3> const&p){
        auto tmp = split(root, p[0]);
        auto tmp2 = split(tmp.second, p[0]+1);
        if(!tmp2.first) return Perstree(merge(merge(tmp.first, make_shared<Persnode>(p)),tmp2.second));
        return Perstree(merge(tmp.first, tmp2.second));
    }
    array<ll, 3> mascal(array<ll, 3> const&dir)const{
        assert(root);
        array<ll, 3> res = root->d;
        stack<shared_ptr<Persnode> > s;
        s.emplace(root);
        while(!s.empty()){
            auto e = s.top();
            s.pop();
            if(!e) continue;
            s.emplace(e->ch[0]);
            s.emplace(e->ch[1]);
            if((e->d[0]-res[0])*dir[0] + (e->d[1]-res[1])*dir[1] + (e->d[2]-res[2])*dir[2]>=0){
                res = e->d;
            }
        }
        return res;
    }
    void print()const{
        stack<pair<shared_ptr<Persnode>, int> > s;
        s.emplace(root, 0);
        while(!s.empty()){
            auto e = s.top();
            s.pop();
            if(e.first==0)continue;
            if(e.second){
                cerr << " " << e.first->d[0] << "/" << e.first->d[1] << "/" << e.first->d[2] << " ";
            } else {
                s.emplace(e.first->ch[0], 0);
                s.emplace(e.first, 1);
                s.emplace(e.first->ch[1], 0);
            }
        }
        cerr << "\n";
    }
};
mt19937 Perstree::RNG;


struct Chull3D{
    using ll = long long;

    struct Frac{
        ll a, b;
        Frac(ll const&x=0, ll const&y=1):a(x), b(y){if(b<0){a=-a; b=-b;}}
        Frac& reduce(){
            ll g = __gcd(a, b);
            a/=g, b/=g;
            if(a<0){ a=-a, b=-b; }
            return *this;
        }
        Frac& operator+(Frac const&o)const{
            ll g = __gcd(b, o.b);
            return Frac(o.b/g*a + b/g*o.a, b/g*o.b).reduce();
        }
        Frac& operator-(Frac const&o)const{
            ll g = __gcd(b, o.b);
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
            return a/(long double)b < o.a/(long double)o.b;
        }
        friend ostream& operator<<(ostream&o, Frac const&f){
            return o << (f.a/(long double)f.b);
        }
    };
    struct Point{
        ll x, y, z;
        Point *prev, *next;
        Point():prev(0), next(0){}
        Point(ll _x, ll _y, ll _z, Point*p, Point*n):x(_x), y(_y), z(_z), prev(p), next(n){}
        void act() {
            if (prev->next != this) prev->next = next->prev = this; // insert
            else { prev->next = next; next->prev = prev; } // delete
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
    ll inf = 1e17;
    Frac INF = Frac(1, 0);
    Point nil = Point(inf, inf, inf, 0, 0);
    Point *NIL = &nil;

    ll turn0(Point *p, Point *q, Point *r) { // <0 iff cw
        if (p == NIL || q == NIL || r == NIL) return 1.0;
        return (q->x-p->x)*(r->y-p->y) - (r->x-p->x)*(q->y-p->y);
    }
    ll turn(Point *p, Point *q, Point *r) { // <0 iff cw
        ll ans = turn0(p, q, r);
        if(ans == 0){
            ans = ((q->x-p->x)*(r->z-p->z) - (r->x-p->x)*(q->z-p->z));
            cerr << "tie break " << ans << "\n";
        }
        return ans;
    }
    Frac time(Point *p, Point *q, Point *r) { // when turn changes
        if (p == NIL || q == NIL || r == NIL) return INF;
        return Frac((q->x-p->x)*(r->z-p->z) - (r->x-p->x)*(q->z-p->z) , turn(p,q,r));
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
        pair<Frac, int> t[6], oldt, newt;
        int i, j, k, l, minl;
        if (n == 1) { A[0] = list->prev = list->next = NIL; return; }

        for (u = list, i = 0; i < n/2-1; u = u->next, i++) ;
        mid = v = u->next;
        hull(list, n/2, B, A); // recurse on left and right sides
        hull(mid, n-n/2, B+n/2*2, A+n/2*2);

        for (;;) // find initial bridge
            if (turn(u, v, v->next) <= 0) v = v->next;
            else if (turn(u->prev, u, v) <= 0) u = u->prev;
            else break;

        //cerr << "::" << n << "\n";
        //cerr << *u << "  " << *v << "\n";
        // merge by tracking bridge uv over time
        int moveu = 0, movev = 0;
        for (i = k = 0, j = n/2*2, oldt = {-INF, -1}; ; oldt = newt) {
            // second doesn't matter apparently, just varies degenerate faces????
            t[0] = {time(B[i]->prev, B[i], B[i]->next), 0};
            t[1] = {time(B[j]->prev, B[j], B[j]->next), 1};
            t[2] = {time(u->prev, u, v), 2};
            t[3] = {time(u, v, v->next), 2};
            t[4] = {time(u, u->next, v), 3};
            t[5] = {time(u, v->prev, v), 3};
            if(moveu == 1) t[2].first = INF;
            if(movev == 1) t[3].first = INF;
            if(moveu == -1) t[4].first = INF;
            if(movev == -1) t[5].first = INF;
            for(int l=0;l<6;++l){
                if(t[l].first.a == 0 && t[l].first.b == 0)
                    t[l].first = oldt.first;//, cerr << "failsafe\n";
            }
            for (newt = {INF, -1}, l = 0; l < 6; l++)
            if (t[l] >= oldt && t[l] < newt) { minl = l; newt = t[l]; }
            //for(int l=0;l<6;++l){
            //    cerr << setw(14) << t[l].first;
            //}   cerr << " : " << setw(0) << minl << "\n";
            if (newt.first == INF) break;
            if(oldt.first < newt.first){ moveu = movev = 0; }
            if(minl == 2){ moveu = -1;}
            if(minl == 3){ movev = -1;}
            if(minl == 4){ moveu = 1;}
            if(minl == 5){ movev = 1;}
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
    }
    size_t n;
    vector<array<ll, 3> > points;
    vector<Point> P;
    vector<Point*> A, B;
    vector<pair<Frac, Perstree>> hulls;
    Chull3D(vector<array<ll, 3> > const&p):n(p.size()), points(p){
        std::sort(points.begin(), points.end());
        for(auto const&e:points)
            P.emplace_back(e[0], e[1], e[2], (Point*)0, (Point*)0);
        Point* list = sort(P.data(), n);
        A.resize(2*n, 0); B.resize(2*n, 0);
        hull(list, n, A.data(), B.data());
        build_hull(list);
    }
    void build_hull(Point*list){
        Perstree hul;
        for(Point*l = list;l!=NIL; l = l->next){
            hul = hul.insert_or_delete(l->arr());
        }
        Frac oldt = -INF;
        for(int i=0; A[i]!=NIL; A[i++]->act()){
            hulls.emplace_back(oldt, hul);
            hul = hul.insert_or_delete(A[i]->arr());
            oldt = time(A[i]->prev, A[i], A[i]->next);
            //cerr << A[i]->arr()[0] << "/" << A[i]->arr()[1] << "/" << A[i]->arr()[2] << "\n";
        }
        hulls.emplace_back(oldt, hul);
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
};

ll operator*(array<ll, 3> const&a, array<ll, 3> const&b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
const ll inf = 1e16;
template<typename T>
void pvec(vector<T> const&v, string const&s){
    string c = (s + string(5, ' ')).substr(0, 5);
    cerr << c;
    for(auto const&e:v) cerr << setw(5) << e << " ";
    cerr << setw(0) << "\n";
}

// DP
vector<ll> DP2(vector<ll> const&v){
    int n = v.size();
    vector<ll> a(n);
    partial_sum(v.begin(), v.end(), a.begin());
    vector<ll> DP(n+1, 0);
    vector<array<ll, 3>> tmp; tmp.reserve(n);
    vector<Chull3D> hulls;
    int l = 1+n/2600;
    for(int i=0;i<n;++i){
        if((int)tmp.size()>=l){
            while(hulls.size() && hulls.back().size()<=tmp.size()){
                copy(hulls.back().points.begin(), hulls.back().points.end(), back_inserter(tmp));
                hulls.pop_back();
            }
            hulls.emplace_back(tmp);
            tmp.clear();
        }
        array<ll, 3> pq{-a[i], -i, -1};
        ll val = a[i], tmp2;
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
        val+=a[i]*i;
        DP[i+1] = val;
        tmp.push_back({i, a[i], -(DP[i+1]+a[i]*i)});
    }
    //pvec(opt, "opt2:");
    //pvec(DP, "DP2:");
    return DP;

}

signed main(){
    #ifdef LOCAL_RUN
    freopen("in.txt", "r", stdin);
    #endif // LOCAL_RUN
    cin.tie(0); ios_base::sync_with_stdio(false);
    int TT; cin >> TT;
    while(TT--){
        int n;cin >> n;
        vector<ll> a(n);
        for(auto &e:a) cin >> e;
        auto D2 = DP2(a);
        cout << D2.back() << "\n";
    }
    return 0;
}
